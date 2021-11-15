#include "cpp_macros.h"
module m_multigrid
  use m_data_structures
  use m_prolong
  use m_restrict
  use m_ghost_cells

  implicit none
  private

  integer :: timer_total_vcycle  = -1
  integer :: timer_total_fmg     = -1
  integer :: timer_smoother      = -1
  integer :: timer_smoother_gc   = -1
  integer :: timer_coarse        = -1
  integer :: timer_correct       = -1
  integer :: timer_update_coarse = -1

  ! Public methods
  public :: mg_fas_vcycle
  public :: mg_fas_fmg
  public :: mg_set_methods
  public :: mg_apply_op

contains

  subroutine mg_set_methods(mg)
    use m_laplacian
    use m_vlaplacian
    use m_helmholtz
    use m_vhelmholtz
    use m_cfc_psi
    use m_cfc_alp
    use m_cfc_beta
    type(mg_t), intent(inout) :: mg

    ! Set default prolongation method (routines below can override this)
    mg%box_prolong => mg_prolong_sparse

    select case (mg%operator_type)
    case (mg_laplacian)
       call laplacian_set_methods(mg)
    case (mg_vlaplacian)
       call vlaplacian_set_methods(mg)
    case (mg_helmholtz)
       call helmholtz_set_methods(mg)
    case (mg_vhelmholtz)
       call vhelmholtz_set_methods(mg)
    case (mg_cfc_psi)
       call cfc_psi_set_methods(mg)
    case (mg_cfc_alp)
       call cfc_alp_set_methods(mg)
    case (mg_cfc_beta)
       call cfc_beta_set_methods(mg)
    case default
       error stop "mg_set_methods: unknown operator"
    end select

    ! For red-black, perform two smoothing sub-steps so that all unknowns are
    ! updated per cycle
    if (mg%smoother_type == mg_smoother_gsrb) then
       mg%n_smoother_substeps = 2
    else
       mg%n_smoother_substeps = 1
    end if
  end subroutine mg_set_methods

  subroutine check_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (.not. associated(mg%box_op) .or. &
         .not. associated(mg%box_smoother)) then
       call mg_set_methods(mg)
    end if

  end subroutine check_methods

  subroutine mg_add_timers(mg)
    type(mg_t), intent(inout) :: mg
    timer_total_vcycle  = mg_add_timer(mg, "mg total V-cycle")
    timer_total_fmg     = mg_add_timer(mg, "mg total FMG cycle")
    timer_smoother      = mg_add_timer(mg, "mg smoother")
    timer_smoother_gc   = mg_add_timer(mg, "mg smoother g.c.")
    timer_coarse        = mg_add_timer(mg, "mg coarse")
    timer_correct       = mg_add_timer(mg, "mg correct")
    timer_update_coarse = mg_add_timer(mg, "mg update coarse")
  end subroutine mg_add_timers

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid).
  subroutine mg_fas_fmg(mg, have_guess, max_res)
    type(mg_t), intent(inout)       :: mg
    logical, intent(in)             :: have_guess !< If false, start from phi = 0
    real(dp), intent(out), optional :: max_res    !< Store max(abs(residual))
    integer                         :: lvl, i, id, idir

    call check_methods(mg)
    if (timer_smoother == -1) call mg_add_timers(mg)

    call mg_timer_start(mg%timers(timer_total_fmg))

    if (.not. have_guess) then
       do lvl = mg%highest_lvl, mg%lowest_lvl, -1
          do i = 1, size(mg%lvls(lvl)%my_ids)
             id = mg%lvls(lvl)%my_ids(i)
             if (.not.mg%vector_equation) then
                mg%boxes(id)%cc(DTIMES(:), mg_iphi) = 0.0_dp
             else
                mg%boxes(id)%cc(DTIMES(:), mg_vec_iphi(1:3)) = 0.0_dp
             end if
          end do
       end do
    end if

    ! Ensure ghost cells are filled correctly
    if (.not.mg%vector_equation) then
       call mg_fill_ghost_cells_lvl(mg, mg%highest_lvl, mg_iphi)
    else
       do idir =1, 3
          call mg_fill_ghost_cells_lvl(mg, mg%highest_lvl, mg_vec_iphi(idir))
       end do
    end if

    do lvl = mg%highest_lvl,  mg%lowest_lvl+1, -1
       ! Set rhs on coarse grid and restrict phi
       call mg_timer_start(mg%timers(timer_update_coarse))
       call update_coarse(mg, lvl)
       call mg_timer_end(mg%timers(timer_update_coarse))
    end do

    if (mg%subtract_mean) then
       ! For fully periodic solutions, the mean source term has to be zero
       call subtract_mean(mg, mg_irhs, .false.)
    end if

    do lvl = mg%lowest_lvl, mg%highest_lvl
       ! Store phi_old
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          if (.not.mg%vector_equation) then
             mg%boxes(id)%cc(DTIMES(:), mg_iold) = &
               mg%boxes(id)%cc(DTIMES(:), mg_iphi)
          else
             do idir =1, 3
                mg%boxes(id)%cc(DTIMES(:), mg_vec_iold(idir)) = &
                  mg%boxes(id)%cc(DTIMES(:), mg_vec_iphi(idir))
             end do
          end if
       end do

       if (lvl > mg%lowest_lvl) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call mg_timer_start(mg%timers(timer_correct))
          call correct_children(mg, lvl-1)
          call mg_timer_end(mg%timers(timer_correct))

          ! Update ghost cells
          if (.not.mg%vector_equation) then
             call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
          else
             do idir =1, 3
                call mg_fill_ghost_cells_lvl(mg, lvl, mg_vec_iphi(idir))
             end do
          end if
       end if

       ! Perform V-cycle, possibly set residual on last iteration
       if (lvl == mg%highest_lvl) then
          call mg_fas_vcycle(mg, lvl, max_res, standalone=.false.)
       else
          call mg_fas_vcycle(mg, lvl, standalone=.false.)
       end if
    end do

    call mg_timer_end(mg%timers(timer_total_fmg))
  end subroutine mg_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme).
  subroutine mg_fas_vcycle(mg, highest_lvl, max_res, standalone)
    use mpi
    type(mg_t), intent(inout)       :: mg
    integer, intent(in), optional   :: highest_lvl !< Maximum level for V-cycle
    real(dp), intent(out), optional :: max_res     !< Store max(abs(residual))
    !> Whether the V-cycle is called by itself (default: true)
    logical, intent(in), optional   :: standalone
    integer                         :: lvl, min_lvl, i, max_lvl, ierr, idir
    real(dp)                        :: init_res, res
    logical                         :: is_standalone

    is_standalone = .true.
    if (present(standalone)) is_standalone = standalone

    call check_methods(mg)
    if (timer_smoother == -1) call mg_add_timers(mg)

    if (is_standalone) &
         call mg_timer_start(mg%timers(timer_total_vcycle))

    if (mg%subtract_mean .and. .not. present(highest_lvl)) then
       ! Assume that this is a stand-alone call. For fully periodic solutions,
       ! ensure the mean source term is zero.
       call subtract_mean(mg, mg_irhs, .false.)
    end if

    min_lvl = mg%lowest_lvl
    max_lvl = mg%highest_lvl
    if (present(highest_lvl)) max_lvl = highest_lvl

    ! Ensure ghost cells are filled correctly
    if (is_standalone) then
       if (.not.mg%vector_equation) then
          call mg_fill_ghost_cells_lvl(mg, max_lvl, mg_iphi)
       else
          do idir =1, 3
             call mg_fill_ghost_cells_lvl(mg, max_lvl, mg_vec_iphi(idir))
          end do
       end if
    end if

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call smooth_boxes(mg, lvl, mg%n_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy mg_iphi to mg_iold for the
       ! correction later
       call mg_timer_start(mg%timers(timer_update_coarse))
       call update_coarse(mg, lvl)
       call mg_timer_end(mg%timers(timer_update_coarse))
    end do

    call mg_timer_start(mg%timers(timer_coarse))
    if (.not. all(mg%boxes(mg%lvls(min_lvl)%ids)%rank == &
         mg%boxes(mg%lvls(min_lvl)%ids(1))%rank)) then
       error stop "Multiple CPUs for coarse grid (not implemented yet)"
    end if

    init_res = max_residual_lvl(mg, min_lvl)
    do i = 1, mg%max_coarse_cycles
       call smooth_boxes(mg, min_lvl, mg%n_cycle_up+mg%n_cycle_down)
       res = max_residual_lvl(mg, min_lvl)
       if (res < mg%residual_coarse_rel * init_res .or. &
            res < mg%residual_coarse_abs) exit
    end do
    call mg_timer_end(mg%timers(timer_coarse))

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, max_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call mg_timer_start(mg%timers(timer_correct))
       call correct_children(mg, lvl-1)

       ! Have to fill ghost cells after correction
       if (.not.mg%vector_equation) then
          call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
       else
          do idir =1, 3
             call mg_fill_ghost_cells_lvl(mg, lvl, mg_vec_iphi(idir))
          end do
       end if
       call mg_timer_end(mg%timers(timer_correct))

       ! Upwards relaxation
       call smooth_boxes(mg, lvl, mg%n_cycle_up)
    end do

    if (present(max_res)) then
       init_res = 0.0_dp
       do lvl = min_lvl, max_lvl
          res = max_residual_lvl(mg, lvl)
          init_res = max(res, init_res)
       end do
       call mpi_allreduce(init_res, max_res, 1, &
            mpi_double, mpi_max, mg%comm, ierr)
    end if

    ! Subtract mean(phi) from phi
    if (mg%subtract_mean) then
       call subtract_mean(mg, mg_iphi, .true.)
    end if

    if (is_standalone) &
         call mg_timer_end(mg%timers(timer_total_vcycle))
  end subroutine mg_fas_vcycle

  subroutine subtract_mean(mg, iv, include_ghostcells)
    use mpi
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iv
    logical, intent(in)       :: include_ghostcells
    integer                   :: i, id, lvl, nc, ierr
    real(dp)                  :: sum_iv, mean_iv, volume

    nc = mg%box_size
    sum_iv = get_sum(mg, iv)
    call mpi_allreduce(sum_iv, mean_iv, 1, &
         mpi_double, mpi_sum, mg%comm, ierr)

    ! Divide by total grid volume to get mean
    volume = nc**NDIM * product(mg%dr(:, 1)) * size(mg%lvls(1)%ids)
    mean_iv = mean_iv / volume

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          if (include_ghostcells) then
             mg%boxes(id)%cc(DTIMES(:), iv) = &
                  mg%boxes(id)%cc(DTIMES(:), iv) - mean_iv
          else
             mg%boxes(id)%cc(DTIMES(1:nc), iv) = &
                  mg%boxes(id)%cc(DTIMES(1:nc), iv) - mean_iv
          end if
       end do
    end do
  end subroutine subtract_mean

  real(dp) function get_sum(mg, iv)
    type(mg_t), intent(in) :: mg
    integer, intent(in)    :: iv
    integer                :: lvl, i, id, nc
    real(dp)               :: w

    get_sum = 0.0_dp
    do lvl = 1, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       w  = product(mg%dr(:, lvl)) ! Adjust for non-Cartesian cases
       do i = 1, size(mg%lvls(lvl)%my_leaves)
          id = mg%lvls(lvl)%my_leaves(i)
          get_sum = get_sum + w * &
               sum(mg%boxes(id)%cc(DTIMES(1:nc), iv))
       end do
    end do
  end function get_sum

  real(dp) function max_residual_lvl(mg, lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, nc, idir
    real(dp)                  :: res

    nc           = mg%box_size_lvl(lvl)
    max_residual_lvl = 0.0_dp

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call residual_box(mg, id, nc)
       if (.not.mg%vector_equation) then
          res = maxval(abs(mg%boxes(id)%cc(DTIMES(1:nc), mg_ires)))
       else
          res = 0.0_dp
          do idir = 1,3
             res = max(res, &
               maxval(abs(mg%boxes(id)%cc(DTIMES(1:nc), mg_vec_ires(idir)))) )
          end do
       end if
       max_residual_lvl = max(max_residual_lvl, res)
    end do
  end function max_residual_lvl

  !   subroutine solve_coarse_grid(mg)
  !     use m_fishpack
  !     type(mg_t), intent(inout) :: mg

  !     real(dp) :: rhs(DTIMES(mg%box_size))
  !     real(dp) :: rmin(NDIM), rmax(NDIM)
  !     integer  :: nc, nx(NDIM), my_boxes, total_boxes

  !     my_boxes    = size(mg%lvls(1)%my_ids)
  !     total_boxes = size(mg%lvls(1)%ids)
  !     nc          = mg%box_size

  !     if (my_boxes == total_boxes) then
  !        nx(:) = nc
  !        rmin  = [DTIMES(0.0_dp)]
  !        rmax  = mg%dr(1) * [DTIMES(nc)]
  !        rhs   = mg%boxes(1)%cc(DTIMES(1:nc), mg_irhs)

  ! #if NDIM == 2
  !        call fishpack_2d(nx, rhs, mg%bc, rmin, rmax)
  ! #elif NDIM == 3
  !        call fishpack_3d(nx, rhs, mg%bc, rmin, rmax)
  ! #endif

  !        mg%boxes(1)%cc(DTIMES(1:nc), mg_iphi) = rhs
  !     else if (my_boxes > 0) then
  !        error stop "Boxes at level 1 at different processors"
  !     end if

  !     call fill_ghost_cells_lvl(mg, 1)
  !   end subroutine solve_coarse_grid

  ! Set rhs on coarse grid, restrict phi, and copy mg_iphi to mg_iold for the
  ! correction later
  subroutine update_coarse(mg, lvl)
    type(mg_t), intent(inout) :: mg     !< Tree containing full grid
    integer, intent(in)       :: lvl !< Update coarse values at lvl-1
    integer                   :: i, id, nc, nc_c, idir

    nc   = mg%box_size_lvl(lvl)
    nc_c = mg%box_size_lvl(lvl-1)

    ! Compute residual
    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call residual_box(mg, id, nc)
    end do

    if (.not.mg%vector_equation) then
       ! Restrict phi and the residual
       call mg_restrict_lvl(mg, mg_iphi, lvl)
       call mg_restrict_lvl(mg, mg_ires, lvl)
   
       ! Restrict temp variables 1 and 2
       if ((mg%operator_type==mg_cfc_alp) .or.&
           (mg%operator_type==mg_cfc_psi) ) then
          call mg_restrict_lvl(mg, mg_itmp1, lvl)
          call mg_restrict_lvl(mg, mg_itmp2, lvl)
       end if
   
       call mg_fill_ghost_cells_lvl(mg, lvl-1, mg_iphi)
   
       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
       ! store current coarse phi in old.
       do i = 1, size(mg%lvls(lvl-1)%my_parents)
          id = mg%lvls(lvl-1)%my_parents(i)
   
          ! Set rhs = L phi
          call mg%box_op(mg, id, nc_c, mg_irhs)
   
          ! Add the fine grid residual to rhs
          mg%boxes(id)%cc(DTIMES(1:nc_c), mg_irhs) = &
               mg%boxes(id)%cc(DTIMES(1:nc_c), mg_irhs) + &
               mg%boxes(id)%cc(DTIMES(1:nc_c), mg_ires)
   
          ! Story a copy of phi
          mg%boxes(id)%cc(DTIMES(:), mg_iold) = &
               mg%boxes(id)%cc(DTIMES(:), mg_iphi)
       end do

    else
       ! Restrict phi and the residual
       do idir = 1,3
          call mg_restrict_lvl(mg, mg_vec_iphi(idir), lvl)
          call mg_restrict_lvl(mg, mg_vec_ires(idir), lvl)
       end do
   
       do idir = 1,3
          call mg_fill_ghost_cells_lvl(mg, lvl-1, mg_vec_iphi(idir))
       end do
   
       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
       ! store current coarse phi in old.
       do i = 1, size(mg%lvls(lvl-1)%my_parents)
          id = mg%lvls(lvl-1)%my_parents(i)
   
          do idir = 1,3
             ! Set rhs = L phi
             call mg%box_vec_op(mg, id, nc_c, idir, mg_vec_irhs(idir))
          end do
   
          do idir = 1,3
             ! Add the fine grid residual to rhs
             mg%boxes(id)%cc(DTIMES(1:nc_c), mg_vec_irhs(idir)) = &
                  mg%boxes(id)%cc(DTIMES(1:nc_c), mg_vec_irhs(idir)) + &
                  mg%boxes(id)%cc(DTIMES(1:nc_c), mg_vec_ires(idir))
          end do
      
          do idir = 1,3
             ! Story a copy of phi
             mg%boxes(id)%cc(DTIMES(:), mg_vec_iold(idir)) = &
                  mg%boxes(id)%cc(DTIMES(:), mg_vec_iphi(idir))
          end do
   
       end do
    end if
  end subroutine update_coarse

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(mg, lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, idir

    if (.not.mg%vector_equation) then
       do i = 1, size(mg%lvls(lvl)%my_parents)
          id = mg%lvls(lvl)%my_parents(i)
   
          ! Store the correction in mg_ires
          mg%boxes(id)%cc(DTIMES(:), mg_ires) = &
               mg%boxes(id)%cc(DTIMES(:), mg_iphi) - &
               mg%boxes(id)%cc(DTIMES(:), mg_iold)
       end do
       call mg_prolong(mg, lvl, mg_ires, mg_iphi, mg%box_prolong, add=.true.)
    else
       do i = 1, size(mg%lvls(lvl)%my_parents)
          id = mg%lvls(lvl)%my_parents(i)
   
          ! Store the correction in mg_ires
          do idir = 1, 3
             mg%boxes(id)%cc(DTIMES(:), mg_vec_ires(idir)) = &
                  mg%boxes(id)%cc(DTIMES(:), mg_vec_iphi(idir)) - &
                  mg%boxes(id)%cc(DTIMES(:), mg_vec_iold(idir))
          end do
       end do
       do idir = 1, 3
          call mg_prolong(mg, lvl, mg_vec_ires(idir), mg_vec_iphi(idir), mg%box_prolong, add=.true.)
       end do
    end if
  end subroutine correct_children

  subroutine smooth_boxes(mg, lvl, n_cycle)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer, intent(in)       :: n_cycle !< Number of cycles to perform
    integer                   :: n, i, id, nc, idir

    nc = mg%box_size_lvl(lvl)

    do n = 1, n_cycle * mg%n_smoother_substeps
       call mg_timer_start(mg%timers(timer_smoother))
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call mg%box_smoother(mg, id, nc, n)
       end do
       call mg_timer_end(mg%timers(timer_smoother))

       call mg_timer_start(mg%timers(timer_smoother_gc))
       if (.not.mg%vector_equation) then
          call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
       else
          do idir =1, 3
             call mg_fill_ghost_cells_lvl(mg, lvl, mg_vec_iphi(idir))
          end do
       end if
       call mg_timer_end(mg%timers(timer_smoother_gc))
    end do
  end subroutine smooth_boxes

  subroutine residual_box(mg, id, nc)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer                   :: idir 

    if (.not.mg%vector_equation) then
       call mg%box_op(mg, id, nc, mg_ires)
   
       mg%boxes(id)%cc(DTIMES(1:nc), mg_ires) = &
            mg%boxes(id)%cc(DTIMES(1:nc), mg_irhs) &
            - mg%boxes(id)%cc(DTIMES(1:nc), mg_ires)
    else
       do idir = 1, 3
          call mg%box_vec_op(mg, id, nc, idir , mg_vec_ires(idir))
      
          mg%boxes(id)%cc(DTIMES(1:nc), mg_vec_ires(idir)) = &
               mg%boxes(id)%cc(DTIMES(1:nc), mg_vec_irhs(idir)) &
               - mg%boxes(id)%cc(DTIMES(1:nc), mg_vec_ires(idir))
       end do
    end if
  end subroutine residual_box

  !> Apply operator to the tree and store in variable i_out
  subroutine mg_apply_op(mg, i_out, op)
    type(mg_t), intent(inout)      :: mg
    integer, intent(in)            :: i_out
    procedure(mg_box_op), optional :: op
    integer                        :: lvl, i, id, nc

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          if (present(op)) then
             call op(mg, id, nc, i_out)
          else
             call mg%box_op(mg, id, nc, i_out)
          end if
       end do
    end do
  end subroutine mg_apply_op

end module m_multigrid
