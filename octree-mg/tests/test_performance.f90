#include "../src/cpp_macros.h"
program test_performance
  use mpi
#ifndef SINGLE_MODULE
  use m_octree_mg
#elif NDIM == 1
  use m_octree_mg_1d
#elif NDIM == 2
  use m_octree_mg_2d
#elif NDIM == 3
  use m_octree_mg_3d
#endif

  implicit none

  integer             :: box_size
  integer             :: domain_size(NDIM)
  real(dp)            :: dr(NDIM), r_min(NDIM) = 0.0_dp
  logical             :: periodic(NDIM)  = .false.
  logical             :: fmg_cycle       = .true.
  integer             :: n_finer         = 0
  real(dp), parameter :: pi              = acos(-1.0_dp)
  character(len=40)   :: arg_string
  integer             :: n, ierr, n_args, n_its
  real(dp)            :: t0, t1, n_unknowns
  type(mg_t)          :: mg

  n_args = command_argument_count()
  if (n_args < NDIM+2 .or. n_args > NDIM+3) then
     error stop "Usage: ./test_uniform_grid box_size domain_size(NDIM) n_its [FMG?]"
  end if

  call get_command_argument(1, arg_string)
  read(arg_string, *) box_size

  do n = 1, NDIM
     call get_command_argument(1+n, arg_string)
     read(arg_string, *) domain_size(n)
  end do

  call get_command_argument(NDIM+2, arg_string)
  read(arg_string, *) n_its

  if (n_args > NDIM+2) then
     call get_command_argument(NDIM+3, arg_string)
     read(arg_string, *) fmg_cycle
  end if

  dr =  1.0_dp / domain_size

  mg%geometry_type = mg_cartesian
  mg%operator_type = mg_laplacian
  mg%smoother_type = mg_smoother_gs

  mg%bc(:, mg_iphi)%bc_type = mg_bc_dirichlet
  mg%bc(:, mg_iphi)%bc_value = 0.0_dp

  call mg_set_methods(mg)
  call mg_comm_init(mg)
  call mg_build_rectangle(mg, domain_size, box_size, dr, r_min, &
       periodic, n_finer)
  call mg_load_balance(mg)
  call mg_load_balance_parents(mg)

  call mg_allocate_storage(mg)
  call set_rhs(mg)

  t0 = mpi_wtime()
  do n = 1, n_its
     if (fmg_cycle) then
        call mg_fas_fmg(mg, n > 1)
     else
        call mg_fas_vcycle(mg)
     end if
  end do
  t1 = mpi_wtime()

  if (mg%my_rank == 0) then
     if (fmg_cycle) then
        print *, "cycle type        FMG"
     else
        print *, "cycle type        V-cycle"
     end if
     print *, "n_cpu            ", mg%n_cpu
     print *, "problem_size     ", domain_size
     print *, "box_size         ", box_size
     print *, "n_iterations     ", n_its
     print *, "time/iteration   ", (t1-t0) / n_its
     print *, "total_time(s)    ", (t1-t0)
     n_unknowns = real(mg%n_boxes, dp) * mg%box_size**NDIM
     print *, "unknowns/microsec", 1e-6_dp * n_its * &
          n_unknowns / (t1-t0)
     print *, ""
  end if
  call mg_timers_show(mg)

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

contains

  subroutine set_rhs(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc, IJK

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
             mg%boxes(id)%cc(IJK, mg_irhs) = 1.0_dp
          end do; CLOSE_DO
       end do
    end do
  end subroutine set_rhs

end program test_performance
