!> Module for cfc
!< Note that this metric solver support with mod_grhd only!
module mod_cfc
  use mod_cfc_parameters
  use mod_cfc_alp
  use mod_cfc_beta
  use mod_cfc_psi
  implicit none
  public

  ! public methods
  public :: cfc_metric_init
  public :: cfc_solve_metric
  public :: cfc_solver_activate
  public :: cfc_update_flag
  public :: cfc_metric_interpolation

  contains

  logical function cfc_update_flag()
    use mod_global_parameters
    cfc_update_flag = .False.

    if (use_cfc) then

       if ( it >= cfc_it_last_update + cfc_dit_update ) then
          ! note that metric will not be updated if the time difference is smaller than cfc_smallest_dt.
          if ( global_time >= cfc_t_last_update + cfc_smallest_dt ) then
             cfc_update_flag=.True.

          end if
       end if

       if ( global_time >= cfc_t_last_update + cfc_dt_update ) then
          cfc_update_flag=.True.
        if (mype == 0) write(*,*) global_time, "u passed cfc_dt update"
       end if

       if (cfc_update_flag) then
          cfc_t_last_update = global_time
          cfc_it_last_update = it
       end if

    end if
  end function cfc_update_flag

  subroutine cfc_solver_activate()
    use mod_global_parameters
    use mod_multigrid_coupling
    integer                      :: n

    namelist /cfc_solver_list/ cfc_tol, cfc_it_max, &
             cfc_print, cfc_psi_tol_init, &
             cfc_smallest_dt, cfc_dt_update, cfc_dit_update, &
             cfc_n_interpolation, reconstruct_cfc, &
             cfc_n_cycle, cfc_redblack

    do n = 1, size(par_files)
       ! Try to read in the namelists. They can be absent or in a different
       ! order, since we rewind before each read.
       rewind(unitpar)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, cfc_solver_list, end=111)
111    close(unitpar)
    end do
    
    if ( cfc_dt_update < 0.0d0 ) then
       call mpistop(" cfc_dt_update can not be negative")
    end if
    if ( cfc_n_interpolation < 2 ) then
       call mpistop(" cfc_n_interpolation < 2")
    end if

    if ( mod(cfc_n_interpolation,2) == 0 ) then
       n = cfc_n_interpolation / 2
    else
       n = (cfc_n_interpolation-1) / 2
    endif
    ! fixme: this seems no needed, as the code not yet read nghostcells
    nghostcells = max(nghostcells, n)
!    if ( n > nghostcells) then
!       write(*,*) "cfc_n_interpolation = ", cfc_n_interpolation
!       write(*,*) "number of ghostcells needed = ", n
!       write(*,*) "nghostcells = ", nghostcells
!       call mpistop(" cfc_n_interpolation is not match with nghostcells")
!    end if

    cfc_t_last_update = global_time
    cfc_it_last_update = it
    
    ! activate the solver
    use_multigrid = .True.
    use_cfc = .True.
  end subroutine cfc_solver_activate

  subroutine cfc_metric_init()
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics
    use mod_ghostcells_update, only: getbc

    integer                        :: id, iigrid, igrid
    double precision               :: psi_err = huge(0.0d0)
    double precision               :: psi_err_old = huge(0.0d0)
    double precision               :: local_err = 0.0d0
    integer, parameter             :: mg_it_max = 10000
    integer                        :: nc, lvl, mg_it
    type(tree_node), pointer       :: pnode

    double precision, allocatable  :: Aij(:^D&,:,:,:)
    double precision, allocatable  :: A2(:^D&,:)

    double precision               :: psi_xns_err, psi_xns_err_local
    double precision               :: alp_xns_err, alp_xns_err_local

    allocate(Aij(ixG^T,1:3,1:3,1:igridstail))
    allocate(A2(ixG^T,1:igridstail))
    Aij = 0.0d0
    A2 = 0.0d0

    ! initialize the solver
    if (cfc_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if

    mg%n_cycle_up   = cfc_n_cycle(1)
    mg%n_cycle_down = cfc_n_cycle(2)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! save the hydro variables
       !pso(igrid)%prim(ixG^T,nhydro_lo:nhydro_hi) = ps(igrid)%prim(ixG^T,nhydro_lo:nhydro_hi)
       pso(igrid)%prim(ixG^T,1:nprim) = ps(igrid)%prim(ixG^T,1:nprim)
       ! make sure all psi and alp are physical
       ps(igrid)%prim(ixG^T, psi_) = max( ps(igrid)%prim(ixG^T, psi_), 1.0d0 )
       ps(igrid)%prim(ixG^T, alp_) = min( max( ps(igrid)%prim(ixG^T, alp_), 0.0d0 ), 1.0d0)
    end do
    !$OMP END PARALLEL DO

    if (mype==0) then
      write(*, '(A,ES9.2,A)') ' Start initializating metric'
      write(*, '(A4,A10,A12,A12,A12)') '  #', 'mg_it', 'psi_err'
    end if

    do mg_it = 1, mg_it_max

       psi_err_old = psi_err
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid)
          ! copy psi
          ps1(igrid)%prim(ixM^T, psi_) = ps(igrid)%prim(ixM^T, psi_)
          ! get the new conserved variables with current psi
          call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),&
                    ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))
       end do
       !$OMP END PARALLEL DO
        ! is not solving beta,  is solving vecX:
       call cfc_solve_beta(.False.,Aij(ixG^T,1:3,1:3,1:igridstail))
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid = 1, igridstail;  igrid = igrids(iigrid);
          ! Geometry subroutines expect this to be set
          block => ps(igrid)
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call cfc_get_Aij_grid(ps(igrid), ixG^LL, ixM^LL, &
                 A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
       end do
       !$OMP END PARALLEL DO
       call cfc_solve_psi(A2(ixG^T,1:igridstail))

       local_err = 0.0d0
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid)
          local_err = max( local_err, maxval( &
                    dabs( (ps(igrid)%prim(ixM^T, psi_) - ps1(igrid)%prim(ixM^T, psi_) ) &
                          /ps1(igrid)%prim(ixM^T, psi_) ) ) )
       end do
       !$OMP END PARALLEL DO
       call mpi_allreduce(local_err, psi_err, 1, mpi_double, mpi_max, icomm, ierrmpi)

       if ( (psi_err <= cfc_psi_tol_init)  ) then
          if (mype==0) then
             write(*,*) "! psi initialization completed! "
             write(*,*) "! mg_it = ", mg_it, "psi_err = ", psi_err
             write(*,*) " "
             print*,'-------------------------------------------------------------------------------'
             write(*,*) " "
          end if
          exit
       else if ( psi_err >= psi_err_old ) then
          if (mype==0) then
             write(*,*) "! Warning: metric is not converging."
             write(*,*) "! mg_it = ", mg_it, "psi_err = ", psi_err
             write(*,*) "! Stop initializating here."
          end if
          exit
       else if (mg_it >= mg_it_max) then
          if (mype==0) then
             write(*,*) "! Warning: metric failed to converge. mg_it = ", mg_it, "psi_err = ", psi_err
             write(*,*) "! Stop initializating here."
          end if
          exit
       end if

       if (mype==0) then
          write(*, '(A4,I10,ES12.3,ES12.3,ES12.3)') " #", &
                   mg_it, psi_err
       end if

    end do

    ! solve rest of the metric variables
    call cfc_solve_alp(A2(ixG^T,1:igridstail))
    call cfc_solve_beta(.True., Aij(ixG^T,1:3,1:3,1:igridstail))

    ! once all the metric variables are solved, Aij is not needed anymore
    deallocate(Aij)
    deallocate(A2)

    psi_xns_err = 0.0d0
    psi_xns_err_local = 0.0d0
    alp_xns_err = 0.0d0
    alp_xns_err_local = 0.0d0
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid)
       psi_xns_err_local = maxval(dabs((ps(igrid)%prim(ixM^T, psi_) - pso(igrid)%prim(ixM^T, psi_) ) &
                       /pso(igrid)%prim(ixM^T, psi_) ))
       alp_xns_err_local = maxval(dabs((ps(igrid)%prim(ixM^T, alp_) - pso(igrid)%prim(ixM^T, alp_) ) &
                       /pso(igrid)%prim(ixM^T, alp_) ))
    end do
    !$OMP END PARALLEL DO

    call mpi_allreduce(psi_xns_err_local, psi_xns_err, 1, mpi_double, mpi_max, icomm, ierrmpi)
    call mpi_allreduce(alp_xns_err_local, alp_xns_err, 1, mpi_double, mpi_max, icomm, ierrmpi)


    if (mype==0) then
       write(*,*) "! metric initialization completed! "
       write(*,*) " "
       write(*,*) " psi_xns_err = ", psi_xns_err
       write(*,*) " alp_xns_err = ", alp_xns_err
       write(*,*) " "
       print*,'-------------------------------------------------------------------------------'
       write(*,*) " "
    end if

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! restore the hydro variables
       ps(igrid)%prim(ixM^T,nhydro_lo:nhydro_hi) = pso(igrid)%prim(ixM^T,nhydro_lo:nhydro_hi)
    end do
    !$OMP END PARALLEL DO

    ! update the boundaries
    call getbc(global_time,0.0d0,ps,1,nprim,.false.)

  end subroutine cfc_metric_init

  subroutine cfc_solve_metric()
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics

    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl
    type(tree_node), pointer       :: pnode

    double precision, allocatable  :: Aij(:^D&,:,:,:)
    double precision, allocatable  :: A2(:^D&,:)

    allocate(Aij(ixG^T,1:3,1:3,1:igridstail))
    allocate(A2(ixG^T,1:igridstail))
    Aij = 0.0d0
    A2 = 0.0d0

    ! initialize the solver
    if (cfc_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if

    mg%n_cycle_up   = cfc_n_cycle(1)
    mg%n_cycle_down = cfc_n_cycle(2)

    !-----------------------------------------------------------------------
    ! Step 1: Solve the PDE for X and thus the conformal extrinsic curvature
    !-----------------------------------------------------------------------
    call cfc_solve_beta(.False.,Aij(ixG^T,1:3,1:3,1:igridstail))

    ! obtain the conformal extrinsic curvature
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);
       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call cfc_get_Aij_grid(ps(igrid), ixG^LL, ixM^LL, &
                 A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    ! Step 2: After having the extrinsic curvature, we can solve psi
    !-----------------------------------------------------------------------
    call cfc_solve_psi(A2(ixG^T,1:igridstail))

    !-----------------------------------------------------------------------
    ! Step 3: Before solving the alp, we need to know S_star first 
    !             (That means update the prim variables)
    !-----------------------------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call phys_to_primitive(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),&
                 ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    ! Step 4: Solve alp
    !-----------------------------------------------------------------------
    call cfc_solve_alp(A2(ixG^T,1:igridstail))

    !-----------------------------------------------------------------------
    ! Step 5: Solve shift vector beta
    !-----------------------------------------------------------------------
    call cfc_solve_beta(.True., Aij(ixG^T,1:3,1:3,1:igridstail))

    ! once all the metric variables are solved, Aij is not needed anymore
    deallocate(Aij)
    deallocate(A2)
  end subroutine cfc_solve_metric

  !> calculate Aij on a single grid
  subroutine cfc_get_Aij_grid( ps_in, ixI^L, ixO^L, A2, Aij)
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: vecX
    type(state), intent(in)         :: ps_in
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out)   :: A2(ixI^S)
    double precision, intent(out)   :: Aij(ixI^S,1:3,1:3)

    integer                         :: jxO^L, hxO^L
    integer                         :: idir, jdir, kdir, ldir
    double precision                :: fij(ixI^S, 1:3, 1:3)
    double precision                :: dvecX(ixI^S,1:3,1:3)
    {^NOONED
    double precision                :: cot_theta(ixI^S)
    }

    associate(prim=>ps_in%prim, x=>ps_in%x)

    call get_gamma_ij_hat(x(ixI^S, 1:^ND),ixI^L,ixO^L,fij(ixI^S, 1:3, 1:3))

    dvecX = 0.0d0
    do idir = 1, ndim
       do jdir = 1, ndir
          call partial_d( prim(ixI^S,vecX(jdir)), ixI^L, ixO^L, idir, dvecX(ixI^S, jdir, idir))
       end do
    end do

    Aij = 0.0d0
    select case (coordinate)
    case (cartesian)
       Aij(ixO^S, 1,1) = 2.0d0/3.0d0 * ( 2.0d0 * dvecX(ixO^S,1,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } & 
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 2,2) = 2.0d0/3.0d0 * ( - dvecX(ixO^S,1,1) &
                         {^NOONED + 2.0d0 * dvecX(ixO^S,2,2) } &
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 3,3) = 2.0d0/3.0d0 * ( - dvecX(ixO^S,1,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } &
                         {^IFTHREED + 2.0d0 * dvecX(ixO^S,3,3) } )
    case (cylindrical)
       Aij(ixO^S, 1,1) = 2.0d0/3.0d0 * ( 2.0d0 * dvecX(ixO^S,1,1) - prim(ixO^S,vecX(1))/x(ixO^S,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } & 
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 2,2) = 2.0d0/3.0d0 * ( -dvecX(ixO^S,1,1) - prim(ixO^S,vecX(1))/x(ixO^S,1) &
                         {^NOONED + 2.0d0 * dvecX(ixO^S,2,2) } &
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 3,3) = 2.0d0/3.0d0/fij(ixO^S,3,3) * ( -dvecX(ixO^S,1,1) + 2.0d0 * prim(ixO^S,vecX(1))/x(ixO^S,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } &
                         {^IFTHREED + 2.0d0 * dvecX(ixO^S,3,3) } )
    case (spherical)
       {^NOONED   cot_theta(ixO^S) = dcos(x(ixO^S,2))/dsin(x(ixO^S,2))    }
       Aij(ixO^S, 1,1) = 2.0d0/3.0d0 * ( 2.0d0 * ( dvecX(ixO^S,1,1) - prim(ixO^S,vecX(1))/x(ixO^S,1) )&
                         {^NOONED - dvecX(ixO^S,2,2) - cot_theta(ixO^S) * prim(ixO^S,vecX(2)) } & 
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 2,2) = 2.0d0/3.0d0/fij(ixO^S,2,2) * ( -dvecX(ixO^S,1,1) + prim(ixO^S,vecX(1))/x(ixO^S,1) &
                         {^NOONED + 2.0d0 * dvecX(ixO^S,2,2) - cot_theta(ixO^S) * prim(ixO^S,vecX(2)) } &
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 3,3) = 2.0d0/3.0d0/fij(ixO^S,3,3) * ( -dvecX(ixO^S,1,1) + prim(ixO^S,vecX(1))/x(ixO^S,1) &
                         {^NOONED - dvecX(ixO^S,2,2) + 2.0d0 * cot_theta(ixO^S) * prim(ixO^S,vecX(2)) } &
                         {^IFTHREED + 2.0d0 * dvecX(ixO^S,3,3) } )
   
    end select
    {^NOONED
    Aij(ixO^S, 1,2) =  dvecX(ixO^S,1,2)/fij(ixO^S,2,2) + dvecX(ixO^S,2,1) 
    Aij(ixO^S, 1,3) =  dvecX(ixO^S,3,1)                + dvecX(ixO^S,1,3)/fij(ixO^S,3,3)
    Aij(ixO^S, 2,3) =  dvecX(ixO^S,3,2)/fij(ixO^S,2,2) + dvecX(ixO^S,2,3)/fij(ixO^S,3,3)
    do idir = 1, 2
       do jdir = idir+1, 3
          Aij(ixO^S, jdir, idir) = Aij(ixO^S, idir, jdir) 
       end do
    end do
    }

    ! get A2
    A2=0.0d0
    do jdir = 1, 3; do idir = 1, 3
       do kdir = 1, 3
          do ldir = 1, 3
             A2(ixO^S) = A2(ixO^S) &
              + fij(ixO^S,idir,kdir) * fij(ixO^S,jdir,ldir) &
                  * Aij(ixO^S,idir,jdir) * Aij(ixO^S,kdir,ldir) 
          end do
       end do
    end do; end do

    end associate
  end subroutine cfc_get_Aij_grid

  subroutine cfc_metric_interpolation(ixI^L,ixO^L,idims,prim,x,primi,xi)
    use mod_global_parameters
    use mod_physics
    use mod_interpolation

    integer, intent(in) :: ixI^L, ixO^L, idims
    double precision, dimension(ixI^S,1:nprim) :: prim
    double precision, dimension(ixI^S,1:nprim) :: primi
    double precision, dimension(ixI^S,1:ndim) :: x, xi

    ! local vars
    integer :: ix^D
    integer :: iw
    integer :: n_lo, n_hi
    double precision, allocatable :: xx(:), yy(:)

    if ( .not. use_CFC ) return
 
    if ( mod(cfc_n_interpolation,2) == 0 ) then
       n_lo = cfc_n_interpolation / 2 - 1
       n_hi = cfc_n_interpolation / 2
    else
       n_lo = (cfc_n_interpolation+1) / 2 - 1
       n_hi = (cfc_n_interpolation-1) / 2
    endif
    allocate(xx(cfc_n_interpolation))
    allocate(yy(cfc_n_interpolation))

    do iw = nmetric_lo, nmetric_hi

       {^IFONED
       do ix1 = ixOmin1, ixOmax1
          xx = x(ix1-n_lo:ix1+n_hi, 1)
          yy = prim(ix1-n_lo:ix1+n_hi, iw)
          call lagrange_interpolation(xx,yy,&
                  xi(ix1,1),primi(ix1, iw))
       end do
       }
       {^IFTWOD
       select case (idims)
       case (1)
          do ix2 = ixOmin2, ixOmax2
             do ix1 = ixOmin1, ixOmax1
                xx = x(ix1-n_lo:ix1+n_hi, ix2, 1)
                yy = prim(ix1-n_lo:ix1+n_hi, ix2, iw)
                call lagrange_interpolation(xx,yy,&
                        xi(ix1, ix2, 1), primi(ix1, ix2, iw))
             end do
          end do
       case (2)
          do ix1 = ixOmin1, ixOmax1
             do ix2 = ixOmin2, ixOmax2
                xx = x(ix1, ix2-n_lo:ix2+n_hi, 2)
                yy = prim(ix1, ix2-n_lo:ix2+n_hi, iw)
                call lagrange_interpolation(xx,yy,&
                        xi(ix1, ix2, 2), primi(ix1, ix2, iw))
             end do
          end do
       case default
          call mpistop(" idims can only be 1 or 2 in 2D")
       end select
       }
       {^IFTHREED
       select case (idims)
       case (1)
          do ix3 = ixOmin3, ixOmax3
             do ix2 = ixOmin2, ixOmax2
                do ix1 = ixOmin1, ixOmax1
                   xx = x(ix1-n_lo:ix1+n_hi, ix2, ix3, 1)
                   yy = prim(ix1-n_lo:ix1+n_hi, ix2, ix3, iw)
                   call lagrange_interpolation(xx,yy,&
                           xi(ix^D, 1), primi(ix^D, iw))
                end do
             end do
          end do
       case (2)
          do ix3 = ixOmin3, ixOmax3
             do ix1 = ixOmin1, ixOmax1
                do ix2 = ixOmin2, ixOmax2
                   xx = x(ix1, ix2-n_lo:ix2+n_hi, ix3, 2)
                   yy = prim(ix1, ix2-n_lo:ix2+n_hi, ix3, iw)
                   call lagrange_interpolation(xx,yy,&
                           xi(ix^D, 2), primi(ix^D, iw))
                end do
             end do
          end do
       case (3)
          do ix1 = ixOmin1, ixOmax1
             do ix2 = ixOmin2, ixOmax2
                do ix3 = ixOmin3, ixOmax3
                   xx = x(ix1, ix2, ix3-n_lo:ix3+n_hi, 3)
                   yy = prim(ix1, ix2, ix3-n_lo:ix3+n_hi, iw)
                   call lagrange_interpolation(xx,yy,&
                           xi(ix^D, 3), primi(ix^D, iw))
                end do
             end do
          end do
       case default
          call mpistop(" idims can only be 1, 2 or 3 in 3D")
       end select
       }

    end do

    deallocate(xx)
    deallocate(yy)
  end subroutine cfc_metric_interpolation

end module mod_cfc
