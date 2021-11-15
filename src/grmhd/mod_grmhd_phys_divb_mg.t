module mod_grmhd_phys_divb_mg
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_divb_mg_init
  public :: grmhd_clean_divb_multigrid

contains

  !> Initialize the module
  subroutine grmhd_phys_divb_mg_init()
    use mod_global_parameters
    use_multigrid = .true.
  end subroutine grmhd_phys_divb_mg_init

  ! elliptic cleaning on Bvec
  subroutine grmhd_clean_divb_multigrid()
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ix^L, idir
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixG^T), grad(ixG^T, ndir)
    double precision             :: res = huge(1.0d0)
    double precision             :: max_divb

    mg%operator_type = mg_laplacian
    if (divB_mg_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if
    mg%n_cycle_up   = divB_mg_n_cycle(1)
    mg%n_cycle_down = divB_mg_n_cycle(2)
    call mg_set_methods(mg)

    ! Set boundary conditions
    do n = 1, 2*ndim
       idir = (n+1)/2
       select case (typeboundary(Bvec(idir), n))
       case ('symm')
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('special')
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case default
          !print *, "divb_multigrid warning: unknown b.c.: ", &
          !     trim(typeboundary(Bvec(idir), n))
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    if (divB_4thorder) then
       ix^L=ixM^LL^LADD2;
    else
       ix^L=ixM^LL^LADD1;
    end if
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! initial guess of (phi)
       mg%boxes(id)%cc(1:nc^D&, mg_iphi) = 0.0d0

       ! note that the input B field are Bcons
       do idir = 1, ndim
          grad(ix^S,idir) = ps(igrid)%prim(ix^S, psi_)**6 * ps(igrid)%prim(ix^S, Bvec(idir))
       end do
       call grmhd_get_divb(ixG^LL, ixM^LL, grad(ixG^T,1:ndir), tmp(ixG^T), divb_4thorder)
       ! divB as the RHS
       mg%boxes(id)%cc({1:nc}, mg_irhs) = tmp(ixM^T)

       ! fixme:
       max_divb = max(max_divb, maxval(abs(tmp(ixM^T))))
       !write(*,*) max_divb
    end do

    ! fixme: maybe we could do laplician first and see if res < tol
    ! Solve laplacian(phi) = divB
    if(stagger_grid) then
      call mpistop("not ready yet")
    else
      do n = 1, divB_mg_it_max
         !call mg_fas_vcycle(mg, max_res=res)
         call mg_fas_fmg(mg, .True., max_res=res)
         !write(*,*) 'divB cleaning :', n, res
         if (res < divB_mg_tol) exit
      end do
      if (res > divB_mg_tol) then 
         !call mpistop("divb_multigrid: no convergence")
         write(*,*) "Warning: divb_multigrid no convergence, res = ", res
      end if
    end if

    ix^L=ixM^LL^LADD1;

    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! Compute the gradient of phi, we need ghost zones as well
       tmp(ix^S) = mg%boxes(id)%cc({0:nc+1,}, mg_iphi)
       if(stagger_grid) then
         call mpistop("not ready yet")
       else
         do idir = 1, ndim
            call gradient(tmp(ixG^T),ixG^LL,ixM^LL,idir,grad(ixG^T, idir))
            ! Apply the correction B* = B - gradient(phi)
            ps(igrid)%prim(ixM^T, Bvec(idir)) = ps(igrid)%prim(ixM^T, Bvec(idir)) &
                                - grad(ixM^T, idir) / ps(igrid)%prim(ixM^T, psi_)**6 
         end do
       end if
    end do
  end subroutine grmhd_clean_divb_multigrid

end module mod_grmhd_phys_divb_mg
