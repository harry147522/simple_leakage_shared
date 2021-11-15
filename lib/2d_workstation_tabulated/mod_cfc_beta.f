module mod_cfc_beta
  use mod_cfc_parameters
  implicit none
  private

  ! public methods
  public :: cfc_solve_beta

  contains

  subroutine cfc_solve_beta(solve_beta, Aij)
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics
    use mod_geometry

    logical, intent(in)            :: solve_beta
    double precision, intent(in)   :: Aij(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3,1:3,&
       1:igridstail)

    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0), res_old = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, idir, nb
    type(tree_node), pointer       :: pnode
    double precision               :: rhs(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3)
    double precision               :: n2r(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3)

    character(len=128)             :: phi_name = ""
    integer                        :: nvec(1:ndir)
    integer                        :: ixmin1,ixmin2,ixmax1,ixmax2

    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;

    if (solve_beta) then 
       phi_name = 'beta'
       nvec(1:ndir) = beta(1:ndir)
    else
       phi_name = 'X'
       nvec(1:ndir) = vecX(1:ndir)
    end if

    mg%operator_type = mg_cfc_beta
    call mg_set_methods(mg)

    select case (coordinate)
    case (cartesian)
       do nb = 1, 2*ndim
          do idir = 1, ndir
             select case (typeboundary(beta(idir), nb))
             case ('symm')
                ! inner boundary, but not reflecting:
                mg%bc(nb, mg_vec_iphi(idir))%bc_type = mg_bc_neumann
                mg%bc(nb, mg_vec_iphi(idir))%bc_value = 0.0d0
             case default
                ! outer boundary, or reflecting inner boundary:
                mg%bc(nb, mg_vec_iphi(idir))%bc_type = mg_bc_dirichlet
                mg%bc(nb, mg_vec_iphi(idir))%bc_value = 0.0d0
             end select
             !write(*,*) nb, idir, typeboundary(beta(idir), nb)
          end do
       end do
    case (cylindrical)
       ! B.C. for r=0
       mg%bc(1, mg_vec_iphi(1))%bc_type = mg_bc_dirichlet
       mg%bc(1, mg_vec_iphi(1))%bc_value = 0.0d0
       mg%bc(1, mg_vec_iphi(2))%bc_type = mg_bc_neumann
       mg%bc(1, mg_vec_iphi(2))%bc_value = 0.0d0
       mg%bc(1, mg_vec_iphi(3))%bc_type = mg_bc_dirichlet
       mg%bc(1, mg_vec_iphi(3))%bc_value = 0.0d0
   
       ! B.C. for r=rmax
       mg%bc(2, mg_vec_iphi(1:3))%bc_type = mg_bc_dirichlet
       mg%bc(2, mg_vec_iphi(1:3))%bc_value = 0.0d0
   
       
       ! B.C. for z=zmin
       if( abs(xprobmin2) < smalldouble ) then
          mg%bc(3, mg_vec_iphi(1))%bc_type = mg_bc_neumann
          mg%bc(3, mg_vec_iphi(2))%bc_type = mg_bc_dirichlet
          mg%bc(3, mg_vec_iphi(3))%bc_type = mg_bc_neumann
          mg%bc(3, mg_vec_iphi(1:3))%bc_value = 0.0d0
       else
          mg%bc(3, mg_vec_iphi(1:3))%bc_type = mg_bc_dirichlet
          mg%bc(3, mg_vec_iphi(1:3))%bc_value = 0.0d0
       end if
       ! B.C. for z = zmax
       mg%bc(4, mg_vec_iphi(1:3))%bc_type = mg_bc_dirichlet
       mg%bc(4, mg_vec_iphi(1:3))%bc_value = 0.0d0
      
       
    case (spherical)
       ! B.C. for r=0
       mg%bc(1, mg_vec_iphi(1:3))%bc_type = mg_bc_dirichlet
       mg%bc(1, mg_vec_iphi(1:3))%bc_value = 0.0d0
   
       ! B.C. for r=rmax
       mg%bc(2, mg_vec_iphi(1:3))%bc_type = mg_bc_dirichlet
       mg%bc(2, mg_vec_iphi(1:3))%bc_value = 0.0d0
   
       
       ! B.C. for theta=0
       mg%bc(3, mg_vec_iphi(1))%bc_type = mg_bc_neumann
       mg%bc(3, mg_vec_iphi(1))%bc_value = 0.0d0
       mg%bc(3, mg_vec_iphi(2))%bc_type = mg_bc_dirichlet
       mg%bc(3, mg_vec_iphi(2))%bc_value = 0.0d0
       mg%bc(3, mg_vec_iphi(3))%bc_type = mg_bc_dirichlet
       mg%bc(3, mg_vec_iphi(3))%bc_value = 0.0d0
   
       ! B.C. for theta=theta_max
       mg%bc(4, mg_vec_iphi(1))%bc_type = mg_bc_neumann
       mg%bc(4, mg_vec_iphi(1))%bc_value = 0.0d0
       mg%bc(4, mg_vec_iphi(2))%bc_type = mg_bc_dirichlet
       mg%bc(4, mg_vec_iphi(2))%bc_value = 0.0d0
       if(abs(xprobmax2-dpi)<=smalldouble) then
          ! if theta_max = pi
          mg%bc(4, mg_vec_iphi(3))%bc_type = mg_bc_dirichlet 
       else
          mg%bc(4, mg_vec_iphi(3))%bc_type = mg_bc_neumann 
       end if
       mg%bc(4, mg_vec_iphi(3))%bc_value = 0.0d0
      
       
    end select

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);

       ! make sure all the data is clean
       

       call cfc_beta_get_rhs(ps(igrid), ixGlo1,ixGlo2,ixGhi1,ixGhi2, ixMlo1,&
          ixMlo2,ixMhi1,ixMhi2, Aij(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3,1:3,&
          iigrid), rhs(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3), solve_beta)
       call get_natural2orthonormal(ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
          1:ndim),ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
          n2r(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3))
       ! store RHS, in orthonormal form
       do idir = 1, ndir
          mg%boxes(id)%cc(1:nc,1:nc, mg_vec_irhs(idir)) = rhs(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2, idir) * n2r(ixMlo1:ixMhi1,ixMlo2:ixMhi2, idir) 
       end do

       ! initial guess, in orthonormal form 
       do idir = 1, ndir
          mg%boxes(id)%cc(0:nc+1,0:nc+1, mg_vec_iphi(idir)) = &
             ps(igrid)%prim(ixmin1:ixmax1,ixmin2:ixmax2,&
              nvec(idir)) * n2r(ixmin1:ixmax1,ixmin2:ixmax2, idir)
       end do

       if ( coordinate /= cartesian ) then
          ! rescale the rhs: rhs = rhs * r**2
          do idir = 1,ndir
             mg%boxes(id)%cc(1:nc,1:nc, mg_vec_irhs(idir)) = &
                mg%boxes(id)%cc(1:nc,1:nc,&
                 mg_vec_irhs(idir)) * ps(igrid)%x(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 1)**2 
          end do
          
          if ( coordinate == spherical ) then
             ! rescale the rhs: rhs_final = rhs * r**2 sin_theta**2
             do idir = 1,ndir
                mg%boxes(id)%cc(1:nc,1:nc,&
                    mg_vec_irhs(idir)) = mg%boxes(id)%cc(1:nc,1:nc,&
                    mg_vec_irhs(idir)) * dsin(ps(igrid)%x(ixMlo1:ixMhi1,&
                   ixMlo2:ixMhi2, 2))**2 
             end do
          end if
         
       end if

    end do

    do mg_it = 1, cfc_it_max
       res_old = res
       call mg_fas_fmg(mg, .true., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       if ((res <= cfc_tol(3))) exit
       !if ((res <= cfc_tol(3)).or.(res == res_old)) exit
       if (mod(mg_it,cfc_print)==0) then
               if (mype==0) write(*,*) 'solving '//trim(adjustl(&
                  phi_name))//':',it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge "//trim(adjustl(&
             phi_name))//"."
          write(*,*) "! it=", it,"N = ", mg_it, " Res = ", res
          write(*,*)&
              "! Maybe the cfc_tol is too small or cfc_it_max is too small"
       end if
    end if

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
       call get_natural2orthonormal(ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
          1:ndim),ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
          n2r(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:3))
       ! copy the solution nvec from MG solver
       do idir = 1, ndir
          ps(igrid)%prim(ixmin1:ixmax1,ixmin2:ixmax2,&
              nvec(idir)) = mg%boxes(id)%cc(0:nc+1,0:nc+1,&
              mg_vec_iphi(idir)) / (n2r(ixmin1:ixmax1,ixmin2:ixmax2, idir))
       end do
    end do

  end subroutine cfc_solve_beta

  subroutine cfc_beta_get_rhs(ps_in, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, Aij, rhs, solve_beta)
    use mod_global_parameters
    use mod_physics, only: mom, alp_, beta, psi_
    use mod_geometry, only: partial_d, get_gamma_ij_hat
    logical, intent(in)             :: solve_beta
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state), intent(in)         :: ps_in
    double precision, intent(in)    :: Aij(ixImin1:ixImax1,ixImin2:ixImax2,1:3,&
       1:3)
    double precision, intent(out)   :: rhs(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:3)

    integer                         :: ixmin1,ixmin2,ixmax1,ixmax2
    integer                         :: idir, jdir
    double precision                :: fij(ixImin1:ixImax1,ixImin2:ixImax2,1:3,&
       1:3)
    double precision                :: alp_over_psi6(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                :: d_alp_over_psi6(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3)

    associate(prim=>ps_in%prim, cons=>ps_in%cons, x=>ps_in%x)

    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;

    rhs(ixImin1:ixImax1,ixImin2:ixImax2, 1:3) = 0.0d0
    call get_gamma_ij_hat(x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       fij(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3))

    ! this is 8 pi f^{ij}S_j
    do idir = 1, ndir
       rhs(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idir) = 8.0d0 * dpi * cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) / fij(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir, idir)
    end do

    if (.not.solve_beta) return

    ! source terms for beta equation
    alp_over_psi6(ixmin1:ixmax1,ixmin2:ixmax2) = prim(ixmin1:ixmax1,&
       ixmin2:ixmax2, alp_) / prim(ixmin1:ixmax1,ixmin2:ixmax2, psi_)**6

    do jdir = 1, ndir
       rhs(ixOmin1:ixOmax1,ixOmin2:ixOmax2,jdir) = rhs(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,jdir) * 2.0d0 * alp_over_psi6(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
    end do

    d_alp_over_psi6(ixImin1:ixImax1,ixImin2:ixImax2,1:3) = 0.0d0
    do idir = 1, ndim
       call partial_d(alp_over_psi6(ixImin1:ixImax1,ixImin2:ixImax2), ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir,&
           d_alp_over_psi6(ixImin1:ixImax1,ixImin2:ixImax2,idir))
    end do

    do idir = 1, ndir
       do jdir = 1, ndim
          rhs(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = rhs(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,idir) + 2.0d0 * d_alp_over_psi6(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, jdir) * Aij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             jdir)
       end do
    end do

    end associate
  end subroutine cfc_beta_get_rhs

end module mod_cfc_beta
