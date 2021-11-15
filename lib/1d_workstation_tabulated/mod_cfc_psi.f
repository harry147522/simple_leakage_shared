module mod_cfc_psi
  use mod_cfc_parameters
  implicit none
  private

  ! public methods
  public :: cfc_solve_psi

  contains

  subroutine cfc_solve_psi(A2)
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics
!    use mod_grhd
!    use mod_grhd_ccsn
!    use mod_grhd_ccsn_leakage
    use mod_geometry
    use mod_eos, only: small_rho_thr

    double precision, intent(in)   :: A2(ixGlo1:ixGhi1,1:igridstail)
    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: src1(ixGlo1:ixGhi1)
    integer                        :: ixmin1,ixmax1

    ixmin1=ixMlo1-1;ixmax1=ixMhi1+1;

    mg%operator_type = mg_cfc_psi
    call mg_set_methods(mg)

    select case (coordinate)
    case (cartesian)
       do nb=1, 2*ndim
          select case (typeboundary(psi_, nb)) 
          case ('cont')
             ! outer boundary, use Robin B. C.
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case ('symm')
             ! inner boundary, use neumann
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             call mpistop&
("Error: The boundary conditions of psi/alp are not set correctly.")
          end select
       end do
    case (cylindrical)
       ! Neumann boundary condition for r=0
       mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
       mg%bc(1, mg_iphi)%bc_value = 0.0d0
       ! Robin boundary condition for r=rmax
       mg%bc(2, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(2, mg_iphi)%bc_value = 0.0d0
       
       
    case (spherical)
       ! Neumann boundary condition for r=0
       mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
       mg%bc(1, mg_iphi)%bc_value = 0.0d0
       ! Robin boundary condition for r=rmax
       mg%bc(2, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(2, mg_iphi)%bc_value = 0.0d0
       
       
    end select

    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! layer of ghost cells on grid leaves are not included here

       ! initial guess of (psi-1)
       mg%boxes(id)%cc(0:nc+1, mg_iphi) = (ps(igrid)%prim(ixmin1:ixmax1,&
           psi_) - 1.0d0)

       ! source term 1
       src1(ixMlo1:ixMhi1) = 2.0d0 * dpi * ( ps(igrid)%cons(ixMlo1:ixMhi1,&
           tau_) + ps(igrid)%cons(ixMlo1:ixMhi1, D_) ) 
       mg%boxes(id)%cc(1:nc, mg_itmp1) = src1(ixMlo1:ixMhi1)

       ! source term 2, A2
       mg%boxes(id)%cc(1:nc, mg_itmp2) = 1.25d-1 * A2(ixMlo1:ixMhi1, iigrid)

       ! No RHS 
       mg%boxes(id)%cc(1:nc, mg_irhs) = 0.0d0
    end do

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, .True., max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       !write(*,*) 'solving psi:', mg_it, res
       if (res <= cfc_tol(2)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving psi:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge psi."
          write(*,*) "! it=", it,"N = ", mg_it, " Res = ", res
          write(*,*)&
              "! Maybe the cfc_tol is too small or cfc_it_max is too small"
       end if
    end if

    ! copy the data from MG solver
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       ! Include one layer of ghost cells on grid leaves
       ! the extra layer is needed, as we need to workout d/dr (alp/psi**6) !
       ps(igrid)%prim(ixmin1:ixmax1, psi_) = mg%boxes(id)%cc(0:nc+1,&
           mg_iphi) + 1.0d0
    end do

  end subroutine cfc_solve_psi

end module mod_cfc_psi
