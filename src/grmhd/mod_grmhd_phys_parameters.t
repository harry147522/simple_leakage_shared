module mod_grmhd_phys_parameters
  use mod_global_parameters, only: std_len
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: use_GR = .false.

  !-------------------------------------------------------------------!
  ! Parameters for con2prim
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tolerance = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_max = 50000
  !> maximum Lorentz factor
  double precision             :: lfac_max = 1.1d1
  double precision             :: v_max

  !-------------------------------------------------------------------!
  ! Parameters for divergence handle
  !-------------------------------------------------------------------!
  logical                      :: divb_4thorder = .False.
  integer                      :: type_divb = 0
  ! DivB cleaning methods
  integer, parameter           :: none = 0
  integer, parameter           :: divb_multigrid = 1
  integer, parameter           :: divb_GLM = 2
  integer, parameter           :: divb_ct = 3

  !-------------------------------------------------------------------!
  ! Parameters for divergence with Multigrid
  !-------------------------------------------------------------------!
  ! N smooths will be applied for 1: cycling up, 2: cycling down
  integer, dimension(1:2)                   :: divB_mg_n_cycle = (/5,5/)
  logical                                   :: divB_mg_redblack = .True.
  double precision                          :: divB_mg_tol = 1.0d-8
  integer                                   :: divB_mg_it_max = 50

  !-------------------------------------------------------------------!
  ! Parameters for divergence with GLM
  !-------------------------------------------------------------------!
  ! N smooths will be applied for 1: cycling up, 2: cycling down
  double precision                          :: divB_glm_kappa = 1.0d0

  !-------------------------------------------------------------------!
  ! Parameters for divergence with constrant transport
  !-------------------------------------------------------------------!
  character(len=std_len)                    :: type_ct = ""


  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grmhd_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_physics
    use mod_eos
    double precision, intent(inout)           :: prim(1:nprim)
    integer                                   :: ix^D

    ! rho and eps are given, update the rest of the primitive variables
    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_))
    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_))
    ! strictly require cs2 is physical
    if ( prim(cs2_) >= 1.0d0 ) then
       prim(cs2_) = 0.0d0
    else
       prim(cs2_) = max( prim(cs2_) , 0.0d0)
    end if
  end subroutine grmhd_update_eos_one_point

  !> get some useful variables from primitive
  subroutine grmhd_get_intermediate_variables(ixI^L, ixO^L, prim, x, &
             gamma, lfac2, lfac, v2, v_hat, b2, B_dot_v, Bvec2, bmu, Ptot, htot )
    use mod_global_parameters
    use mod_physics
    use mod_geometry
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: prim(ixI^S, 1:nprim)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(out), optional :: B_dot_v(ixI^S)   
    double precision, intent(out), optional :: v_hat(ixI^S,1:ndir)   
    double precision, intent(out), optional :: bmu(ixI^S,0:ndir)   ! projection of B^mu along fluid four-velocity u^nu
    double precision, intent(out), optional :: b2(ixI^S)        ! 
    double precision, intent(out), optional :: Bvec2(ixI^S)        ! 
    double precision, intent(out), optional :: gamma(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: v2(ixI^S)        ! 
    double precision, intent(out), optional :: lfac2(ixI^S) ! Lorentz factor square
    double precision, intent(out), optional :: lfac(ixI^S) ! Lorentz factor
    double precision, intent(out), optional :: htot(ixI^S) ! modified enthalpy:(h + b2/rho)
    double precision, intent(out), optional :: Ptot(ixI^S) ! total pressure

    integer                                 :: idir
    double precision                        :: B_dot_v_tmp(ixI^S)   
    double precision                        :: bmu_tmp(ixI^S,0:ndir)   ! projection of B^mu along fluid four-velocity u^nu
    double precision                        :: b2_tmp(ixI^S)        ! 
    double precision                        :: W2v2(ixI^S)        ! 
    double precision                        :: v_hat_tmp(ixI^S,1:ndir) 
    double precision                        :: lfac2_tmp(ixI^S) ! Lorentz factor square
    double precision                        :: lfac_tmp(ixI^S) ! Lorentz factor
    double precision                        :: gamma_tmp(ixI^S,1:3,1:3)

    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    if ( present(gamma) ) then
       gamma(ixO^S,1:3,1:3) = gamma_tmp(ixO^S,1:3,1:3) 
    end if

    W2v2 = 0.0d0
    ! calculate W^2 * v^2 first
    do idir = 1, ndir
       W2v2(ixO^S) = W2v2(ixO^S) + gamma_tmp(ixO^S,idir,idir)*prim(ixO^S, W_vel(idir))**2
    end do

    ! Calculate the Lorentz factor from velocity
    lfac2_tmp(ixO^S) = W2v2(ixO^S) + 1.0d0 
    lfac_tmp(ixO^S) = dsqrt( lfac2_tmp(ixO^S) )
    if ( present(lfac2) ) lfac2(ixO^S) = lfac2_tmp(ixO^S)
    if ( present(lfac) ) lfac(ixO^S) = lfac_tmp(ixO^S)
    if ( present(v2) ) v2(ixO^S) = W2v2(ixO^S) / lfac2_tmp(ixO^S)

    ! v_hat^i = v^i * alp - beta
    do idir = 1, ndir
       v_hat_tmp(ixO^S, idir) = prim(ixO^S, alp_) * prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) - prim(ixO^S, beta(idir))
    end do
    if ( present(v_hat) ) v_hat(ixO^S,1:ndir) = v_hat_tmp(ixO^S,1:ndir)

    B_dot_v_tmp(ixO^S) = 0.0d0
    do idir = 1, ndir
       B_dot_v_tmp(ixO^S) = B_dot_v_tmp(ixO^S) & 
              + prim(ixO^S, Bvec(idir)) * gamma_tmp(ixO^S,idir,idir) & 
               * prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) 
    end do
    if ( present(B_dot_v) ) then
       B_dot_v(ixO^S) = B_dot_v_tmp(ixO^S)
    end if

    b2_tmp(ixO^S) = 0.0d0
    do idir = 1, ndir
       b2_tmp(ixO^S) = b2_tmp(ixO^S) & 
              + prim(ixO^S, Bvec(idir))**2 * gamma_tmp(ixO^S,idir,idir)
    end do
    if ( present(Bvec2) ) then
       Bvec2(ixO^S) = b2_tmp(ixO^S)
    end if
    b2_tmp(ixO^S) = b2_tmp(ixO^S) / lfac2_tmp(ixO^S) + B_dot_v_tmp(ixO^S)**2
    if ( present(b2) ) then
       b2(ixO^S) = b2_tmp(ixO^S)
    end if

    if ( present(bmu) ) then
       ! Calculate the projection of B^mu along fluid four-velocity u^nu
       bmu(ixO^S,0) = B_dot_v_tmp(ixO^S) * lfac_tmp(ixO^S) / prim(ixO^S, alp_)
       do idir = 1, ndir
          bmu(ixO^S,idir) = prim(ixO^S, Bvec(idir)) / lfac_tmp(ixO^S) &
                     + bmu(ixO^S,0) * v_hat_tmp(ixO^S, idir)
       end do
    end if

    if ( present(Ptot) ) then
       ! Calculate the total pressure
       Ptot(ixO^S) = prim(ixO^S, press_) + 0.5d0 * b2_tmp(ixO^S)
    end if

    if ( present(htot) ) then
       ! Calculate the magnetically modified specific enthalpy
       htot(ixO^S) = 1.0d0 + prim(ixO^S, eps_) & 
                + ( prim(ixO^S, press_) + b2_tmp(ixO^S) ) / prim(ixO^S, rho_) 
    end if
  end subroutine grmhd_get_intermediate_variables

  subroutine grmhd_get_lambda(ixI^L, ixO^L, idim, prim, x, lambda)
     use mod_global_parameters
     use mod_physics
     integer, intent(in)             :: ixI^L, ixO^L, idim
     double precision, intent(in)    :: prim(ixI^S, 1:nprim)
     double precision, intent(in)    :: x(ixI^S, 1:ndim)
     double precision, intent(out)   :: lambda(ixI^S, 1:2)

     double precision                :: gamma(ixI^S,1:3,1:3)
     double precision                :: htot(ixI^S)    ! modified enthalpy:(h + b2/rho)
     double precision                :: cs2(ixI^S), ca2(ixI^S)
     double precision                :: v2(ixI^S)
     double precision                :: lfac(ixI^S)
     double precision                :: b2(ixI^S)
     double precision                :: a2(ixI^S) ! upper bound for the fast wave speed
     double precision                :: tmp1(ixI^S), tmp2(ixI^S)
     double precision                :: vel(ixI^S, 1:ndir)
     double precision                :: clight(ixI^S)
     integer                         :: idir

     ! sound speed
     cs2(ixO^S) = prim(ixO^S, cs2_)

    call grmhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                v2=v2(ixI^S), lfac=lfac(ixI^S), &
                b2=b2(ixI^S), &
                htot=htot(ixI^S) )

     do idir = 1, ndir
        vel(ixO^S,idir) = prim(ixO^S, W_vel(idir)) / lfac(ixO^S)
     end do

     ! Calculate Alfven speed
     ca2(ixO^S) =  b2(ixO^S) / ( prim(ixO^S, rho_) * htot(ixO^S) )
 
     ! upper bound for the fast wave speed
     a2(ixO^S) = cs2(ixO^S) + ca2(ixO^S) - ca2(ixO^S) * cs2(ixO^S)
 
     tmp1(ixO^S) = ( 1.0d0 - a2(ixO^S) ) * vel(ixO^S, idim) 
     tmp2(ixO^S) = dsqrt( a2(ixO^S) * ( 1.0d0 - v2(ixO^S) ) * &
                         ( ( 1.0d0 - v2(ixO^S) * a2(ixO^S) ) / gamma(ixO^S,idim,idim) &
                            - ( 1.0d0 - a2(ixO^S) ) * vel(ixO^S, idim)**2 ) )

     lambda(ixO^S,1) = ( tmp1(ixO^S) + tmp2(ixO^S) ) / ( 1.0d0 - v2(ixO^S) * a2(ixO^S) ) 
     lambda(ixO^S,2) = ( tmp1(ixO^S) - tmp2(ixO^S) ) / ( 1.0d0 - v2(ixO^S) * a2(ixO^S) ) 

     ! limit with speed of light
     clight(ixO^S) = dsqrt( 1.0d0 / gamma(ixO^S, idim, idim) )
     lambda(ixO^S,1) = max( min( lambda(ixO^S,1), clight(ixO^S) ), -clight(ixO^S) )
     lambda(ixO^S,2) = max( min( lambda(ixO^S,2), clight(ixO^S) ), -clight(ixO^S) )

     lambda(ixO^S,1) = prim(ixO^S, alp_) * lambda(ixO^S,1) - prim(ixO^S, beta(idim))
     lambda(ixO^S,2) = prim(ixO^S, alp_) * lambda(ixO^S,2) - prim(ixO^S, beta(idim))

     ! when using GLM, we need to include two additional modes with speed of
     ! light, to without violating causality
     if (type_divb == divb_glm) then
        ! calculate sqrt(g^ii)
        tmp1(ixO^S) = dsqrt( 1.0d0 / gamma(ixO^S, idim, idim) &
                     - prim(ixO^S,beta(idim))**2 / prim(ixO^S, alp_)  )
        lambda(ixO^S,1) = max( lambda(ixO^S,1), &
                           tmp1(ixO^S) + prim(ixO^S,beta(idim)) / prim(ixO^S, alp_) )
        lambda(ixO^S,2) = min( lambda(ixO^S,2), &
                          -tmp1(ixO^S) + prim(ixO^S,beta(idim)) / prim(ixO^S, alp_) )
     end if
  end subroutine grmhd_get_lambda

  !> Calculate div B within ixO
  subroutine grmhd_get_divb(ixI^L,ixO^L, Bvector, divb, fourthorder)
    use mod_global_parameters
    use mod_physics
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: bvector(ixI^S,1:ndir)
    double precision, intent(inout) :: divb(ixI^S)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: divb_corner(ixI^S), sign
    double precision                   :: aux_vol(ixI^S)
    integer                            :: ixC^L, idir, ic^D, ix^L

    if ( stagger_grid ) then
      stop "fixme"
      divb=0.d0
      do idir=1,ndim
        ixC^L=ixO^L-kr(idir,^D);
        divb(ixO^S)=divb(ixO^S)+block%prims(ixO^S,idir)*block%surfaceC(ixO^S,idir)-&
                                block%prims(ixC^S,idir)*block%surfaceC(ixC^S,idir)
      end do
      divb(ixO^S)=divb(ixO^S)/block%dvolume(ixO^S)
    else
      !bvector(ixI^S,:)=cons(ixI^S,Bcons(:))
      !select case(typediv)
      !case("central")
        call divvector(bvector(ixI^S,1:ndir),ixI^L,ixO^L,divb(ixI^S),fourthorder)
      !case("limited")
        !call divvectorS(bvec,ixI^L,ixO^L,divb)
      !end select
    end if
  end subroutine grmhd_get_divb

end module mod_grmhd_phys_parameters
