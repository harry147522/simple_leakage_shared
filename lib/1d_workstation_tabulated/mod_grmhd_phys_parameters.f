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
    integer                                   :: ix1

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
  subroutine grmhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
      prim, x, gamma, lfac2, lfac, v2, v_hat, b2, B_dot_v, Bvec2, bmu, Ptot,&
      htot )
    use mod_global_parameters
    use mod_physics
    use mod_geometry
    integer, intent(in)                     :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(in)            :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)            :: x(ixImin1:ixImax1, 1:ndim)

    double precision, intent(out), optional :: B_dot_v(ixImin1:ixImax1)   
    double precision, intent(out), optional :: v_hat(ixImin1:ixImax1,&
       1:ndir)   
    double precision, intent(out), optional :: bmu(ixImin1:ixImax1,0:ndir) !projection of Bmu along fluid four-velocity unu
    double precision, intent(out), optional :: b2(ixImin1:ixImax1)        ! 
    double precision, intent(out), optional :: Bvec2(ixImin1:ixImax1) !
    double precision, intent(out), optional :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision, intent(out), optional :: v2(ixImin1:ixImax1)        ! 
    double precision, intent(out), optional :: lfac2(ixImin1:ixImax1) !Lorentz factor square
    double precision, intent(out), optional :: lfac(ixImin1:ixImax1) !Lorentz factor
    double precision, intent(out), optional :: htot(ixImin1:ixImax1) !modified enthalpy:(h + b2/rho)
    double precision, intent(out), optional :: Ptot(ixImin1:ixImax1) !total pressure

    integer                                 :: idir
    double precision                        :: B_dot_v_tmp(ixImin1:ixImax1)   
    double precision                        :: bmu_tmp(ixImin1:ixImax1,0:ndir) !projection of Bmu along fluid four-velocity unu
    double precision                        :: b2_tmp(ixImin1:ixImax1) !
    double precision                        :: W2v2(ixImin1:ixImax1)        ! 
    double precision                        :: v_hat_tmp(ixImin1:ixImax1,&
       1:ndir) 
    double precision                        :: lfac2_tmp(ixImin1:ixImax1) !Lorentz factor square
    double precision                        :: lfac_tmp(ixImin1:ixImax1) !Lorentz factor
    double precision                        :: gamma_tmp(ixImin1:ixImax1,1:3,&
       1:3)

    ! get the metric
    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma_tmp(ixImin1:ixImax1,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixOmin1:ixOmax1,idir,idir) = gamma_tmp(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, psi_)**4 
    end do
    if ( present(gamma) ) then
       gamma(ixOmin1:ixOmax1,1:3,1:3) = gamma_tmp(ixOmin1:ixOmax1,1:3,1:3) 
    end if

    W2v2 = 0.0d0
    ! calculate W^2 * v^2 first
    do idir = 1, ndir
       W2v2(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) + &
          gamma_tmp(ixOmin1:ixOmax1,idir,idir)*prim(ixOmin1:ixOmax1,&
           W_vel(idir))**2
    end do

    ! Calculate the Lorentz factor from velocity
    lfac2_tmp(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) + 1.0d0 
    lfac_tmp(ixOmin1:ixOmax1) = dsqrt( lfac2_tmp(ixOmin1:ixOmax1) )
    if ( present(lfac2) ) lfac2(ixOmin1:ixOmax1) = lfac2_tmp(ixOmin1:ixOmax1)
    if ( present(lfac) ) lfac(ixOmin1:ixOmax1) = lfac_tmp(ixOmin1:ixOmax1)
    if ( present(v2) ) v2(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) / &
       lfac2_tmp(ixOmin1:ixOmax1)

    ! v_hat^i = v^i * alp - beta
    do idir = 1, ndir
       v_hat_tmp(ixOmin1:ixOmax1, idir) = prim(ixOmin1:ixOmax1,&
           alp_) * prim(ixOmin1:ixOmax1, W_vel(idir)) / &
          lfac_tmp(ixOmin1:ixOmax1) - prim(ixOmin1:ixOmax1, beta(idir))
    end do
    if ( present(v_hat) ) v_hat(ixOmin1:ixOmax1,&
       1:ndir) = v_hat_tmp(ixOmin1:ixOmax1,1:ndir)

    B_dot_v_tmp(ixOmin1:ixOmax1) = 0.0d0
    do idir = 1, ndir
       B_dot_v_tmp(ixOmin1:ixOmax1) = B_dot_v_tmp(ixOmin1:ixOmax1) + &
          prim(ixOmin1:ixOmax1, Bvec(idir)) * gamma_tmp(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, W_vel(idir)) / &
          lfac_tmp(ixOmin1:ixOmax1) 
    end do
    if ( present(B_dot_v) ) then
       B_dot_v(ixOmin1:ixOmax1) = B_dot_v_tmp(ixOmin1:ixOmax1)
    end if

    b2_tmp(ixOmin1:ixOmax1) = 0.0d0
    do idir = 1, ndir
       b2_tmp(ixOmin1:ixOmax1) = b2_tmp(ixOmin1:ixOmax1) + &
          prim(ixOmin1:ixOmax1, Bvec(idir))**2 * gamma_tmp(ixOmin1:ixOmax1,&
          idir,idir)
    end do
    if ( present(Bvec2) ) then
       Bvec2(ixOmin1:ixOmax1) = b2_tmp(ixOmin1:ixOmax1)
    end if
    b2_tmp(ixOmin1:ixOmax1) = b2_tmp(ixOmin1:ixOmax1) / &
       lfac2_tmp(ixOmin1:ixOmax1) + B_dot_v_tmp(ixOmin1:ixOmax1)**2
    if ( present(b2) ) then
       b2(ixOmin1:ixOmax1) = b2_tmp(ixOmin1:ixOmax1)
    end if

    if ( present(bmu) ) then
       ! Calculate the projection of B^mu along fluid four-velocity u^nu
       bmu(ixOmin1:ixOmax1,0) = B_dot_v_tmp(ixOmin1:ixOmax1) * &
          lfac_tmp(ixOmin1:ixOmax1) / prim(ixOmin1:ixOmax1, alp_)
       do idir = 1, ndir
          bmu(ixOmin1:ixOmax1,idir) = prim(ixOmin1:ixOmax1,&
              Bvec(idir)) / lfac_tmp(ixOmin1:ixOmax1) + bmu(ixOmin1:ixOmax1,&
             0) * v_hat_tmp(ixOmin1:ixOmax1, idir)
       end do
    end if

    if ( present(Ptot) ) then
       ! Calculate the total pressure
       Ptot(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
           press_) + 0.5d0 * b2_tmp(ixOmin1:ixOmax1)
    end if

    if ( present(htot) ) then
       ! Calculate the magnetically modified specific enthalpy
       htot(ixOmin1:ixOmax1) = 1.0d0 + prim(ixOmin1:ixOmax1,&
           eps_) + ( prim(ixOmin1:ixOmax1,&
           press_) + b2_tmp(ixOmin1:ixOmax1) ) / prim(ixOmin1:ixOmax1, rho_) 
    end if
  end subroutine grmhd_get_intermediate_variables

  subroutine grmhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim, prim, x,&
      lambda)
     use mod_global_parameters
     use mod_physics
     integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
     double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
     double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
     double precision, intent(out)   :: lambda(ixImin1:ixImax1, 1:2)

     double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
     double precision                :: htot(ixImin1:ixImax1) !modified enthalpy:(h + b2/rho)
     double precision                :: cs2(ixImin1:ixImax1),&
         ca2(ixImin1:ixImax1)
     double precision                :: v2(ixImin1:ixImax1)
     double precision                :: lfac(ixImin1:ixImax1)
     double precision                :: b2(ixImin1:ixImax1)
     double precision                :: a2(ixImin1:ixImax1) !upper bound for the fast wave speed
     double precision                :: tmp1(ixImin1:ixImax1),&
         tmp2(ixImin1:ixImax1)
     double precision                :: vel(ixImin1:ixImax1, 1:ndir)
     double precision                :: clight(ixImin1:ixImax1)
     integer                         :: idir

     ! sound speed
     cs2(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1, cs2_)

    call grmhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3), v2=v2(ixImin1:ixImax1),&
        lfac=lfac(ixImin1:ixImax1), b2=b2(ixImin1:ixImax1),&
        htot=htot(ixImin1:ixImax1) )

     do idir = 1, ndir
        vel(ixOmin1:ixOmax1,idir) = prim(ixOmin1:ixOmax1,&
            W_vel(idir)) / lfac(ixOmin1:ixOmax1)
     end do

     ! Calculate Alfven speed
     ca2(ixOmin1:ixOmax1) =  b2(ixOmin1:ixOmax1) / ( prim(ixOmin1:ixOmax1,&
         rho_) * htot(ixOmin1:ixOmax1) )
 
     ! upper bound for the fast wave speed
     a2(ixOmin1:ixOmax1) = cs2(ixOmin1:ixOmax1) + ca2(ixOmin1:ixOmax1) - &
        ca2(ixOmin1:ixOmax1) * cs2(ixOmin1:ixOmax1)
 
     tmp1(ixOmin1:ixOmax1) = ( 1.0d0 - a2(ixOmin1:ixOmax1) ) * &
        vel(ixOmin1:ixOmax1, idim) 
     tmp2(ixOmin1:ixOmax1) = dsqrt( a2(ixOmin1:ixOmax1) * ( 1.0d0 - &
        v2(ixOmin1:ixOmax1) ) * ( ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
        a2(ixOmin1:ixOmax1) ) / gamma(ixOmin1:ixOmax1,idim,&
        idim) - ( 1.0d0 - a2(ixOmin1:ixOmax1) ) * vel(ixOmin1:ixOmax1,&
         idim)**2 ) )

     lambda(ixOmin1:ixOmax1,1) = ( tmp1(ixOmin1:ixOmax1) + &
        tmp2(ixOmin1:ixOmax1) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
        a2(ixOmin1:ixOmax1) ) 
     lambda(ixOmin1:ixOmax1,2) = ( tmp1(ixOmin1:ixOmax1) - &
        tmp2(ixOmin1:ixOmax1) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
        a2(ixOmin1:ixOmax1) ) 

     ! limit with speed of light
     clight(ixOmin1:ixOmax1) = dsqrt( 1.0d0 / gamma(ixOmin1:ixOmax1, idim,&
         idim) )
     lambda(ixOmin1:ixOmax1,1) = max( min( lambda(ixOmin1:ixOmax1,1),&
         clight(ixOmin1:ixOmax1) ), -clight(ixOmin1:ixOmax1) )
     lambda(ixOmin1:ixOmax1,2) = max( min( lambda(ixOmin1:ixOmax1,2),&
         clight(ixOmin1:ixOmax1) ), -clight(ixOmin1:ixOmax1) )

     lambda(ixOmin1:ixOmax1,1) = prim(ixOmin1:ixOmax1,&
         alp_) * lambda(ixOmin1:ixOmax1,1) - prim(ixOmin1:ixOmax1, beta(idim))
     lambda(ixOmin1:ixOmax1,2) = prim(ixOmin1:ixOmax1,&
         alp_) * lambda(ixOmin1:ixOmax1,2) - prim(ixOmin1:ixOmax1, beta(idim))

     ! when using GLM, we need to include two additional modes with speed of
     ! light, to without violating causality
     if (type_divb == divb_glm) then
        ! calculate sqrt(g^ii)
        tmp1(ixOmin1:ixOmax1) = dsqrt( 1.0d0 / gamma(ixOmin1:ixOmax1, idim,&
            idim) - prim(ixOmin1:ixOmax1,beta(idim))**2 / prim(ixOmin1:ixOmax1,&
            alp_)  )
        lambda(ixOmin1:ixOmax1,1) = max( lambda(ixOmin1:ixOmax1,1),&
            tmp1(ixOmin1:ixOmax1) + prim(ixOmin1:ixOmax1,&
           beta(idim)) / prim(ixOmin1:ixOmax1, alp_) )
        lambda(ixOmin1:ixOmax1,2) = min( lambda(ixOmin1:ixOmax1,2),&
            -tmp1(ixOmin1:ixOmax1) + prim(ixOmin1:ixOmax1,&
           beta(idim)) / prim(ixOmin1:ixOmax1, alp_) )
     end if
  end subroutine grmhd_get_lambda

  !> Calculate div B within ixO
  subroutine grmhd_get_divb(ixImin1,ixImax1,ixOmin1,ixOmax1, Bvector, divb,&
      fourthorder)
    use mod_global_parameters
    use mod_physics
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: bvector(ixImin1:ixImax1,1:ndir)
    double precision, intent(inout) :: divb(ixImin1:ixImax1)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: divb_corner(ixImin1:ixImax1), sign
    double precision                   :: aux_vol(ixImin1:ixImax1)
    integer                            :: ixCmin1,ixCmax1, idir, ic1, ixmin1,&
       ixmax1

    if ( stagger_grid ) then
      stop "fixme"
      divb=0.d0
      do idir=1,ndim
        ixCmin1=ixOmin1-kr(idir,1);ixCmax1=ixOmax1-kr(idir,1);
        divb(ixOmin1:ixOmax1)=divb(ixOmin1:ixOmax1)+&
           block%prims(ixOmin1:ixOmax1,idir)*block%surfaceC(ixOmin1:ixOmax1,&
           idir)-block%prims(ixCmin1:ixCmax1,&
           idir)*block%surfaceC(ixCmin1:ixCmax1,idir)
      end do
      divb(ixOmin1:ixOmax1)=divb(ixOmin1:ixOmax1)/block%dvolume(&
         ixOmin1:ixOmax1)
    else
      !bvector(ixI^S,:)=cons(ixI^S,Bcons(:))
      !select case(typediv)
      !case("central")
        call divvector(bvector(ixImin1:ixImax1,1:ndir),ixImin1,ixImax1,ixOmin1,&
           ixOmax1,divb(ixImin1:ixImax1),fourthorder)
      !case("limited")
        !call divvectorS(bvec,ixI^L,ixO^L,divb)
      !end select
    end if
  end subroutine grmhd_get_divb

end module mod_grmhd_phys_parameters
