module mod_grrhd_phys_flux
  use mod_physics
  use mod_grrhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grrhd_phys_flux_init

contains
  !> Initialize the module
  subroutine grrhd_phys_flux_init()
    phys_get_cbounds         => grrhd_get_cbounds
    phys_get_flux            => grrhd_get_flux
  end subroutine grrhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine grrhd_get_cbounds(consL, consR, primL, primR, x, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, idim, cmax, cmin)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative left and right status
    double precision, intent(in)    :: consL(ixImin1:ixImax1, 1:ncons),&
        consR(ixImin1:ixImax1, 1:ncons)
    ! primitive left and right status
    double precision, intent(in)    :: primL(ixImin1:ixImax1, 1:nprim),&
        primR(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1)

    double precision, dimension(ixImin1:ixImax1) :: tmp1,tmp2
    double precision, dimension(ixImin1:ixImax1,1:2) :: lambdaL, lambdaR

    call grrhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        primL(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lambdaL(ixImin1:ixImax1,1:2))
    call grrhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        primR(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lambdaR(ixImin1:ixImax1,1:2))
  
    tmp1(ixOmin1:ixOmax1)=max(0.0d0, lambdaL(ixOmin1:ixOmax1,1),&
        lambdaR(ixOmin1:ixOmax1,1) )
    tmp2(ixOmin1:ixOmax1)=min(0.0d0, lambdaL(ixOmin1:ixOmax1,2),&
        lambdaR(ixOmin1:ixOmax1,2) )

    if(present(cmin)) then
      cmax(ixOmin1:ixOmax1) = tmp1(ixOmin1:ixOmax1)
      cmin(ixOmin1:ixOmax1) = tmp2(ixOmin1:ixOmax1)
    else
      cmax(ixOmin1:ixOmax1) = max(abs(tmp1(ixOmin1:ixOmax1)),&
          abs(tmp2(ixOmin1:ixOmax1)))
    end if

  end subroutine grrhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grrhd_get_flux(cons, prim, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixImin1:ixImax1, 1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1, 1:ncons)

    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision                :: v_hat(ixImin1:ixImax1) 
    double precision                :: lfac(ixImin1:ixImax1) 
    double precision                :: alppsi6(ixImin1:ixImax1) 
    integer                         :: idir, jdir, kdir
    double precision                :: Pij(ixImin1:ixImax1,1:3,1:3)
    double precision                :: P_ij(ixImin1:ixImax1,1:3,1:3)

    call grrhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3),  lfac=lfac(ixImin1:ixImax1),&
        Pij=P_ij(ixImin1:ixImax1,1:3,1:3)  )

    ! v_hat = v * alp - beta
    v_hat(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
       W_vel(idim)) / lfac(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
       alp_) - prim(ixOmin1:ixOmax1,beta(idim))

    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
       alp_) * prim(ixOmin1:ixOmax1,psi_)**6

    if ( evolve_hydro ) then 
       ! Density flux
       f(ixOmin1:ixOmax1, D_) = v_hat(ixOmin1:ixOmax1) * cons(ixOmin1:ixOmax1,&
           D_)
   
       ! Momentum flux
       do idir = 1, ndir
         f(ixOmin1:ixOmax1, mom(idir)) = v_hat(ixOmin1:ixOmax1) * &
            cons(ixOmin1:ixOmax1, mom(idir))
       end do
       f(ixOmin1:ixOmax1, mom(idim)) = f(ixOmin1:ixOmax1,&
           mom(idim)) + alppsi6(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
           press_) 
   
       ! Energy flux 
       f(ixOmin1:ixOmax1, tau_) = v_hat(ixOmin1:ixOmax1) * &
          cons(ixOmin1:ixOmax1, tau_) + alppsi6(ixOmin1:ixOmax1) * &
          prim(ixOmin1:ixOmax1, press_) * prim(ixOmin1:ixOmax1,&
          W_vel(idim)) / lfac(ixOmin1:ixOmax1)
    else
       f(ixOmin1:ixOmax1, D_) = 0.0d0
       f(ixOmin1:ixOmax1, mom(1:ndir)) = 0.0d0
       f(ixOmin1:ixOmax1, tau_) = 0.0d0
    end if

    if ( use_radiation ) then 
       ! first we get P^{ij}, here we denote it as P_ij
       !call grrhd_get_Pij(ixI^L, ixO^L, cons, prim, x, P_ij)
       ! here we need P^j_i, here we denote it as Pij
       Pij = 0.0d0
       do idir = 1, ndir; do jdir = 1, ndir
          do kdir = 1, ndir
             Pij(ixOmin1:ixOmax1,jdir,idir) = Pij(ixOmin1:ixOmax1,jdir,&
                idir) + gamma(ixOmin1:ixOmax1,kdir,&
                idir) * P_ij(ixOmin1:ixOmax1, jdir, kdir)
          end do
       end do; end do
   
       ! Momentum flux
       do idir = 1, ndir
         f(ixOmin1:ixOmax1, nu_Fcons(idir)) = prim(ixOmin1:ixOmax1,&
             alp_) * Pij(ixOmin1:ixOmax1, idim, idir) - prim(ixOmin1:ixOmax1,&
             beta(idim)) * cons(ixOmin1:ixOmax1, nu_Fcons(idir))
       end do
   
       ! Energy flux 
       f(ixOmin1:ixOmax1, nu_Econs) = prim(ixOmin1:ixOmax1,&
          alp_) * cons(ixOmin1:ixOmax1, nu_Fcons(idim)) / &
          gamma(ixOmin1:ixOmax1, idim, idim)  - prim(ixOmin1:ixOmax1,&
           beta(idim)) * cons(ixOmin1:ixOmax1,nu_Econs) 
    end if

  end subroutine grrhd_get_flux

  subroutine grrhd_modify_flux(ixImin1,ixImax1, ixOmin1,ixOmax1, consL, consR,&
      primL, primR, xi, s, idir, f)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(inout) :: consL(ixImin1:ixImax1,1:ncons),&
        consR(ixImin1:ixImax1,1:ncons)
    double precision, intent(inout) :: primL(ixImin1:ixImax1,1:nprim),&
        primR(ixImin1:ixImax1,1:nprim)
    double precision, intent(in)    :: xi(ixImin1:ixImax1, 1:1)
    type(state)                     :: s
    double precision, intent(inout) :: f(ixImin1:ixImax1, ncons)

    double precision                :: weight_factor(ixImin1:ixImax1, 1:2)
    double precision                :: opacities(ixImin1:ixImax1, 1:3) !they are eta, and two kappas
    double precision                :: kappa_total(ixImin1:ixImax1) !total opacity

    !call usr_get_resistivity(ixI^L, ixO^L, cons(ixI^S, 1:ncons), x(ixI^S, 1:ndim), eta(ixI^S))

  end subroutine grrhd_modify_flux

end module mod_grrhd_phys_flux
