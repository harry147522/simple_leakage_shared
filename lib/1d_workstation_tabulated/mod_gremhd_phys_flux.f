module mod_gremhd_phys_flux
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: gremhd_phys_flux_init

contains
  !> Initialize the module
  subroutine gremhd_phys_flux_init()
    phys_get_cbounds         => gremhd_get_cbounds
    phys_get_flux            => gremhd_get_flux
  end subroutine gremhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine gremhd_get_cbounds(consL, consR, primL, primR, x, ixImin1,ixImax1,&
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

    double precision, dimension(ixImin1:ixImax1,1:3,1:3) :: gamma_hat
    double precision, dimension(ixImin1:ixImax1,1:2) :: lambdaL, lambdaR
    double precision, dimension(ixImin1:ixImax1,1:2) :: tmp_c

    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma_hat(ixImin1:ixImax1,1:3,1:3))
    lambdaL(ixOmin1:ixOmax1,1) = dsqrt( 1.0d0 / gamma_hat(ixOmin1:ixOmax1,&
        idim, idim)) / primL(ixOmin1:ixOmax1, psi_)**2
    lambdaL(ixOmin1:ixOmax1,2) = - lambdaL(ixOmin1:ixOmax1,1)
    lambdaL(ixOmin1:ixOmax1,1) = primL(ixOmin1:ixOmax1,&
        alp_) * lambdaL(ixOmin1:ixOmax1,1) - primL(ixOmin1:ixOmax1,&
        beta(idim))
    lambdaL(ixOmin1:ixOmax1,2) = primL(ixOmin1:ixOmax1,&
        alp_) * lambdaL(ixOmin1:ixOmax1,2) - primL(ixOmin1:ixOmax1,&
        beta(idim))

    lambdaR(ixOmin1:ixOmax1,1) = dsqrt( 1.0d0 / gamma_hat(ixOmin1:ixOmax1,&
        idim, idim)) / primR(ixOmin1:ixOmax1, psi_)**2 
    lambdaR(ixOmin1:ixOmax1,2) = - lambdaR(ixOmin1:ixOmax1,1)
    lambdaR(ixOmin1:ixOmax1,1) = primR(ixOmin1:ixOmax1,&
        alp_) * lambdaR(ixOmin1:ixOmax1,1) - primR(ixOmin1:ixOmax1,&
        beta(idim))
    lambdaR(ixOmin1:ixOmax1,2) = primR(ixOmin1:ixOmax1,&
        alp_) * lambdaR(ixOmin1:ixOmax1,2) - primR(ixOmin1:ixOmax1,&
        beta(idim))
    
    tmp_c(ixOmin1:ixOmax1,1)=max(0.0d0, lambdaL(ixOmin1:ixOmax1,1),&
        lambdaR(ixOmin1:ixOmax1,1) ) 
    tmp_c(ixOmin1:ixOmax1,2)=min(0.0d0, lambdaL(ixOmin1:ixOmax1,2),&
        lambdaR(ixOmin1:ixOmax1,2) )

    if(present(cmin)) then
      cmax(ixOmin1:ixOmax1) = tmp_c(ixOmin1:ixOmax1,1)
      cmin(ixOmin1:ixOmax1) = tmp_c(ixOmin1:ixOmax1,2)
    else
      cmax(ixOmin1:ixOmax1) = max(abs(tmp_c(ixOmin1:ixOmax1,1)),&
          abs(tmp_c(ixOmin1:ixOmax1,2)))
    end if
  end subroutine gremhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine gremhd_get_flux(cons, prim, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      idim, f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixImin1:ixImax1, 1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1, ncons)

    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision                :: v_hat(ixImin1:ixImax1, 1:ndir) 
    double precision                :: alppsi6(ixImin1:ixImax1) 
    double precision                :: h(ixImin1:ixImax1) ! 
    double precision                :: lfac(ixImin1:ixImax1) ! Lorentz factor
    double precision                :: sqrt_gamma(ixImin1:ixImax1)
    double precision                :: E2(ixImin1:ixImax1),&
        B2(ixImin1:ixImax1) 
    double precision                :: inv_sqrt_gamma(ixImin1:ixImax1),&
        inv_gamma_ii(ixImin1:ixImax1)
    integer                         :: idir, jdir, kdir

    call gremhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3),&
        sqrt_gamma=sqrt_gamma(ixImin1:ixImax1), lfac=lfac(ixImin1:ixImax1),&
        v_hat=v_hat(ixImin1:ixImax1,1:ndir), h=h(ixImin1:ixImax1),&
        E2=E2(ixImin1:ixImax1), B2=B2(ixImin1:ixImax1))

    where ( sqrt_gamma(ixOmin1:ixOmax1) > smalldouble )
       inv_sqrt_gamma(ixOmin1:ixOmax1) = 1.0d0 / sqrt_gamma(ixOmin1:ixOmax1)
       inv_gamma_ii(ixOmin1:ixOmax1) = 1.0d0 / gamma(ixOmin1:ixOmax1, idim,&
           idim)
    else where
       inv_sqrt_gamma(ixOmin1:ixOmax1) = 0.0d0
       inv_gamma_ii(ixOmin1:ixOmax1) = 0.0d0
    end where

    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
       alp_) * prim(ixOmin1:ixOmax1,psi_)**6

    if ( evolve_hydro ) then
       ! Density flux
       f(ixOmin1:ixOmax1, D_) = v_hat(ixOmin1:ixOmax1,&
           idim) * cons(ixOmin1:ixOmax1, D_)
   
       ! Momentum flux f^idim_idir
       do idir = 1, ndir
         f(ixOmin1:ixOmax1, mom(idir)) =  - prim(ixOmin1:ixOmax1,&
             beta(idim)) * cons(ixOmin1:ixOmax1,&
             mom(idir)) + alppsi6(ixOmin1:ixOmax1) * gamma(ixOmin1:ixOmax1,&
             idir, idir) * ( prim(ixOmin1:ixOmax1,&
             rho_) * h(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
             W_vel(idim)) * prim(ixOmin1:ixOmax1,&
             W_vel(idir))- prim(ixOmin1:ixOmax1,&
             Evec(idim)) * prim(ixOmin1:ixOmax1,&
             Evec(idir)) - prim(ixOmin1:ixOmax1,&
             Bvec(idim)) * prim(ixOmin1:ixOmax1, Bvec(idir)) ) 
       end do
       f(ixOmin1:ixOmax1, mom(idim)) = f(ixOmin1:ixOmax1,&
           mom(idim)) + alppsi6(ixOmin1:ixOmax1) * ( prim(ixOmin1:ixOmax1,&
           press_) + 0.5d0 * ( E2(ixOmin1:ixOmax1) + B2(ixOmin1:ixOmax1)) )
   
       ! Energy flux 
       f(ixOmin1:ixOmax1, tau_) = - prim(ixOmin1:ixOmax1,&
           beta(idim)) * cons(ixOmin1:ixOmax1, tau_) + prim(ixOmin1:ixOmax1,&
           alp_) * ( cons(ixOmin1:ixOmax1,&
           mom(idim)) * inv_gamma_ii(ixOmin1:ixOmax1) - cons(ixOmin1:ixOmax1,&
           D_) * prim(ixOmin1:ixOmax1, W_vel(idim)) / lfac(ixOmin1:ixOmax1) ) 
    else
       f(ixOmin1:ixOmax1, D_) = 0.0d0
       f(ixOmin1:ixOmax1, tau_) = 0.0d0
       f(ixOmin1:ixOmax1, mom(:)) = 0.0d0
    end if

    if ( evolve_EM ) then
       ! flux of E and B field
       do idir = 1, ndir
         if (idir /= idim) then
            f(ixOmin1:ixOmax1, Econs(idir)) = - prim(ixOmin1:ixOmax1,&
                beta(idim)) * cons(ixOmin1:ixOmax1,&
                Econs(idir)) + prim(ixOmin1:ixOmax1,&
                beta(idir)) * cons(ixOmin1:ixOmax1, Econs(idim))
            f(ixOmin1:ixOmax1, Bcons(idir)) = - prim(ixOmin1:ixOmax1,&
                beta(idim)) * cons(ixOmin1:ixOmax1,&
                Bcons(idir)) + prim(ixOmin1:ixOmax1,&
                beta(idir)) * cons(ixOmin1:ixOmax1, Bcons(idim))
         else
            f(ixOmin1:ixOmax1, Econs(idir)) = 0.0d0
            f(ixOmin1:ixOmax1, Bcons(idir)) = 0.0d0
         end if
       end do

       do idir = 1, ndir
          do jdir = 1, ndir
             f(ixOmin1:ixOmax1, Econs(idir)) = f(ixOmin1:ixOmax1,&
                 Econs(idir)) + prim(ixOmin1:ixOmax1,&
                 alp_) * inv_sqrt_gamma(ixOmin1:ixOmax1) * lvc(idim, idir,&
                 jdir) * gamma(ixOmin1:ixOmax1, jdir,&
                 jdir) * cons(ixOmin1:ixOmax1, Bcons(jdir))
             f(ixOmin1:ixOmax1, Bcons(idir)) = f(ixOmin1:ixOmax1,&
                 Bcons(idir)) - prim(ixOmin1:ixOmax1,&
                 alp_) * inv_sqrt_gamma(ixOmin1:ixOmax1) * lvc(idim, idir,&
                 jdir) * gamma(ixOmin1:ixOmax1, jdir,&
                 jdir) * cons(ixOmin1:ixOmax1, Econs(jdir))
          end do
       end do
    else
       f(ixOmin1:ixOmax1, Bcons(:)) = 0.0d0
       f(ixOmin1:ixOmax1, Econs(:)) = 0.0d0
    end if

    ! beware of sqrt_gamma = 0 cases, which means it hits the coordinate singularities!
    where( sqrt_gamma(ixOmin1:ixOmax1) < smalldouble )
       f(ixOmin1:ixOmax1, D_) = 0.0d0
       f(ixOmin1:ixOmax1, tau_) = 0.0d0
       f(ixOmin1:ixOmax1, mom(idim)) = 0.0d0
       f(ixOmin1:ixOmax1, Econs(idim)) = 0.0d0
       f(ixOmin1:ixOmax1, Bcons(idim)) = 0.0d0
    end where
  end subroutine gremhd_get_flux

end module mod_gremhd_phys_flux
