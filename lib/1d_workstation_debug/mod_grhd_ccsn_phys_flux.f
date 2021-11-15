module mod_grhd_ccsn_phys_flux
  use mod_physics
  use mod_grhd_ccsn_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grhd_ccsn_phys_flux_init

contains
  !> Initialize the module
  subroutine grhd_ccsn_phys_flux_init()
    phys_get_cbounds         => grhd_ccsn_get_cbounds
    phys_get_flux            => grhd_ccsn_get_flux
  end subroutine grhd_ccsn_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine grhd_ccsn_get_cbounds(consL, consR, primL, primR, x, ixImin1,&
     ixImax1, ixOmin1,ixOmax1, idim, cmax, cmin)
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

    call grhd_ccsn_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        primL(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lambdaL(ixImin1:ixImax1,1:2))
    call grhd_ccsn_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
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

  end subroutine grhd_ccsn_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grhd_ccsn_get_flux(cons, prim, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixImin1:ixImax1, 1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1, 1:ncons)

    double precision                :: v_hat(ixImin1:ixImax1) 
    double precision                :: lfac(ixImin1:ixImax1) 
    double precision                :: alppsi6(ixImin1:ixImax1) 
    integer                         :: idir

    call grhd_ccsn_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lfac=lfac(ixImin1:ixImax1) )

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

       ! electron fraction  flux
       f(ixOmin1:ixOmax1, ye_con_) = v_hat(ixOmin1:ixOmax1) * &
          cons(ixOmin1:ixOmax1, ye_con_)
   
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
       f(ixOmin1:ixOmax1, ye_con_) = 0.0d0
       f(ixOmin1:ixOmax1, mom(1:ndir)) = 0.0d0
       f(ixOmin1:ixOmax1, tau_) = 0.0d0
    end if

  end subroutine grhd_ccsn_get_flux

end module mod_grhd_ccsn_phys_flux
