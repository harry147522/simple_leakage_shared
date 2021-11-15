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
     ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, cmax,&
      cmin)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative left and right status
    double precision, intent(in)    :: consL(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons), consR(ixImin1:ixImax1,ixImin2:ixImax2, 1:ncons)
    ! primitive left and right status
    double precision, intent(in)    :: primL(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim), primR(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: tmp1,tmp2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:2) :: lambdaL, lambdaR

    call grhd_ccsn_get_lambda(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, primL(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim),&
        x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), lambdaL(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:2))
    call grhd_ccsn_get_lambda(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, primR(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim),&
        x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), lambdaR(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:2))
  
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(0.0d0, lambdaL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1), lambdaR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) )
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(0.0d0, lambdaL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2), lambdaR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) )

    if(present(cmin)) then
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
      cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    else
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(abs(tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)), abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
    end if

  end subroutine grhd_ccsn_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grhd_ccsn_get_flux(cons, prim, x, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)

    double precision                :: v_hat(ixImin1:ixImax1,ixImin2:ixImax2) 
    double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2) 
    double precision                :: alppsi6(ixImin1:ixImax1,&
       ixImin2:ixImax2) 
    integer                         :: idir

    call grhd_ccsn_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        lfac=lfac(ixImin1:ixImax1,ixImin2:ixImax2) )

    ! v_hat = v * alp - beta
    v_hat(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,W_vel(idim)) / lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       alp_) - prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,beta(idim))

    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,alp_) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)**6

    if ( evolve_hydro ) then 
       ! Density flux
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = v_hat(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_)

       ! electron fraction  flux
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, ye_con_) = v_hat(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, ye_con_)
   
       ! Momentum flux
       do idir = 1, ndir
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = v_hat(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) * cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idir))
       end do
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idim)) = f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idim)) + alppsi6(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, press_) 
   
       ! Energy flux 
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = v_hat(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           tau_) + alppsi6(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           press_) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          W_vel(idim)) / lfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    else
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = 0.0d0
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, ye_con_) = 0.0d0
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(1:ndir)) = 0.0d0
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = 0.0d0
    end if

  end subroutine grhd_ccsn_get_flux

end module mod_grhd_ccsn_phys_flux
