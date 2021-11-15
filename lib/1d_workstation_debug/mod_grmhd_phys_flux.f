module mod_grmhd_phys_flux
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_flux_init

contains
  !> Initialize the module
  subroutine grmhd_phys_flux_init()
    phys_get_cbounds         => grmhd_get_cbounds
    phys_get_flux            => grmhd_get_flux
  end subroutine grmhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine grmhd_get_cbounds(consL, consR, primL, primR, x, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, idim, cmax, cmin)
    use mod_global_parameters

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

    double precision, dimension(ixImin1:ixImax1,1:2) :: lambdaL, lambdaR
    double precision, dimension(ixImin1:ixImax1,1:2) :: tmp_c

    call grmhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        primL(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lambdaL(ixImin1:ixImax1,1:2))
    call grmhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        primR(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lambdaR(ixImin1:ixImax1,1:2))

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

  end subroutine grmhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grmhd_get_flux(cons, prim, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      idim, f)
    use mod_global_parameters

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
    double precision                :: B_dot_v(ixImin1:ixImax1)
    double precision                :: Ptot(ixImin1:ixImax1) ! total pressure
    double precision                :: lfac(ixImin1:ixImax1) ! Lorentz factor
    double precision                :: B_over_W(ixImin1:ixImax1, 1:ndir) 
    integer                         :: idir

    call grmhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3), lfac=lfac(ixImin1:ixImax1),&
        v_hat=v_hat(ixImin1:ixImax1,1:ndir), B_dot_v=B_dot_v(ixImin1:ixImax1),&
        Ptot=Ptot(ixImin1:ixImax1) )

    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
       alp_) * prim(ixOmin1:ixOmax1,psi_)**6

    ! B^i / W
    do idir = 1, ndir
       B_over_W(ixOmin1:ixOmax1, idir) = prim(ixOmin1:ixOmax1,&
           Bvec(idir)) / lfac(ixOmin1:ixOmax1)
    end do

    ! Density flux
    f(ixOmin1:ixOmax1, D_) = v_hat(ixOmin1:ixOmax1,&
        idim) * cons(ixOmin1:ixOmax1, D_)

    ! Momentum flux f^idim_idir
    do idir = 1, ndir
      f(ixOmin1:ixOmax1, mom(idir)) = v_hat(ixOmin1:ixOmax1,&
          idim) * cons(ixOmin1:ixOmax1, mom(idir)) - alppsi6(ixOmin1:ixOmax1) &
         * ( B_over_W(ixOmin1:ixOmax1, idir) + B_dot_v(ixOmin1:ixOmax1) * &
         prim(ixOmin1:ixOmax1, W_vel(idir)) ) * gamma(ixOmin1:ixOmax1, idir,&
          idir) * B_over_W(ixOmin1:ixOmax1, idim)
    end do
    f(ixOmin1:ixOmax1, mom(idim)) = f(ixOmin1:ixOmax1,&
        mom(idim)) + alppsi6(ixOmin1:ixOmax1) * Ptot(ixOmin1:ixOmax1) 

    ! Energy flux 
    f(ixOmin1:ixOmax1, tau_) =  v_hat(ixOmin1:ixOmax1,&
        idim) * cons(ixOmin1:ixOmax1, tau_) + alppsi6(ixOmin1:ixOmax1) * ( &
       Ptot(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
        W_vel(idim)) / lfac(ixOmin1:ixOmax1) ) - alppsi6(ixOmin1:ixOmax1) * &
       B_dot_v(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1, Bvec(idim)) 

    ! flux of B field
    do idir = 1, ndir
      f(ixOmin1:ixOmax1, Bcons(idir)) = v_hat(ixOmin1:ixOmax1,&
          idim) * cons(ixOmin1:ixOmax1, Bcons(idir)) - v_hat(ixOmin1:ixOmax1,&
          idir) * cons(ixOmin1:ixOmax1, Bcons(idim))
    end do

    ! flux terms when using GLM 
    if (type_divb == divb_GLM) then
      f(ixOmin1:ixOmax1, Bphi_cons_) = prim(ixOmin1:ixOmax1,&
         alp_) * cons(ixOmin1:ixOmax1, Bcons(idim)) - alppsi6(ixOmin1:ixOmax1) &
         * prim(ixOmin1:ixOmax1,beta(idim)) * cons(ixOmin1:ixOmax1,&
          Bphi_cons_)
      do idir = 1, ndir
        f(ixOmin1:ixOmax1, Bcons(idir)) = f(ixOmin1:ixOmax1,&
            Bcons(idir)) + prim(ixOmin1:ixOmax1,&
            beta(idir)) * cons(ixOmin1:ixOmax1, Bcons(idim))
      end do
    end if

  end subroutine grmhd_get_flux

end module mod_grmhd_phys_flux
