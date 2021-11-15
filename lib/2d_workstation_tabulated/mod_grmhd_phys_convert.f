module mod_grmhd_phys_convert
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  double precision                ::  h0 !lower bound for the relativistic enthalpy over entire validity region of eos
  double precision                ::  one_over_h0  ! 1 / h0

  ! Public methods
  public :: grmhd_phys_convert_init

contains

  !> Initialize the module
  subroutine grmhd_phys_convert_init()
    use mod_global_parameters, only: bigdouble, smalldouble
    use mod_eos, only: eos_hmin

    ! initialize 1/h0
    h0 = eos_hmin
    one_over_h0 = 1.0d0 / h0

    ! initialize the k_max and v_max
    v_max = dsqrt(1.0d0 - 1.0d0 / lfac_max**2)

    phys_to_conserved        => grmhd_to_conserved
    phys_to_primitive        => grmhd_to_primitive
  end subroutine grmhd_phys_convert_init

  !> Transform primitive variables into conservative ones
  subroutine grmhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, cons, prim, x)
    use mod_global_parameters
    use mod_eos, only: small_rho_thr, small_rho, small_eps
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: cons(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)
    double precision, intent(inout) :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    double precision                :: B_dot_v(ixImin1:ixImax1,&
       ixImin2:ixImax2)   
    double precision                :: htot(ixImin1:ixImax1,ixImin2:ixImax2) !modified enthalpy:(h + b2/rho)
    double precision                :: Ptot(ixImin1:ixImax1,ixImin2:ixImax2) !total pressure
    double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2) !Lorentz factor
    double precision                :: psi6(ixImin1:ixImax1,ixImin2:ixImax2) !conformal factor **6
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    integer                         :: idir

    call grmhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        lfac=lfac(ixImin1:ixImax1,ixImin2:ixImax2),&
        B_dot_v=B_dot_v(ixImin1:ixImax1,ixImin2:ixImax2),&
        Ptot=Ptot(ixImin1:ixImax1,ixImin2:ixImax2), htot=htot(ixImin1:ixImax1,&
       ixImin2:ixImax2) )

    ! Conserved density D
    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) 

    ! Convert velocity to covariant momentum
    do idir = 1, ndir
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) = ( ( cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           D_) * htot(ixOmin1:ixOmax1,ixOmin2:ixOmax2) - lfac(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * B_dot_v(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)**2 ) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           W_vel(idir)) - B_dot_v(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir)) ) * gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,idir)
    end do

    ! Conserved energy - D = tau
    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, D_) * ( lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * htot(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - 1.0d0 ) - Ptot(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - ( lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * B_dot_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )**2

    ! conformal transformation
    psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, psi_)**6 
    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, D_) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, tau_) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    do idir = 1, ndir
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))   = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bcons(idir)) = prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir)) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do
    if (type_divb == divb_GLM)  cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        Bphi_cons_) = prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        Bphi_) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
  end subroutine grmhd_to_conserved

  !> Transform conservative variables into primitive ones 
  subroutine grmhd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, cons, prim, x)
    use mod_global_parameters
    use mod_geometry
    use mod_small_values
    use mod_eos
    use mod_rootfinding
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: cons(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)
    double precision, intent(inout) :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    integer                         :: idir, flag(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                :: psi6(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: cons_tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ncons)
    integer                         :: ix1,ix2 
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3) !metric gamma_ij_hat
    ! rescaled variables
    double precision                :: q
    double precision                :: r_i(1:ndir) !r_i
    double precision                :: bi(1:ndir) !bi,bi, note that this is a rescaled Bfield, do not mix with the one in prim2con
    double precision                :: r_dot_b, b2, r2
    double precision                :: b_sqr_r_norm_sqr
    double precision                :: vi_hat(1:ndir)
    double precision                :: rescale_factor
    ! variables needed during the root founding
    logical                         :: adjustment
    integer                         :: error_code
    double precision                :: eps_min, eps_max
    double precision                :: W_hat, eps_hat, a_hat, rho_hat, prs_hat
    double precision                :: v_hat_sqr, v0_sqr
    double precision                :: chi
    double precision                :: r_bar_sqr
    double precision                :: q_bar
    double precision                :: mu_bounds(2) ! bounds of mu_plus
    double precision                :: mu_plus ! upper bound of mu
    double precision                :: mu ! this is the root of the equation
    ! extra variables needed after the root founding
    double precision                :: B_dot_v
    double precision                :: h_hat

    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0
    call get_gamma_ij_hat(x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), ixImin1,&
       ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
        gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3))
    do idir = 1, ndir
       gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) = gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, psi_)**4 
    end do
    psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, psi_)**6 

    ! conformal transformation back to normal conserved variables
    cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, D_) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, tau_) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    do idir = 1, ndir
       cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       ! Nothing to do with B field actually
       prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bcons(idir)) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do
    if (type_divb == divb_GLM) prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        Bphi_) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        Bphi_cons_) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !Note that v is contravariant
    do ix1 = ixOmin1,ixOmax1 
    do ix2 = ixOmin2,ixOmax2 
      adjustment = .False.
      if ( (cons_tmp(ix1,ix2, D_) <= small_D) .or. (cons_tmp(ix1,ix2,&
          tau_) <= small_tau) ) then
         ! atmosphere handling 
         ! skip the primitive recovery
         W_hat = 1.0d0
         vi_hat(1:ndir) = 0.0d0
         rho_hat = small_rho
         eps_hat = small_eps
         adjustment = .True. ! adjust all conserved variables
      else 
         ! Step 0: initialize all the variables
         ! rescale the conserved variables
         q = cons_tmp(ix1,ix2, tau_) / cons_tmp(ix1,ix2, D_)
         do idir = 1, ndir
            r_i(idir) = cons_tmp(ix1,ix2, mom(idir)) / cons_tmp(ix1,ix2, D_)
            bi(idir) = prim(ix1,ix2, Bvec(idir)) / dsqrt( cons_tmp(ix1,ix2,&
                D_) )
         end do
   
         ! some useful relations
         r_dot_b = 0.0d0
         b2 = 0.0d0 
         r2 = 0.0d0
         do idir = 1, ndir
            r_dot_b = r_dot_b + r_i(idir) * bi(idir)
            b2 = b2 + bi(idir)**2 * gamma(ix1,ix2,idir,idir)
            r2 = r2 + r_i(idir)**2 / gamma(ix1,ix2,idir,idir)
         end do
           
         v0_sqr = r2 / ( h0**2 + r2 )
   
         ! decompose the momentum in parts parallel and normal to the magnetic
         ! field, and we need b**2 * r_norm**2 only:
         b_sqr_r_norm_sqr = b2 * r2 - r_dot_b**2
   
         ! Step 1: solve mu+
         ! note that func_for_mu_plus(mu) <= 0 is always true.
         if ( dabs(func_for_mu_plus(one_over_h0)) <= tolerance ) then
            mu_plus = one_over_h0
         else
            ! try constrained newton raphson first
            mu_bounds(1) = 0.0d0
            mu_bounds(2) = one_over_h0 
            call rootfinding_constrained_newton_raphson(mu_plus, mu_bounds(1),&
                mu_bounds(2), tolerance, 50, error_code, func_for_mu_plus,&
                d_mu_plus)

            if ( error_code /= 0 ) then
               !if newton raphson failed, try illinois
               write(*,*)&
"Warning: Newton raphson failed when finding mu+, now try illinois"
               if (any(mu_bounds(:) /= mu_bounds(:))) then
                  mu_bounds(1) = 0.0d0
                  mu_bounds(2) = one_over_h0 
               end if
               call rootfinding_illinois(mu_plus, mu_bounds(1), mu_bounds(2),&
                   tolerance, iter_max, error_code, func_for_mu_plus)
            end if
            ! check mu+
            select case (error_code)
            !case (0) ! root is found
            !   write(*,*) "z= ", z
            case (-1) ! nothing happened
               call mpistop("have you ever attemp to find mu+ in con2prim?")
            case (1) ! failed to find the root
               call mpistop("Fail to find the root for mu_plus")
               !write(*,*) "Warning: Fail to find the root for mu_plus, now use the maximum value"
               !mu_plus = one_over_h0
            case (2) ! root is NaN
               call mpistop("NaN mu_plus")
               !write(*,*) "Warning: mu_plus is NaN, now use the maximum value"
               !mu_plus = one_over_h0
            case (3) ! z is not bracketed
               call mpistop("the root is not bracketed for mu_plus")
            end select
         end if
   
         ! Step 2: solve mu
         ! check if the bound make sense, use mu_bounds(1) as a tmp variable
         mu_bounds(1) = func_for_mu(mu_plus)
         if ( dabs(mu_bounds(1)) <= tolerance ) then
            mu = mu_plus
         else
            if ( mu_bounds(1) < 0.0d0 ) mu_plus = one_over_h0
            call rootfinding_illinois(mu, 0.0d0, mu_plus, tolerance, iter_max,&
                error_code, func_for_mu)
            ! check mu
            select case (error_code)
            !case (0) ! root is found
            !   write(*,*) "z= ", z
            case (-1) ! nothing happened
               call mpistop("have you ever attemp to find mu in con2prim?")
            case (1) ! failed to find the root
               call mpistop("Fail to find the root for mu")
            case (2) ! root is NaN
               ! This special flag is used to indicate NaN error. 
               ! Normally it wont be negative
               flag(ix1,ix2) = -1 
            case (3) ! z is not bracketed
               call mpistop("the root is not bracketed for mu")
            end select
         end if
   
         ! Step 3: calculate all the primitive variables from the solution mu
         chi = chi_of_mu(mu)
         r_bar_sqr = r2 * chi**2 + mu * chi * ( 1.0d0 + chi ) * r_dot_b**2 
         q_bar = q - 0.5d0 * b2 - 0.5d0 * mu**2 * chi**2 * b_sqr_r_norm_sqr
         v_hat_sqr = min( mu**2 * r_bar_sqr , v0_sqr)
         W_hat = 1.0d0 / dsqrt( 1.0d0 - v_hat_sqr )
         rho_hat = cons_tmp(ix1,ix2, D_) / W_hat
         rho_hat = max( min( eos_rhomax, rho_hat ), eos_rhomin )
         eps_hat = W_hat * ( q_bar - mu * r_bar_sqr ) + v_hat_sqr * W_hat**2 / &
            ( 1.0d0 + W_hat )
         ! calculate velocities
         do idir = 1, ndir
            vi_hat(idir) = mu * chi * ( r_i(idir) / gamma(ix1,ix2,idir,&
               idir) + mu * r_dot_b * bi(idir) )
         end do
         
         ! adjust the results if it is invalid
         if ( rho_hat <= small_rho_thr ) then
            ! reset the primitive variables
            W_hat = 1.0d0
            vi_hat(1:ndir) = 0.0d0
            rho_hat = small_rho
            eps_hat = small_eps
            adjustment = .True. ! adjust conserved variables
         else
            ! limit the velocities
            if ( W_hat > lfac_max ) then
               ! rescale the velocity such that v = v_max and keeping D constant
               rescale_factor = v_max / dsqrt( v_hat_sqr )
               vi_hat(1:ndir) = vi_hat(1:ndir) * rescale_factor
               v_hat_sqr = v_max**2
               W_hat = lfac_max
               ! although D is kept constant, density is changed a bit
               rho_hat = cons_tmp(ix1,ix2, D_) / W_hat
               adjustment = .True. ! adjust conserved variables
            end if

            ! check if eps fails into the validity range
            call eos_get_eps_range(rho_hat, eps_min, eps_max)
            if ( eps_hat < eps_min ) then
               eps_hat = max( eps_hat, eps_min )
               ! This special flag is used to indicate that 
               ! we will only update tau without changing D and S
               flag(ix1,ix2) = 1
               adjustment = .True. ! adjust conserved variables
            else if ( eps_hat > eps_max ) then
               eps_hat = min( eps_max, eps_hat )
               adjustment = .True. ! adjust all conserved variables
            end if

         end if
      end if ! end of con2prim
   
      ! update the primitive variables,
      ! note that B field was updated at the begining of the loop
      call eos_get_pressure_one_grid(prs_hat, rho_hat, eps_hat)
      prim(ix1,ix2, rho_)          = rho_hat
      prim(ix1,ix2, eps_)          = eps_hat
      prim(ix1,ix2, press_)        = prs_hat
      prim(ix1,ix2, W_vel(1:ndir)) = vi_hat(1:ndir) * W_hat
      ! update cs2
      call eos_get_cs2_one_grid(prim(ix1,ix2,cs2_),prim(ix1,ix2, rho_),&
         prim(ix1,ix2, eps_))
      ! strictly require cs2 is physical
      !prim(ix^D,cs2_) = max( min( 1.0d0, prim(ix^D,cs2_) ) , 0.0d0)

      ! if any adjustments, recompute pressure and conserved varables consistently
      if ( adjustment ) then
         B_dot_v = 0.0d0
         do idir = 1, ndir
            B_dot_v = B_dot_v + prim(ix1,ix2,&
                Bvec(idir)) * vi_hat(idir) * gamma(ix1,ix2,idir,idir)
         end do
         b2 = B_dot_v**2
         do idir = 1, ndir
            b2 = b2 + prim(ix1,ix2, Bvec(idir))**2 * gamma(ix1,ix2,idir,&
               idir) / W_hat**2
         end do
         h_hat = 1.0d0 + eps_hat + ( prs_hat + b2 ) / rho_hat !magnetrically modified enthalpy
         prs_hat = prs_hat + 0.5d0 * b2 ! magnetrically modified pressure

         ! update conserved variables
         cons_tmp(ix1,ix2, D_) = W_hat * rho_hat
         cons_tmp(ix1,ix2, tau_) = cons_tmp(ix1,ix2,&
             D_) * ( W_hat * h_hat - 1.0d0 ) - prs_hat - ( W_hat * B_dot_v &
            )**2

         ! conformal transformation
         cons(ix1,ix2, tau_) = psi6(ix1,ix2) * cons_tmp(ix1,ix2, tau_) 

         ! we need to update all conserved variables
         if ( flag(ix1,ix2) /= 1 ) then
            do idir = 1, ndir
               cons(ix1,ix2, mom(idir)) = ( ( cons(ix1,ix2,&
                   D_) * h_hat - W_hat * B_dot_v**2 ) * prim(ix1,ix2,&
                   W_vel(idir)) - B_dot_v * prim(ix1,ix2,&
                   Bvec(idir)) ) * gamma(ix1,ix2,idir,idir)
            end do
            ! conformal transformation
            cons(ix1,ix2, D_) = psi6(ix1,ix2) * cons_tmp(ix1,ix2, D_) 
            cons(ix1,ix2, mom(1:ndir)) = psi6(ix1,ix2) * cons_tmp(ix1,ix2,&
                mom(1:ndir)) 
         end if
      end if ! end adjustments for conserved

    enddo
    enddo

    !check if NaN
    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < 0)) then
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
       call small_values_error(prim, x, ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2, flag,&
           'grmhd_to_primitive: (NaN mu)')
    end if

    contains

       double precision function chi_of_mu(mu)
          double precision, intent(in)    :: mu
          chi_of_mu = 1.0d0 / (1.0d0 + mu * b2)
       end function chi_of_mu

       double precision function func_for_mu_plus(mu)
          double precision, intent(in)    :: mu
          double precision                :: chi, r_bar_sqr
          chi = chi_of_mu(mu)
          r_bar_sqr = r2 * chi**2 + mu * chi * ( 1.0d0 + chi ) * r_dot_b**2 
          func_for_mu_plus = mu * dsqrt( h0**2 + r_bar_sqr ) - 1.0d0
          !if (func_for_mu_plus /= func_for_mu_plus ) then
          !   write(*,*) mu, chi, r_bar_sqr
          !   write(*,*) h0**2 + r_bar_sqr
          !   write(*,*) "f_mu_plus NaN!!"
          !end if
       end function func_for_mu_plus

       double precision function d_mu_plus(mu)
          double precision, intent(in)    :: mu
          double precision                :: chi, d_chi
          double precision                :: r_bar_sqr, d_r_bar_sqr
          chi = chi_of_mu(mu)
          d_chi = - b2 / (1.0d0 + mu * b2)**2
          r_bar_sqr = r2 * chi**2 + mu * chi * ( 1.0d0 + chi ) * r_dot_b**2 
          d_r_bar_sqr = 2.0d0 * r2 * chi * d_chi + ( chi + mu * d_chi + chi**2 &
             + 2.0d0 * mu * chi * d_chi ) * r_dot_b**2 
          d_mu_plus = dsqrt( h0**2 + r_bar_sqr )
          d_mu_plus = d_mu_plus + mu * d_r_bar_sqr / d_mu_plus
       end function d_mu_plus

       !> master function f(mu) for finding the root mu
       double precision function func_for_mu(mu)
          double precision, intent(in)    :: mu
          double precision                :: chi
          double precision                :: r_bar_sqr
          double precision                :: q_bar
          double precision                :: W_hat, eps_hat, a_hat, rho_hat,&
              prs_hat
          double precision                :: v_hat_sqr
          double precision                :: mu_hat
          double precision                :: nu_A, nu_B, nu_hat
          chi = chi_of_mu(mu)
          r_bar_sqr = r2 * chi**2 + mu * chi * ( 1.0d0 + chi ) * r_dot_b**2
          q_bar = q - 0.5d0 * b2 - 0.5d0 * mu**2 * chi**2 * b_sqr_r_norm_sqr
          v_hat_sqr = min( mu**2 * r_bar_sqr , v0_sqr)
          W_hat = 1.0d0 / dsqrt( 1.0d0 - v_hat_sqr )
          rho_hat = cons_tmp(ix1,ix2, D_) / W_hat
          rho_hat = max( min( eos_rhomax, rho_hat ), eos_rhomin )
          eps_hat = W_hat * ( q_bar - mu * r_bar_sqr ) + v_hat_sqr * W_hat**2 &
             / ( 1.0d0 + W_hat )
          call eos_get_eps_range(rho_hat,eps_min,eps_max)
          eps_hat = max( min( eps_max, eps_hat ), eps_min )
          call eos_get_pressure_one_grid(prs_hat, rho_hat, eps_hat)
          a_hat = prs_hat / (rho_hat*(1.0d0+eps_hat))
          
          nu_A = (1.0d0 + a_hat) * (1.0d0+eps_hat) / W_hat ! h/W
          nu_B = (1.0d0 + a_hat) * ( 1.0d0 + q_bar - mu * r_bar_sqr ) 
          nu_hat = max(nu_A,nu_B)

          mu_hat = 1.0d0 / ( nu_hat + r_bar_sqr * mu )
          func_for_mu = mu - mu_hat
          !if (func_for_mu /= func_for_mu ) then
          !   write(*,*) mu, chi, r_bar_sqr
          !   write(*,*) nu_hat
          !   write(*,*) "f_mu NaN!!"
          !end if
       end function func_for_mu
  end subroutine grmhd_to_primitive

end module mod_grmhd_phys_convert
