module mod_gremhd_phys_convert
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: gremhd_phys_convert_init

contains

  !> Initialize the module
  subroutine gremhd_phys_convert_init()
    ! initialize the k_max and v_max
    v_max = dsqrt(1.0d0 - 1.0d0 / lfac_max**2)
    k_max = 2.0d0 * v_max / (1.0d0 + v_max**2) ! k_max < 1
    phys_to_conserved        => gremhd_to_conserved
    phys_to_primitive        => gremhd_to_primitive
  end subroutine gremhd_phys_convert_init

  !> Transform primitive variables into conservative ones
  subroutine gremhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, cons, prim, x)
    use mod_global_parameters
    use mod_eos, only: small_tau
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: cons(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)
    double precision, intent(inout) :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    double precision                :: h(ixImin1:ixImax1,ixImin2:ixImax2) !enthalpy
    double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2) !Lorentz factor
    double precision                :: E2(ixImin1:ixImax1,ixImin2:ixImax2),&
        B2(ixImin1:ixImax1,ixImin2:ixImax2) 
    double precision                :: psi6(ixImin1:ixImax1,ixImin2:ixImax2) !conformal factor **6
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    double precision                :: sqrt_gamma(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                         :: idir, jdir, kdir

    if ( fix_small_values ) then
       call gremhd_handle_small_values(prim(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, .True.,&
           'gremhd_to_conserved')
    end if
    psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, psi_)**6 

    call gremhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        sqrt_gamma=sqrt_gamma(ixImin1:ixImax1,ixImin2:ixImax2),&
        E2=E2(ixImin1:ixImax1,ixImin2:ixImax2), B2=B2(ixImin1:ixImax1,&
       ixImin2:ixImax2), lfac=lfac(ixImin1:ixImax1,ixImin2:ixImax2),&
        h=h(ixImin1:ixImax1,ixImin2:ixImax2) )

    ! Conserved density D
    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) 

    ! Convert velocity to covariant momentum
    do idir = 1, ndir
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) = ( cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           D_) * h(ixOmin1:ixOmax1,ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, W_vel(idir)) ) * gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,idir,idir)
       do jdir = 1, ndir
          do kdir = 1, ndir
             cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) + sqrt_gamma(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) * lvc(idir,jdir,kdir) * prim(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2, Evec(jdir)) * prim(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2, Bvec(kdir))
          end do
       end do
    end do

    ! Conserved energy - D = tau
    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, D_) * ( lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * h(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - 1.0d0 ) - prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        press_) + 0.5d0 * ( E2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )
    where ( cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) < small_tau )
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = small_tau
    end where

    ! conformal transformation
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
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Econs(idir)) = prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Evec(idir)) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do
    if (type_divb == divb_GLM)  cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        Bphi_cons_) = prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        Bphi_) * psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
  end subroutine gremhd_to_conserved

  !> Transform conservative variables into primitive ones 
  subroutine gremhd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
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

    integer                         :: idir, jdir, kdir
    integer                         :: flag(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: psi6(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: cons_tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ncons)
    integer                         :: ix1,ix2 
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3) !metric gamma_ij
    double precision                :: sqrt_gamma(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                :: sqrt_gamma_hat(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    ! variables needed for the root founding
    logical                         :: adjustment
    integer                         :: error_code
    ! temp conserved variables (at one grid)
    double precision                :: D, S, tau
    double precision                :: r, q, k

    double precision                :: v_bar ! temp velocity
    double precision                :: lfac ! Lorentz factor
    double precision                :: ui_bar(1:ndir)
    double precision                :: rho_bar, eps_bar
    double precision                :: a_bar, h_bar
    double precision                :: eps_min, eps_max
    double precision                :: rescale_factor
    double precision                :: z ! this is the root of the equation
    double precision                :: zp, zm

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
    call get_sqrt_gamma_hat(x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
        sqrt_gamma_hat(ixImin1:ixImax1,ixImin2:ixImax2))
    sqrt_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = psi6(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * sqrt_gamma_hat(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    ! Here, we assume that both E and B are known 
    do idir = 1, ndir
       prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bcons(idir)) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Evec(idir)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Econs(idir)) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! conformal transformation back to normal conserved variables
    cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2, D_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, D_) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, tau_) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    do idir = 1, ndir
       cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) / psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! remove EM part in the conserved variables
    do idir = 1, ndir
       cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           tau_) = cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           tau_) - 0.5d0 * gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,&
           idir) * ( prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Evec(idir))**2 + prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir))**2 )
    end do
    do idir = 1, ndir
       do jdir = 1, ndir
          do kdir = 1, ndir
             cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) = cons_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) - sqrt_gamma(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) * lvc(idir, jdir,&
                 kdir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 Evec(jdir)) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 Bvec(kdir))
          end do
       end do
    end do

    !Note that v is contravariant
    do ix1 = ixOmin1,ixOmax1 
    do ix2 = ixOmin2,ixOmax2 
      adjustment = .False.
      if ( (cons_tmp(ix1,ix2, D_) <= small_D) .or. (cons_tmp(ix1,ix2,&
          tau_) <= small_tau) ) then
         ! atmosphere handling
         ! skip the primitive recovery
         lfac = 1.0d0
         ui_bar  = 0.0d0
         rho_bar = small_rho
         eps_bar = small_eps
         adjustment = .True. ! adjust all conserved variables
      else
         D = cons_tmp(ix1,ix2, D_)
         S = 0.0d0
         do idir = 1, ndir
            S = S + cons_tmp(ix1,ix2, mom(idir))**2 / gamma(ix1,ix2,idir,idir)
         end do
         S = dsqrt( S )
         tau = cons_tmp(ix1,ix2, tau_)
   
         r = S / D
         q = max( tau / D, 0.0d0 )  ! limit q to be non-negative
         k = S / (tau + D) ! k = r / (1+q)
   
         ! find the root z
         if (dabs(k) <= smalldouble) then
            z = 0.0d0
         else
            if (k > k_max) then 
               ! rescale the momentum such that r=k_max*(1+q)
               rescale_factor = k_max / k
               r = k_max * (1.0d0 + q)
               k = k_max
               cons_tmp(ix1,ix2, mom(:)) = cons_tmp(ix1,ix2,&
                   mom(:)) * rescale_factor
               adjustment = .True. ! adjust all conserved variables
            endif
          
            ! range of the root
            zm = 0.5d0 * k / dsqrt(1.0d0 - k**2 / 4.0d0)
            zp = k / dsqrt(1.0d0 - k**2)
            if ( dabs(func(zm)) <= tolerance ) then
               z = zm
            else if ( dabs(func(zp)) <= tolerance ) then
               z = zp
            else
               call rootfinding_illinois(z, zm, zp, tolerance, iter_max,&
                   error_code, func)
               ! after getting z
               ! check if there are any errors or correction
               select case (error_code)
               !case (0) ! root is found
               !   write(*,*) "z= ", z
               case (-1) ! nothing happened
                  call mpistop&
                     ("have you ever attemp to find the root in con2prim?")
               case (1) ! failed to find the root
                  call mpistop("Fail to find the root in con2prim")
               case (2) ! z is NaN
                  call mpistop("NaN")
                  flag(ix1,ix2) = -1 !This special flag is used to indicate NaN error. Normally it wont be negative
               end select
            endif 
         endif 
   
         ! calculate all the primitive variables from the solution z
         call get_vars(z, W_out=lfac, rho_out=rho_bar, eps_out=eps_bar,&
             h_out=h_bar)
         do idir=1,ndir
            ui_bar(idir) = cons_tmp(ix1,ix2, mom(idir)) / gamma(ix1,ix2,idir,&
               idir) / ( h_bar * D ) 
         enddo
   
         ! check if we need any adjustments
         if ( rho_bar <= small_rho_thr ) then
            ! reset the primitive variables
            lfac = 1.0d0
            ui_bar = 0.0d0
            rho_bar = small_rho
            eps_bar = small_eps
            adjustment = .True. ! adjust conserved variables
         else
            ! check if eps fails into the validity range
            call eos_get_eps_range(rho_bar, eps_min, eps_max)
            if ( eps_bar < eps_min .or. eps_bar > eps_max ) then
               eps_bar = max( min( eps_max, eps_bar ), eps_min )
               adjustment = .True. ! adjust conserved variables
            end if
      
            ! limit the velocities
            v_bar = z / lfac
            if (v_bar > v_max) then 
               ! rescale the velocity such that v_bar = v_max and keeping D constant
               rescale_factor = ( v_max / v_bar ) * ( lfac_max / lfac )
               ui_bar = ui_bar * rescale_factor
               v_bar = v_max
               lfac = lfac_max
      
               ! although D is kept constant, density is changed a bit
               rho_bar = D / lfac
               ! limit eps again based on the updated rho
               call eos_get_eps_range(rho_bar, eps_min, eps_max)
               eps_bar = max( min( eps_max, eps_bar ), eps_min )
               adjustment = .True. ! adjust conserved variables
            end if
         end if
      end if
   
      ! update all the hydro primitive variables here
      prim(ix1,ix2, rho_) = rho_bar
      prim(ix1,ix2, eps_) = eps_bar
      do idir=1,ndir
         prim(ix1,ix2, W_vel(idir)) = ui_bar(idir)
      enddo
      call gremhd_update_eos_one_point(prim(ix1,ix2, 1:nprim))       

      ! if any adjustments, recompute conserved varables consistently
      if ( adjustment ) then
         h_bar = 1.0d0 + eps_bar + prim(ix1,ix2, press_) / rho_bar
         cons(ix1,ix2, D_) = lfac * rho_bar
         ! Convert velocity to covariant momentum
         do idir = 1, ndir
            cons(ix1,ix2, mom(idir)) = ( cons(ix1,ix2, D_) * h_bar * prim(ix1,&
               ix2, W_vel(idir)) ) * gamma(ix1,ix2,idir,idir)
            do jdir = 1, ndir
               do kdir = 1, ndir
                  cons(ix1,ix2, mom(idir)) = cons(ix1,ix2,&
                      mom(idir)) + sqrt_gamma(ix1,ix2) * lvc(idir,jdir,&
                     kdir) * prim(ix1,ix2, Evec(jdir)) * prim(ix1,ix2,&
                      Bvec(kdir))
               end do
            end do
         end do
     
         ! Conserved energy - D = tau
         cons(ix1,ix2, tau_) = cons(ix1,ix2,&
             D_) * ( lfac * h_bar - 1.0d0 ) - prim(ix1,ix2, press_)
         do idir = 1, ndir
            cons(ix1,ix2, tau_) = cons(ix1,ix2, tau_) + 0.5d0 * ( prim(ix1,ix2,&
                Evec(idir))**2 + prim(ix1,ix2, Bvec(idir))**2 ) * gamma(ix1,&
               ix2, idir, idir)
         end do
     
         ! conformal transformation
         cons(ix1,ix2, D_) = cons(ix1,ix2, D_) * psi6(ix1,ix2)
         cons(ix1,ix2, tau_) = cons(ix1,ix2, tau_) * psi6(ix1,ix2)
         do idir = 1, ndir
            cons(ix1,ix2, mom(idir))   = cons(ix1,ix2, mom(idir)) * psi6(ix1,&
               ix2)
         end do
      end if ! end adjustments for conserved

    end do
    end do

    !check if NaN
    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < 0)) then
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
       call small_values_error(prim, x, ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2, flag,&
           'gremhd_to_primitive: (NaN z)')
    end if

    contains
       subroutine get_vars(z, W_out, rho_out, eps_out, h_out)
          implicit none
          double precision, intent(in)    :: z
          double precision, intent(out), optional   :: W_out
          double precision, intent(out), optional   :: rho_out
          double precision, intent(out), optional   :: eps_out
          double precision, intent(out), optional   :: h_out

          double precision                :: W_bar, eps_bar, a_bar, rho_bar,&
              prs_bar
          W_bar = dsqrt(1.0d0+z**2)
          if (present(W_out)) W_out = W_bar
          rho_bar = D / W_bar
          rho_bar = max( min( eos_rhomax, rho_bar ), eos_rhomin )
          if (present(rho_out)) rho_out = rho_bar
          eps_bar = W_bar * q - z * r + z**2/(1.0d0+W_bar)
          call eos_get_eps_range(rho_bar,eps_min,eps_max)
          eps_bar = max( min( eps_max, eps_bar ), eps_min )
          if (present(eps_out)) eps_out = eps_bar
          if (present(h_out)) then
             call eos_get_pressure_one_grid(prs_bar, rho_bar, eps_bar)
             a_bar = prs_bar / (rho_bar*(1.0d0+eps_bar))
             a_bar = max( min( 1.0d0, a_bar ) , 0.0d0)
             h_out = (1.0d0+eps_bar)*(1.0d0+a_bar); !h_bar = max( hbar, 0.0d0 )
          end if
       end subroutine get_vars          

       !> master function f(z) for finding the root z
       double precision function func(z)
          double precision, intent(in)    :: z
          double precision                :: h_bar
          call get_vars(z, h_out=h_bar)
          func = z - r / h_bar
       end function func
  end subroutine gremhd_to_primitive

end module mod_gremhd_phys_convert
