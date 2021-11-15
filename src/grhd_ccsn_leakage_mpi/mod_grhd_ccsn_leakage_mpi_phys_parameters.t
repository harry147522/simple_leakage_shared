module mod_grhd_ccsn_leakage_mpi_phys_parameters
  use mod_physics
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: evolve_hydro = .True.
  logical                      :: use_process_adv_global = .false.
  logical                      :: use_GR = .false.

  !-------------------------------------------------------------------!
  ! Parameters for con2prim
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tolerance = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_max = 50000
  !> maximum Lorentz factor
 ! double precision             :: lfac_max = 1.1d1
  double precision             :: lfac_max = 1.6d1
  double precision             :: v_max, k_max

  !-------------------------------------------------------------------!
  ! Parameters for 1D core treatment
  !-------------------------------------------------------------------!
  logical                      :: oneDcore = .False.
  double precision             :: r_core = 0.0d0

  !-------------------------------------------------------------------!

  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grhd_ccsn_leakage_mpi_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)

    ! rho and eps are given, update the rest of the primitive variables
!    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_), ye=prim( ye_))
!    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_), ye=prim( ye_))

    call eos_eps_get_all_one_grid(rho = prim(rho_), temp = prim(temp_), ye = prim(ye_), &
                              eps = prim(eps_), ent = prim(ent_), cs2 = prim(cs2_) , prs = prim(press_))

       if ( prim(temp_) < eos_tempmin) then
           prim(temp_) = eos_tempmin*1.05d0
       endif

    ! strictly require cs2 is physical
    if ( prim(cs2_) >= 1.0d0 ) then
       prim(cs2_) = 0.0d0
    else
       prim(cs2_) = max( prim(cs2_) , 0.0d0)
    end if
  end subroutine grhd_ccsn_leakage_mpi_update_eos_one_point

  !> get some useful variables from primitive
  subroutine grhd_ccsn_leakage_mpi_get_intermediate_variables(ixI^L, ixO^L, prim, x, &
             gamma, lfac2, lfac, v2, h )
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: prim(ixI^S, 1:nprim)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(out), optional :: gamma(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: v2(ixI^S)        ! 
    double precision, intent(out), optional :: lfac2(ixI^S) ! Lorentz factor square
    double precision, intent(out), optional :: lfac(ixI^S) ! Lorentz factor
    double precision, intent(out), optional :: h(ixI^S) ! enthalpy: h

    integer                                 :: idir
    double precision                        :: W2v2(ixI^S)        ! 
    double precision                        :: v_bar(ixI^S,1:3) 
    double precision                        :: lfac2_tmp(ixI^S) ! Lorentz factor square
    double precision                        :: lfac_tmp(ixI^S) ! Lorentz factor
    double precision                        :: gamma_tmp(ixI^S,1:3,1:3)

    if ( present(gamma) .or. present(lfac2) .or. present(lfac) .or. present(v2) ) then
       ! get the metric
       call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
       do idir = 1, ndir
          gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
       end do
       if ( present(gamma) ) then
          gamma(ixO^S,1:3,1:3) = gamma_tmp(ixO^S,1:3,1:3) 
       end if
   
       if ( present(lfac2) .or. present(lfac) .or. present(v2) ) then
          W2v2 = 0.0d0
          ! calculate W^2 * v^2 first
          do idir = 1, ndir
             W2v2(ixO^S) = W2v2(ixO^S) + gamma_tmp(ixO^S,idir,idir)*prim(ixO^S, W_vel(idir))**2
          end do
          ! Calculate the Lorentz factor from velocity
          lfac2_tmp(ixO^S) = W2v2(ixO^S) + 1.0d0 
          if ( present(lfac2) ) lfac2(ixO^S) = lfac2_tmp(ixO^S)
          if ( present(lfac) ) lfac(ixO^S) = dsqrt( lfac2_tmp(ixO^S) )
          if ( present(v2) ) v2(ixO^S) = W2v2(ixO^S) / lfac2_tmp(ixO^S)
       end if
    end if

    if ( present(h) ) then
       ! Calculate the magnetically modified specific enthalpy
       h(ixO^S) = 1.0d0 + prim(ixO^S, eps_) & 
                + prim(ixO^S, press_) / prim(ixO^S, rho_) 
    end if
  end subroutine grhd_ccsn_leakage_mpi_get_intermediate_variables

  subroutine grhd_ccsn_leakage_mpi_get_lambda(ixI^L, ixO^L, idim, prim, x, lambda)
     use mod_global_parameters
     integer, intent(in)             :: ixI^L, ixO^L, idim
     double precision, intent(in)    :: prim(ixI^S, 1:nprim)
     double precision, intent(in)    :: x(ixI^S, 1:ndim)
     double precision, intent(out)   :: lambda(ixI^S, 1:2)

     double precision                :: gamma(ixI^S,1:3,1:3)
     double precision                :: h(ixI^S)
     double precision                :: cs2(ixI^S)
     double precision                :: v2(ixI^S)
     double precision                :: lfac(ixI^S)
     double precision                :: tmp1(ixI^S), tmp2(ixI^S)
     double precision                :: vel(ixI^S, 1:ndir)
     double precision                :: clight(ixI^S)
     double precision                :: inv_gamma_ii(ixI^S)
     integer                         :: idir


     call grhd_ccsn_leakage_mpi_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), v2=v2(ixI^S), lfac=lfac(ixI^S))

     do idir = 1, ndir
        vel(ixO^S,idir) = prim(ixO^S, W_vel(idir)) / lfac(ixO^S)
     end do

     ! beware of coordinate singularities
     where ( gamma(ixO^S, idim, idim) > smalldouble )
        ! sound speed square
        cs2(ixO^S) = prim(ixO^S, cs2_)
        inv_gamma_ii(ixO^S) = 1.0d0 / gamma(ixO^S, idim, idim)
     else where
        inv_gamma_ii(ixO^S) = 0.0d0
        vel(ixO^S, idim) = 0.0d0
        cs2(ixO^S) = 0.0d0
        v2(ixO^S) = 0.0d0
     end where

     tmp1(ixO^S) = vel(ixO^S, idim) * ( 1.0d0 - cs2(ixO^S) )
     tmp2(ixO^S) = dsqrt( cs2(ixO^S) * ( 1.0d0 - v2(ixO^S) ) * &
                    ( ( 1.0d0 - v2(ixO^S) * cs2(ixO^S) ) * inv_gamma_ii(ixO^S) &
                    - vel(ixO^S, idim)**2 * ( 1.0d0 - cs2(ixO^S) ) ) )

     lambda(ixO^S,1) = ( tmp1(ixO^S) + tmp2(ixO^S) ) / ( 1.0d0 - v2(ixO^S) * cs2(ixO^S) )
     lambda(ixO^S,2) = ( tmp1(ixO^S) - tmp2(ixO^S) ) / ( 1.0d0 - v2(ixO^S) * cs2(ixO^S) )

     ! limit with speed of light
     clight(ixO^S) = dsqrt( inv_gamma_ii(ixO^S) )
     lambda(ixO^S,1) = max( min( lambda(ixO^S,1), clight(ixO^S) ), -clight(ixO^S) )
     lambda(ixO^S,2) = max( min( lambda(ixO^S,2), clight(ixO^S) ), -clight(ixO^S) )

     !lambda(ixO^S,0) = prim(ixO^S, alp_) * prim(ixO^S, vel(idim)) - prim(ixO^S, beta(idim))  
     lambda(ixO^S,1) = prim(ixO^S, alp_) * lambda(ixO^S,1) - prim(ixO^S, beta(idim)) 
     lambda(ixO^S,2) = prim(ixO^S, alp_) * lambda(ixO^S,2) - prim(ixO^S, beta(idim)) 
  end subroutine grhd_ccsn_leakage_mpi_get_lambda

  !> This subroutine fix the abnormal values in conserved variables !
  subroutine grhd_ccsn_leakage_mpi_handle_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    use mod_geometry
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: prim(ixI^S,1:nprim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                      :: idir
    integer                      :: ix^D
    double precision             :: eps_min, eps_max
    logical                      :: need_to_fix_eos

 ! avoid coordinate singularities
    ! fixme: this might depends on different bc, but in general this should work
    if ( coordinate /= cartesian ) then
       where ( dabs(x(ixO^S,1)) < smalldouble ) 
          prim(ixO^S, W_vel(1)) = 0.0d0
          !prim(ixO^S, veloc(1)) = 0.0d0
          prim(ixO^S, beta(1)) = 0.0d0
       end where
       {^NOONED
       ! fixme: need to check
       if ( coordinate == spherical ) then
          where ( dabs(x(ixO^S,2)) < smalldouble ) 
             prim(ixO^S, W_vel(2)) = 0.0d0
             !prim(ixO^S, veloc(2)) = 0.0d0
             prim(ixO^S, beta(2)) = 0.0d0
          end where
       end if
       }
    end if



    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       {do ix^D = ixO^LIM^D \}
          need_to_fix_eos = .False.
          if ( prim(ix^D, rho_) < small_rho_thr ) then
             ! atmosphere handling
             call eos_get_eps_range( prim(ix^D, rho_), eps_min, eps_max)
!             prim(ix^D, rho_) = max(small_rho, eos_rhomin)
             prim(ix^D, rho_) = small_rho
!             prim(ix^D, eps_) = max(small_eps, eps_min)
             prim(ix^D, eps_) = small_eps
             prim(ix^D, W_vel(:)) = 0.0d0
             need_to_fix_eos = .True.
          else
             ! this is not atmosphere
             call eos_get_eps_range( prim(ix^D, rho_), eps_min, eps_max)
             if ( ( prim(ix^D, eps_) < eps_min ) .or. ( prim(ix^D, eps_) > eps_max ) ) then
                prim(ix^D, eps_) = max( min( eps_max, prim(ix^D, eps_) ), eps_min )
                need_to_fix_eos = .True.
             end if
          end if
          if ( need_to_fix_eos .and. update_eos ) then
             call grhd_ccsn_leakage_mpi_update_eos_one_point(prim(ix^D,1:nprim))
          end if
       {enddo^D&\}
       !where ( dabs(prim(ixO^S, W_vel(ndir))) <= smalldouble )
       !   prim(ixO^S, W_vel(ndir)) = 0.0d0
       !end where
    case default
       ! nothing to do here
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       return
    end select
  end subroutine grhd_ccsn_leakage_mpi_handle_small_values

end module mod_grhd_ccsn_leakage_mpi_phys_parameters
