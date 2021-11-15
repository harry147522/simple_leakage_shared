module mod_grhd_ccsn_phys_parameters
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
  subroutine grhd_ccsn_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)

    ! rho and eps are given, update the rest of the primitive variables
!    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_), ye=prim( ye_))
!    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_), ye=prim( ye_))

    call eos_eps_get_all_one_grid(rho = prim(rho_), temp = prim(temp_),&
        ye = prim(ye_), eps = prim(eps_), ent = prim(ent_), cs2 = prim(cs2_) ,&
        prs = prim(press_))

       if ( prim(temp_) < eos_tempmin) then
           prim(temp_) = eos_tempmin*1.05d0
       endif

    ! strictly require cs2 is physical
    if ( prim(cs2_) >= 1.0d0 ) then
       prim(cs2_) = 0.0d0
    else
       prim(cs2_) = max( prim(cs2_) , 0.0d0)
    end if
  end subroutine grhd_ccsn_update_eos_one_point

  !> get some useful variables from primitive
  subroutine grhd_ccsn_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,&
     ixOmax1, prim, x, gamma, lfac2, lfac, v2, h )
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                     :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(in)            :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)            :: x(ixImin1:ixImax1, 1:ndim)

    double precision, intent(out), optional :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision, intent(out), optional :: v2(ixImin1:ixImax1)        ! 
    double precision, intent(out), optional :: lfac2(ixImin1:ixImax1) !Lorentz factor square
    double precision, intent(out), optional :: lfac(ixImin1:ixImax1) !Lorentz factor
    double precision, intent(out), optional :: h(ixImin1:ixImax1) !enthalpy: h

    integer                                 :: idir
    double precision                        :: W2v2(ixImin1:ixImax1)        ! 
    double precision                        :: v_bar(ixImin1:ixImax1,1:3) 
    double precision                        :: lfac2_tmp(ixImin1:ixImax1) !Lorentz factor square
    double precision                        :: lfac_tmp(ixImin1:ixImax1) !Lorentz factor
    double precision                        :: gamma_tmp(ixImin1:ixImax1,1:3,&
       1:3)

    if ( present(gamma) .or. present(lfac2) .or. present(lfac) .or. &
       present(v2) ) then
       ! get the metric
       call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1,&
           ixOmin1,ixOmax1, gamma_tmp(ixImin1:ixImax1,1:3,1:3))
       do idir = 1, ndir
          gamma_tmp(ixOmin1:ixOmax1,idir,idir) = gamma_tmp(ixOmin1:ixOmax1,&
             idir,idir) * prim(ixOmin1:ixOmax1, psi_)**4 
       end do
       if ( present(gamma) ) then
          gamma(ixOmin1:ixOmax1,1:3,1:3) = gamma_tmp(ixOmin1:ixOmax1,1:3,1:3) 
       end if
   
       if ( present(lfac2) .or. present(lfac) .or. present(v2) ) then
          W2v2 = 0.0d0
          ! calculate W^2 * v^2 first
          do idir = 1, ndir
             W2v2(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) + &
                gamma_tmp(ixOmin1:ixOmax1,idir,idir)*prim(ixOmin1:ixOmax1,&
                 W_vel(idir))**2
          end do
          ! Calculate the Lorentz factor from velocity
          lfac2_tmp(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) + 1.0d0 
          if ( present(lfac2) ) lfac2(ixOmin1:ixOmax1) = &
             lfac2_tmp(ixOmin1:ixOmax1)
          if ( present(lfac) ) lfac(ixOmin1:ixOmax1) = dsqrt( &
             lfac2_tmp(ixOmin1:ixOmax1) )
          if ( present(v2) ) v2(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) / &
             lfac2_tmp(ixOmin1:ixOmax1)
       end if
    end if

    if ( present(h) ) then
       ! Calculate the magnetically modified specific enthalpy
       h(ixOmin1:ixOmax1) = 1.0d0 + prim(ixOmin1:ixOmax1,&
           eps_) + prim(ixOmin1:ixOmax1, press_) / prim(ixOmin1:ixOmax1,&
           rho_) 
    end if
  end subroutine grhd_ccsn_get_intermediate_variables

  subroutine grhd_ccsn_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim, prim,&
      x, lambda)
     use mod_global_parameters
     integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
     double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
     double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
     double precision, intent(out)   :: lambda(ixImin1:ixImax1, 1:2)

     double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
     double precision                :: h(ixImin1:ixImax1)
     double precision                :: cs2(ixImin1:ixImax1)
     double precision                :: v2(ixImin1:ixImax1)
     double precision                :: lfac(ixImin1:ixImax1)
     double precision                :: tmp1(ixImin1:ixImax1),&
         tmp2(ixImin1:ixImax1)
     double precision                :: vel(ixImin1:ixImax1, 1:ndir)
     double precision                :: clight(ixImin1:ixImax1)
     double precision                :: inv_gamma_ii(ixImin1:ixImax1)
     integer                         :: idir


     call grhd_ccsn_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,&
        ixOmax1, prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
         gamma=gamma(ixImin1:ixImax1,1:3,1:3), v2=v2(ixImin1:ixImax1),&
         lfac=lfac(ixImin1:ixImax1))

     do idir = 1, ndir
        vel(ixOmin1:ixOmax1,idir) = prim(ixOmin1:ixOmax1,&
            W_vel(idir)) / lfac(ixOmin1:ixOmax1)
     end do

     ! beware of coordinate singularities
     where ( gamma(ixOmin1:ixOmax1, idim, idim) > smalldouble )
        ! sound speed square
        cs2(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1, cs2_)
        inv_gamma_ii(ixOmin1:ixOmax1) = 1.0d0 / gamma(ixOmin1:ixOmax1, idim,&
            idim)
     else where
        inv_gamma_ii(ixOmin1:ixOmax1) = 0.0d0
        vel(ixOmin1:ixOmax1, idim) = 0.0d0
        cs2(ixOmin1:ixOmax1) = 0.0d0
        v2(ixOmin1:ixOmax1) = 0.0d0
     end where

     tmp1(ixOmin1:ixOmax1) = vel(ixOmin1:ixOmax1,&
         idim) * ( 1.0d0 - cs2(ixOmin1:ixOmax1) )
     tmp2(ixOmin1:ixOmax1) = dsqrt( cs2(ixOmin1:ixOmax1) * ( 1.0d0 - &
        v2(ixOmin1:ixOmax1) ) * ( ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
        cs2(ixOmin1:ixOmax1) ) * inv_gamma_ii(ixOmin1:ixOmax1) - &
        vel(ixOmin1:ixOmax1, idim)**2 * ( 1.0d0 - cs2(ixOmin1:ixOmax1) ) ) )

     lambda(ixOmin1:ixOmax1,1) = ( tmp1(ixOmin1:ixOmax1) + &
        tmp2(ixOmin1:ixOmax1) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
        cs2(ixOmin1:ixOmax1) )
     lambda(ixOmin1:ixOmax1,2) = ( tmp1(ixOmin1:ixOmax1) - &
        tmp2(ixOmin1:ixOmax1) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
        cs2(ixOmin1:ixOmax1) )

     ! limit with speed of light
     clight(ixOmin1:ixOmax1) = dsqrt( inv_gamma_ii(ixOmin1:ixOmax1) )
     lambda(ixOmin1:ixOmax1,1) = max( min( lambda(ixOmin1:ixOmax1,1),&
         clight(ixOmin1:ixOmax1) ), -clight(ixOmin1:ixOmax1) )
     lambda(ixOmin1:ixOmax1,2) = max( min( lambda(ixOmin1:ixOmax1,2),&
         clight(ixOmin1:ixOmax1) ), -clight(ixOmin1:ixOmax1) )

     !lambda(ixO^S,0) = prim(ixO^S, alp_) * prim(ixO^S, vel(idim)) - prim(ixO^S, beta(idim))  
     lambda(ixOmin1:ixOmax1,1) = prim(ixOmin1:ixOmax1,&
         alp_) * lambda(ixOmin1:ixOmax1,1) - prim(ixOmin1:ixOmax1,&
         beta(idim)) 
     lambda(ixOmin1:ixOmax1,2) = prim(ixOmin1:ixOmax1,&
         alp_) * lambda(ixOmin1:ixOmax1,2) - prim(ixOmin1:ixOmax1,&
         beta(idim)) 
  end subroutine grhd_ccsn_get_lambda

  !> This subroutine fix the abnormal values in conserved variables !
  subroutine grhd_ccsn_handle_small_values(prim, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(inout) :: prim(ixImin1:ixImax1,1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                      :: idir
    integer                      :: ix1
    double precision             :: eps_min, eps_max
    logical                      :: need_to_fix_eos

 ! avoid coordinate singularities
    ! fixme: this might depends on different bc, but in general this should work
    if ( coordinate /= cartesian ) then
       where ( dabs(x(ixOmin1:ixOmax1,1)) < smalldouble ) 
          prim(ixOmin1:ixOmax1, W_vel(1)) = 0.0d0
          !prim(ixO^S, veloc(1)) = 0.0d0
          prim(ixOmin1:ixOmax1, beta(1)) = 0.0d0
       end where
       
    end if



    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       do ix1 = ixOmin1,ixOmax1 
          need_to_fix_eos = .False.
          if ( prim(ix1, rho_) < small_rho_thr ) then
             ! atmosphere handling
             call eos_get_eps_range( prim(ix1, rho_), eps_min, eps_max)
!             prim(ix^D, rho_) = max(small_rho, eos_rhomin)
             prim(ix1, rho_) = small_rho
!             prim(ix^D, eps_) = max(small_eps, eps_min)
             prim(ix1, eps_) = small_eps
             prim(ix1, W_vel(:)) = 0.0d0
             need_to_fix_eos = .True.
          else
             ! this is not atmosphere
             call eos_get_eps_range( prim(ix1, rho_), eps_min, eps_max)
             if ( ( prim(ix1, eps_) < eps_min ) .or. ( prim(ix1,&
                 eps_) > eps_max ) ) then
                prim(ix1, eps_) = max( min( eps_max, prim(ix1, eps_) ),&
                    eps_min )
                need_to_fix_eos = .True.
             end if
          end if
          if ( need_to_fix_eos .and. update_eos ) then
             call grhd_ccsn_update_eos_one_point(prim(ix1,1:nprim))
          end if
       enddo
       !where ( dabs(prim(ixO^S, W_vel(ndir))) <= smalldouble )
       !   prim(ixO^S, W_vel(ndir)) = 0.0d0
       !end where
    case default
       ! nothing to do here
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       return
    end select
  end subroutine grhd_ccsn_handle_small_values

end module mod_grhd_ccsn_phys_parameters
