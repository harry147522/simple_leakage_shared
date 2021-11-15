module mod_grhd_phys_parameters
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
  double precision             :: lfac_max = 1.1d1
  double precision             :: v_max, k_max

  !-------------------------------------------------------------------!
  ! Parameters for 1D core treatment
  !-------------------------------------------------------------------!
  logical                      :: oneDcore = .False.
  double precision             :: r_core = 0.0d0

  !-------------------------------------------------------------------!

  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grhd_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)

    ! rho and eps are given, update the rest of the primitive variables
    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_))
    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_))
    ! strictly require cs2 is physical
    if ( prim(cs2_) > 1.0d0 ) then
       prim(cs2_) = 0.0d0
    else
       prim(cs2_) = max( prim(cs2_) , 0.0d0)
    end if
  end subroutine grhd_update_eos_one_point

  !> get some useful variables from primitive
  subroutine grhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim, x, gamma, lfac2, lfac, v2, h )
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                     :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)            :: prim(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nprim)
    double precision, intent(in)            :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:ndim)

    double precision, intent(out), optional :: gamma(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3,1:3)
    double precision, intent(out), optional :: v2(ixImin1:ixImax1,&
       ixImin2:ixImax2) !
    double precision, intent(out), optional :: lfac2(ixImin1:ixImax1,&
       ixImin2:ixImax2) !Lorentz factor square
    double precision, intent(out), optional :: lfac(ixImin1:ixImax1,&
       ixImin2:ixImax2) !Lorentz factor
    double precision, intent(out), optional :: h(ixImin1:ixImax1,&
       ixImin2:ixImax2) !enthalpy: h

    integer                                 :: idir
    double precision                        :: W2v2(ixImin1:ixImax1,&
       ixImin2:ixImax2) !
    double precision                        :: v_bar(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3) 
    double precision                        :: lfac2_tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2) !Lorentz factor square
    double precision                        :: lfac_tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2) !Lorentz factor
    double precision                        :: gamma_tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3,1:3)

    if ( present(gamma) .or. present(lfac2) .or. present(lfac) .or. &
       present(v2) ) then
       ! get the metric
       call get_gamma_ij_hat(x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
           ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
           gamma_tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3))
       do idir = 1, ndir
          gamma_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idir) = gamma_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, psi_)**4 
       end do
       if ( present(gamma) ) then
          gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3,&
             1:3) = gamma_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3,1:3) 
       end if
   
       if ( present(lfac2) .or. present(lfac) .or. present(v2) ) then
          W2v2 = 0.0d0
          ! calculate W^2 * v^2 first
          do idir = 1, ndir
             W2v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = W2v2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) + gamma_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                idir,idir)*prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 W_vel(idir))**2
          end do
          ! Calculate the Lorentz factor from velocity
          lfac2_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = W2v2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) + 1.0d0 
          if ( present(lfac2) ) lfac2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = lfac2_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          if ( present(lfac) ) lfac(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = dsqrt( lfac2_tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) )
          if ( present(v2) ) v2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = W2v2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) / lfac2_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    end if

    if ( present(h) ) then
       ! Calculate the magnetically modified specific enthalpy
       h(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.0d0 + prim(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, eps_) + prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           press_) / prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) 
    end if
  end subroutine grhd_get_intermediate_variables

  subroutine grhd_get_lambda(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, idim, prim, x, lambda)
     use mod_global_parameters
     integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
     double precision, intent(in)    :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nprim)
     double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndim)
     double precision, intent(out)   :: lambda(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:2)

     double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:3,1:3)
     double precision                :: h(ixImin1:ixImax1,ixImin2:ixImax2)
     double precision                :: cs2(ixImin1:ixImax1,ixImin2:ixImax2)
     double precision                :: v2(ixImin1:ixImax1,ixImin2:ixImax2)
     double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2)
     double precision                :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2),&
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
     double precision                :: vel(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndir)
     double precision                :: clight(ixImin1:ixImax1,&
        ixImin2:ixImax2)
     double precision                :: inv_gamma_ii(ixImin1:ixImax1,&
        ixImin2:ixImax2)
     integer                         :: idir


     call grhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
         gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
         v2=v2(ixImin1:ixImax1,ixImin2:ixImax2), lfac=lfac(ixImin1:ixImax1,&
        ixImin2:ixImax2))

     do idir = 1, ndir
        vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = prim(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, W_vel(idir)) / lfac(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)
     end do

     ! beware of coordinate singularities
     where ( gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idim,&
         idim) > smalldouble )
        ! sound speed square
        cs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prim(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, cs2_)
        inv_gamma_ii(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.0d0 / &
           gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idim, idim)
     else where
        inv_gamma_ii(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
        vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idim) = 0.0d0
        cs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
        v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
     end where

     tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = vel(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, idim) * ( 1.0d0 - cs2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) )
     tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dsqrt( cs2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) * ( 1.0d0 - v2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) ) * ( ( 1.0d0 - v2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) * cs2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) ) * inv_gamma_ii(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) - vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idim)**2 * ( 1.0d0 - cs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) ) ) )

     lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) = ( tmp1(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) + tmp2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) * cs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )
     lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) = ( tmp1(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) - tmp2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) * cs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )

     ! limit with speed of light
     clight(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dsqrt( &
        inv_gamma_ii(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )
     lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1) = max( min( lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1),&
         clight(ixOmin1:ixOmax1,ixOmin2:ixOmax2) ), -clight(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) )
     lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2) = max( min( lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2),&
         clight(ixOmin1:ixOmax1,ixOmin2:ixOmax2) ), -clight(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) )

     !lambda(ixO^S,0) = prim(ixO^S, alp_) * prim(ixO^S, vel(idim)) - prim(ixO^S, beta(idim))  
     lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) = prim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, alp_) * lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1) - prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, beta(idim)) 
     lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) = prim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, alp_) * lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2) - prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, beta(idim)) 
  end subroutine grhd_get_lambda

  !> This subroutine fix the abnormal values in conserved variables !
  subroutine grhd_handle_small_values(prim, x, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                      :: idir
    integer                      :: ix1,ix2
    double precision             :: eps_min, eps_max
    logical                      :: need_to_fix_eos

    ! avoid coordinate singularities
    ! fixme: this might depends on different bc, but in general this should work
    if ( coordinate /= cartesian ) then
       where ( dabs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)) < smalldouble ) 
          prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, W_vel(1)) = 0.0d0
          !prim(ixO^S, veloc(1)) = 0.0d0
          prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, beta(1)) = 0.0d0
       end where
    end if

    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       do ix1 = ixOmin1,ixOmax1 
       do ix2 = ixOmin2,ixOmax2 
          need_to_fix_eos = .False.
          if ( prim(ix1,ix2, rho_) < small_rho_thr ) then
             ! atmosphere handling
             prim(ix1,ix2, rho_) = small_rho
             prim(ix1,ix2, eps_) = small_eps
             prim(ix1,ix2, W_vel(:)) = 0.0d0
             need_to_fix_eos = .True.
          else
             ! this is not atmosphere
             call eos_get_eps_range( prim(ix1,ix2, rho_), eps_min, eps_max)
             if ( ( prim(ix1,ix2, eps_) < eps_min ) .or. ( prim(ix1,ix2,&
                 eps_) > eps_max ) ) then
                prim(ix1,ix2, eps_) = max( min( eps_max, prim(ix1,ix2, eps_) ),&
                    eps_min )
                need_to_fix_eos = .True.
             end if
          end if
          if ( need_to_fix_eos .and. update_eos ) then
             call grhd_update_eos_one_point(prim(ix1,ix2,1:nprim))
          end if
       enddo
       enddo
       !where ( dabs(prim(ixO^S, W_vel(ndir))) <= smalldouble )
       !   prim(ixO^S, W_vel(ndir)) = 0.0d0
       !end where
    case default
       ! nothing to do here
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       return
    end select
  end subroutine grhd_handle_small_values

end module mod_grhd_phys_parameters
