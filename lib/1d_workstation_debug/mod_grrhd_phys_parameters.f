module mod_grrhd_phys_parameters
  use mod_physics
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: evolve_hydro = .True.
  logical                      :: use_process_adv_global = .false.
  logical                      :: use_GR = .false.
  logical                      :: use_radiation = .false.

  !-------------------------------------------------------------------!
  ! Parameters for con2prim
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tolerance = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_max = 50000
  !> maximum Lorentz factor
  double precision             :: lfac_max = 11.0d0
  double precision             :: v_max, k_max

  !-------------------------------------------------------------------!
  ! Parameters for radiation transport
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tol_rad = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_rad_max = 50000

  !-------------------------------------------------------------------!

  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grrhd_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)

    ! rho and eps are given, update the rest of the primitive variables
    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_))
    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_))
    ! strictly require cs2 is physical
    prim(cs2_) = max( min( 1.0d0, prim(cs2_) ) , 0.0d0)
  end subroutine grrhd_update_eos_one_point

  subroutine grrhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim, prim, x,&
      lambda)
     use mod_global_parameters
     integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
     double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
     double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
     double precision, intent(out)   :: lambda(ixImin1:ixImax1, 1:2)

     integer                         :: idir
     integer                         :: ix1
     double precision                :: tmp1(ixImin1:ixImax1),&
         tmp2(ixImin1:ixImax1)
     double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
     double precision                :: cs2(ixImin1:ixImax1),&
         v2(ixImin1:ixImax1)
     double precision                :: lfac(ixImin1:ixImax1),&
         lfac2(ixImin1:ixImax1)
     double precision                :: pi(ixImin1:ixImax1)
     double precision                :: chi(ixImin1:ixImax1) !closure function
     double precision                :: coeff_thin(ixImin1:ixImax1),&
         coeff_thick(ixImin1:ixImax1)
     double precision                :: lambda_thin(ixImin1:ixImax1, 1:2),&
         lambda_thick(ixImin1:ixImax1, 1:2)
     double precision                :: vel(ixImin1:ixImax1, 1:ndir)

     call grrhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
         prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
         gamma=gamma(ixImin1:ixImax1,1:3,1:3), v2=v2(ixImin1:ixImax1),&
         lfac2=lfac2(ixImin1:ixImax1), lfac=lfac(ixImin1:ixImax1))
     if (.not.use_radiation) then
        ! normal GRHD lambda
        ! sound speed square
        cs2(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1, cs2_)
   
        do idir = 1, ndir
           vel(ixOmin1:ixOmax1,idir) = prim(ixOmin1:ixOmax1,&
               W_vel(idir)) / lfac(ixOmin1:ixOmax1)
        end do
   
        tmp1(ixOmin1:ixOmax1) = vel(ixOmin1:ixOmax1,&
            idim) * ( 1.0d0 - cs2(ixOmin1:ixOmax1) )
        tmp2(ixOmin1:ixOmax1) = dsqrt( cs2(ixOmin1:ixOmax1) * ( 1.0d0 - &
           v2(ixOmin1:ixOmax1) ) * ( ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
           cs2(ixOmin1:ixOmax1) ) / gamma(ixOmin1:ixOmax1, idim,&
            idim) - vel(ixOmin1:ixOmax1, idim)**2 * ( 1.0d0 - &
           cs2(ixOmin1:ixOmax1) ) ) )
   
        lambda(ixOmin1:ixOmax1,1) = ( tmp1(ixOmin1:ixOmax1) + &
           tmp2(ixOmin1:ixOmax1) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
           cs2(ixOmin1:ixOmax1) )
        lambda(ixOmin1:ixOmax1,2) = ( tmp1(ixOmin1:ixOmax1) - &
           tmp2(ixOmin1:ixOmax1) ) / ( 1.0d0 - v2(ixOmin1:ixOmax1) * &
           cs2(ixOmin1:ixOmax1) )
   
        ! limit with speed of light, here we use tmp1 to store it
        tmp1(ixOmin1:ixOmax1) = dsqrt( 1.0d0 / gamma(ixOmin1:ixOmax1, idim,&
            idim) )
        lambda(ixOmin1:ixOmax1,1) = max( min( lambda(ixOmin1:ixOmax1,1),&
            tmp1(ixOmin1:ixOmax1) ), -tmp1(ixOmin1:ixOmax1) )
        lambda(ixOmin1:ixOmax1,2) = max( min( lambda(ixOmin1:ixOmax1,2),&
            tmp1(ixOmin1:ixOmax1) ), -tmp1(ixOmin1:ixOmax1) )
   
        !lambda(ixO^S,0) = prim(ixO^S, alp_) * prim(ixO^S, vel(idim)) - prim(ixO^S, beta(idim))  
        lambda(ixOmin1:ixOmax1,1) = prim(ixOmin1:ixOmax1,&
            alp_) * lambda(ixOmin1:ixOmax1,1) - prim(ixOmin1:ixOmax1,&
            beta(idim)) 
        lambda(ixOmin1:ixOmax1,2) = prim(ixOmin1:ixOmax1,&
            alp_) * lambda(ixOmin1:ixOmax1,2) - prim(ixOmin1:ixOmax1,&
            beta(idim)) 
     else
        ! get F -> F * E
        do idir = 1, ndir
           vel(ixOmin1:ixOmax1, idir) = prim(ixOmin1:ixOmax1,&
               nu_F_over_E(idir)) * prim(ixOmin1:ixOmax1, nu_E) 
        end do
        ! get sqrt F
        tmp2 = 0.0d0
        do idir = 1, ndir
           tmp2(ixOmin1:ixOmax1) = tmp2(ixOmin1:ixOmax1) + &
              gamma(ixOmin1:ixOmax1,idir,idir) * vel(ixOmin1:ixOmax1, idir)**2
        end do
        tmp2(ixOmin1:ixOmax1) = dsqrt( tmp2(ixOmin1:ixOmax1) )
    
        ! fixme: need to up index
        tmp1(ixOmin1:ixOmax1) = dabs( vel(ixOmin1:ixOmax1,&
            idim) ) / tmp2(ixOmin1:ixOmax1)
    
        lambda_thin(ixOmin1:ixOmax1,1) = - prim(ixOmin1:ixOmax1,&
            beta(idim)) + prim(ixOmin1:ixOmax1, alp_) * tmp1(ixOmin1:ixOmax1) 
        lambda_thin(ixOmin1:ixOmax1,2) = - prim(ixOmin1:ixOmax1,&
            beta(idim)) - prim(ixOmin1:ixOmax1, alp_) * tmp1(ixOmin1:ixOmax1) 

        ! now work out lambda_thick
        pi(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
            alp_) * prim(ixOmin1:ixOmax1,&
            W_vel(idim)) / lfac2(ixOmin1:ixOmax1)
        tmp1(ixOmin1:ixOmax1) = 2.0d0 * prim(ixOmin1:ixOmax1,&
            alp_) * prim(ixOmin1:ixOmax1, W_vel(idim))
        tmp2(ixOmin1:ixOmax1) = dsqrt( prim(ixOmin1:ixOmax1,&
            alp_) * ( 2.0d0 * lfac2(ixOmin1:ixOmax1) + 1.0d0 ) / &
           gamma(ixOmin1:ixOmax1, idim, idim) - 2.0d0 * lfac2(ixOmin1:ixOmax1) &
           * pi(ixOmin1:ixOmax1)**2 ) 

        lambda_thick(ixOmin1:ixOmax1,1) = - prim(ixOmin1:ixOmax1,&
            beta(idim)) + ( tmp1(ixOmin1:ixOmax1) + tmp2(ixOmin1:ixOmax1) ) / &
           ( 2.0d0 * lfac2(ixOmin1:ixOmax1) + 1.0d0 )
        lambda_thick(ixOmin1:ixOmax1,2) = - prim(ixOmin1:ixOmax1,&
            beta(idim)) + ( tmp1(ixOmin1:ixOmax1) - tmp2(ixOmin1:ixOmax1) ) / &
           ( 2.0d0 * lfac2(ixOmin1:ixOmax1) + 1.0d0 )  

        lambda_thick(ixOmin1:ixOmax1,1) = min( - prim(ixOmin1:ixOmax1,&
            beta(idim)) + pi(ixOmin1:ixOmax1), lambda_thick(ixOmin1:ixOmax1,&
           1) )
        lambda_thick(ixOmin1:ixOmax1,2) = min( - prim(ixOmin1:ixOmax1,&
            beta(idim)) + pi(ixOmin1:ixOmax1), lambda_thick(ixOmin1:ixOmax1,&
           2) )
    
        ! now combine two lambdas
        ! loop over points
        do ix1 = ixOmin1,ixOmax1 
           chi(ix1) = minerbo_closure_func( prim(ix1, nu_zeta) )
           coeff_thin(ix1) = c_thin(chi(ix1))
           coeff_thick(ix1) = 1.0d0 - coeff_thin(ix1)
        end do ! end loop for every points
        lambda(ixOmin1:ixOmax1,1) = coeff_thin(ixOmin1:ixOmax1) * &
           lambda_thin(ixOmin1:ixOmax1,1) + coeff_thick(ixOmin1:ixOmax1) * &
           lambda_thick(ixOmin1:ixOmax1,1)
        lambda(ixOmin1:ixOmax1,2) = coeff_thin(ixOmin1:ixOmax1) * &
           lambda_thin(ixOmin1:ixOmax1,2) + coeff_thick(ixOmin1:ixOmax1) * &
           lambda_thick(ixOmin1:ixOmax1,2)
     end if
  end subroutine grrhd_get_lambda

  !> This subroutine fix the abnormal values in conserved variables !
  subroutine grrhd_handle_small_values(prim, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(inout) :: prim(ixImin1:ixImax1,1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                      :: idir
    integer                      :: ix1
    double precision             :: eps_min, eps_max
    logical                      :: need_to_fix_eos

    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       do ix1 = ixOmin1,ixOmax1 
          need_to_fix_eos = .False.
          if ( prim(ix1, rho_) < small_rho_thr ) then
             ! atmosphere handling
             prim(ix1, rho_) = small_rho
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
             call grrhd_update_eos_one_point(prim(ix1,1:nprim))
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
  end subroutine grrhd_handle_small_values

  !> get some useful variables from primitive
  subroutine grrhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
      prim, x, Pij, gamma, lfac2, lfac, v2, h )
    use mod_global_parameters
    use mod_geometry
    use mod_rootfinding
    integer, intent(in)                     :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(in)            :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)            :: x(ixImin1:ixImax1, 1:ndim)

    double precision, intent(out), optional :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision, intent(out), optional :: v2(ixImin1:ixImax1)        ! 
    double precision, intent(out), optional :: lfac2(ixImin1:ixImax1) !Lorentz factor square
    double precision, intent(out), optional :: lfac(ixImin1:ixImax1) !Lorentz factor
    double precision, intent(out), optional :: h(ixImin1:ixImax1) !enthalpy: h
    double precision, intent(out), optional :: Pij(ixImin1:ixImax1,1:3,1:3)

    integer                                 :: idir, jdir
    double precision                        :: W2v2(ixImin1:ixImax1),&
        v2_tmp(ixImin1:ixImax1) !
    !double precision                        :: v_bar(ixI^S,1:3) 
    double precision                        :: lfac2_tmp(ixImin1:ixImax1) !Lorentz factor square
    double precision                        :: lfac_tmp(ixImin1:ixImax1) !Lorentz factor
    double precision                        :: gamma_tmp(ixImin1:ixImax1,1:3,&
       1:3)
    double precision                        :: E(ixImin1:ixImax1),&
        F2(ixImin1:ixImax1), Fi(ixImin1:ixImax1,1:ndir) !F_i
    double precision                        :: v_dot_F(ixImin1:ixImax1) !vi F_i
    double precision                        :: gamma_H(ixImin1:ixImax1,1:ndir)
    integer                                 :: ix1 
    ! --------------   variables needed to get Pij  ------------- 
    double precision                        :: chi(ixImin1:ixImax1) !closure function
    double precision                        :: coeff_thin(ixImin1:ixImax1),&
        coeff_thick(ixImin1:ixImax1)
    double precision                        :: J_over_3(ixImin1:ixImax1)

    if ( present(h) ) then
       ! Calculate the magnetically modified specific enthalpy
       h(ixOmin1:ixOmax1) = 1.0d0 + prim(ixOmin1:ixOmax1,&
           eps_) + prim(ixOmin1:ixOmax1, press_) / prim(ixOmin1:ixOmax1,&
           rho_) 
    end if

    if ( .not. ( present(gamma) .or. present(lfac2) .or. present(lfac) .or. &
       present(v2) .or. present(Pij) ) ) then
       return
    end if

    ! get the metric
    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma_tmp(ixImin1:ixImax1,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixOmin1:ixOmax1,idir,idir) = gamma_tmp(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, psi_)**4 
    end do
    if ( present(gamma) ) then
       gamma(ixOmin1:ixOmax1,1:3,1:3) = gamma_tmp(ixOmin1:ixOmax1,1:3,1:3) 
    end if

    W2v2 = 0.0d0
    ! calculate W^2 * v^2 first
    do idir = 1, ndir
       W2v2(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) + &
          gamma_tmp(ixOmin1:ixOmax1,idir,idir)*prim(ixOmin1:ixOmax1,&
           W_vel(idir))**2
    end do
    ! Calculate the Lorentz factor from velocity
    lfac2_tmp(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) + 1.0d0 
    if ( present(lfac2) ) lfac2(ixOmin1:ixOmax1) = lfac2_tmp(ixOmin1:ixOmax1)
    lfac_tmp(ixOmin1:ixOmax1) = dsqrt( lfac2_tmp(ixOmin1:ixOmax1) )
    if ( present(lfac) ) lfac(ixOmin1:ixOmax1) = lfac_tmp(ixOmin1:ixOmax1)
    v2_tmp(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) / &
       lfac2_tmp(ixOmin1:ixOmax1)
    if ( present(v2) ) v2(ixOmin1:ixOmax1) = v2_tmp(ixOmin1:ixOmax1)

    if ( .not. present(Pij)  ) then
       return
    end if

    ! some useful relations
    E(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1, nu_E)
    ! Note that the momentum density in prim is F_i / E
    do idir = 1, ndir
       Fi(ixOmin1:ixOmax1, idir) = prim(ixOmin1:ixOmax1,&
           nu_F_over_E(idir)) * E(ixOmin1:ixOmax1)
    end do
    ! calculate the F^i F_i
    F2 = 0.0d0
    do idir = 1, ndir
       F2(ixOmin1:ixOmax1) = F2(ixOmin1:ixOmax1) + Fi(ixOmin1:ixOmax1,&
           idir)**2 / gamma(ixOmin1:ixOmax1, idir, idir)
    end do
    ! v^i * F_i
    v_dot_F = 0.0d0
    do idir = 1, ndir
       v_dot_F(ixOmin1:ixOmax1) = v_dot_F(ixOmin1:ixOmax1) + &
          prim(ixOmin1:ixOmax1, W_vel(idir)) / lfac_tmp(ixOmin1:ixOmax1) * &
          Fi(ixOmin1:ixOmax1, idir)
    end do

    if ( present(Pij) ) then
       ! loop over points
       do ix1 = ixOmin1,ixOmax1 
          chi(ix1) = minerbo_closure_func( prim(ix1, nu_zeta) )
          coeff_thin(ix1) = c_thin(chi(ix1))
          coeff_thick(ix1) = 1.0d0 - coeff_thin(ix1)
       end do ! end loop for every points

       ! J and H can be computed here if needed

       ! optically thin part of pressure tensor
       do idir = 1, ndir; do jdir = 1, ndir
          Pij(ixOmin1:ixOmax1,idir,jdir) = coeff_thin(ixOmin1:ixOmax1) * &
             E(ixOmin1:ixOmax1)/F2(ixOmin1:ixOmax1) * Fi(ixOmin1:ixOmax1,&
             idir) * Fi(ixOmin1:ixOmax1,jdir)
       end do; end do

       ! optically thick part of pressure tensor
       J_over_3(ixOmin1:ixOmax1) = 1.0d0 / (2.0d0 * lfac2_tmp(ixOmin1:ixOmax1) &
          + 1.0d0) * ( (2.0d0 * lfac2_tmp(ixOmin1:ixOmax1) - 1.0d0) * &
          E(ixOmin1:ixOmax1) - 2.0d0 * lfac2_tmp(ixOmin1:ixOmax1) * &
          v_dot_F(ixOmin1:ixOmax1) )
       do idir = 1, ndir
          gamma_H(ixOmin1:ixOmax1,idir) = Fi(ixOmin1:ixOmax1,&
             idir)/gamma(ixOmin1:ixOmax1,idir,&
             idir)/lfac_tmp(ixOmin1:ixOmax1) + prim(ixOmin1:ixOmax1,&
             W_vel(idir)) / (2.0d0 * lfac2_tmp(ixOmin1:ixOmax1) + 1.0d0) * ( &
             (4.0d0 * lfac2_tmp(ixOmin1:ixOmax1) + 1.0d0) * &
             v_dot_F(ixOmin1:ixOmax1) - 4.0d0 * lfac2_tmp(ixOmin1:ixOmax1) * &
             E(ixOmin1:ixOmax1) )
       end do

       do idir = 1, ndir
          Pij(ixOmin1:ixOmax1,idir,idir) = Pij(ixOmin1:ixOmax1,idir,&
             idir) + coeff_thick(ixOmin1:ixOmax1) * ( &
             J_over_3(ixOmin1:ixOmax1) * ( 4.0d0 * prim(ixOmin1:ixOmax1,&
             W_vel(idir)) * prim(ixOmin1:ixOmax1,&
             W_vel(idir)) + 1.0d0 / gamma(ixOmin1:ixOmax1,idir,idir)) )
       end do

       do idir = 1, ndir; do jdir = 1, ndir
          Pij(ixOmin1:ixOmax1,idir,jdir) = Pij(ixOmin1:ixOmax1,idir,&
             jdir) + coeff_thick(ixOmin1:ixOmax1) * ( &
             lfac_tmp(ixOmin1:ixOmax1) * (gamma_H(ixOmin1:ixOmax1,&
             idir) * prim(ixOmin1:ixOmax1,&
              W_vel(jdir)) / lfac_tmp(ixOmin1:ixOmax1) &
             +gamma_H(ixOmin1:ixOmax1,jdir) * prim(ixOmin1:ixOmax1,&
              W_vel(idir)) / lfac_tmp(ixOmin1:ixOmax1) ) ) 
       end do; end do
    end if
  end subroutine grrhd_get_intermediate_variables

  !> the implementation here follows SpECTRE https://spectre-code.org/
  subroutine grrhd_get_nu_zeta(ixImin1,ixImax1, ixOmin1,ixOmax1, prim, x,&
      v2_in, lfac_in, lfac2_in, gamma_in )
    use mod_global_parameters
    use mod_geometry
    use mod_rootfinding
    implicit none
    integer, intent(in)                     :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(inout)         :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)            :: x(ixImin1:ixImax1, 1:ndim)

    double precision, intent(in), optional  :: gamma_in(ixImin1:ixImax1,1:3,&
       1:3)
    double precision, intent(in), optional  :: v2_in(ixImin1:ixImax1) !
    double precision, intent(in), optional  :: lfac2_in(ixImin1:ixImax1) !Lorentz factor square
    double precision, intent(in), optional  :: lfac_in(ixImin1:ixImax1) !Lorentz factor

    integer                                 :: ix1 
    integer                                 :: idir, jdir
    double precision                        :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision                        :: v2(ixImin1:ixImax1),&
        lfac2(ixImin1:ixImax1), lfac(ixImin1:ixImax1) 
    integer                                 :: error_code = -1
    ! --------------   variables needed to solve zeta  ------------- 
    double precision                        :: E, F2, Fi(1:ndir) ! F_i
    double precision                        :: v_dot_F  !  vi F_i
    double precision                        :: zeta !variable Eddington factor
    double precision                        :: zeta_m, zeta_p
    double precision                        :: chi ! closure function    
    double precision                        :: coeff_thin, coeff_thick

    ! fluid-frame energy density:
    double precision                        :: J_0, J_thin, J_thick
    ! fluid-frame momentum density:
    double precision                        :: H_0_t, H_thin_t, H_thick_t
    double precision                        :: H_0_v, H_thin_v, H_thick_v
    double precision                        :: H_0_F, H_thin_F, H_thick_F
    ! fluid-frame momentum density square and the relevent variables:
    double precision                        :: H2_0
    double precision                        :: H2_thin, H2_thick
    double precision                        :: H2_thin_thin, H2_thick_thick
    double precision                        :: H2_thin_thick

    if ( .not. ( present(gamma_in) .or. present(lfac_in) .or. &
       present(lfac2_in) .or. present(v2_in) ) ) then
     call grrhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
         prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
         gamma=gamma(ixImin1:ixImax1,1:3,1:3), v2=v2(ixImin1:ixImax1),&
         lfac2=lfac2(ixImin1:ixImax1), lfac=lfac(ixImin1:ixImax1))
    else
       gamma(ixImin1:ixImax1,1:3,1:3) = gamma_in(ixImin1:ixImax1,1:3,1:3) 
       lfac(ixImin1:ixImax1) = lfac_in(ixImin1:ixImax1) 
       lfac2(ixImin1:ixImax1) = lfac2_in(ixImin1:ixImax1) 
       v2(ixImin1:ixImax1) = v2_in(ixImin1:ixImax1) 
    end if

    ! loop over points
    do ix1 = ixOmin1,ixOmax1 
       E = prim(ix1, nu_E)
       ! Note that the momentum density in prim is F_i / E
       Fi(1:ndir) = prim(ix1, nu_F_over_E(1:ndir)) * E
   
       ! calculate the F^i F_i
       F2 = 0.0d0
       do idir = 1, ndir
          F2 = F2 + Fi(idir)**2 / gamma(ix1,idir,idir)
       end do
       ! if the velocity can be ignored, assume that
       ! fluid frame = inertial frame
       if ( v2(ix1) <= smalldouble ) then
          zeta = min( max( dsqrt( F2 )/ E, 0.0d0), 1.0d0)
      
       ! if the velocity can NOT be ignored, we need to go through a more
       ! expensive closure calculation
       else
          ! initial guess of the closure factor
          zeta = min( max( prim(ix1, nu_zeta), 0.0d0), 1.0d0)
      
          ! v^i * F_i
          v_dot_F = 0.0d0
          do idir = 1, ndir
             v_dot_F = v_dot_F + prim(ix1, W_vel(idir)) / lfac(ix1) * Fi(idir)
          end do
          ! Decomposition of the fluid-frame energy density:
          ! J = J_0 + coeff_thin * J_Thin + coeff_thick * J_Thick
          J_0 = lfac2(ix1) * ( E - 2.0d0 * v_dot_F )
          J_thin = lfac2(ix1) * E * v_dot_F**2 / F2 
          J_thick = ( lfac2(ix1) - 1.0d0 ) / ( 1.0d0 + 2.0d0 * lfac2(ix1) ) * &
             ( 4.0d0 * lfac2(ix1) * v_dot_F + E * ( 3.0d0 - 2.0d0 * lfac2(ix1) &
             ) )
      
          ! Decomposition of the fluid-frame momentum density:
          ! H_a = -( H_0_t + coeff_thick H_thick_t + coeff_thin H_thin_t) t_a
          !  - ( H_0_v + coeff_thick H_thick_v + coeff_thin H_thin_v) v_a
          !  - ( H_0_F + coeff_thick H_thick_F + coeff_thin H_thin_F) F_a
          ! with t_a the unit normal, v_a the 3-velocity, and F_a the
          ! inertial frame momentum density. 
          !
          ! This is a decomposition of convenience, 
          ! which is not unique: F_a and v_a are not
          ! orthogonal vectors, but both are normal to t_a.
          H_0_t = lfac(ix1) * ( J_0 + v_dot_F - E )
          H_0_v = lfac(ix1) * J_0 
          H_0_F = - lfac(ix1)
          H_thin_t = lfac(ix1) * J_thin
          H_thin_v = H_thin_t
          H_thin_F = lfac(ix1) * E * v_dot_F / F2
          H_thick_t = lfac(ix1) * J_thick
          H_thick_v = H_thick_t + lfac(ix1) / ( 2.0d0 * lfac2(ix1) + 1.0d0 ) * &
             ( (3.0d0 - 2.0d0 * lfac2(ix1)) * E - (1.0d0 - 2.0d0 * lfac2(ix1)) &
             * v_dot_F )
          H_thick_F = lfac2(ix1) * v2(ix1)
          ! Quantities needed for the computation of H^2 = H^a H_a,
          ! independent of zeta. We write:
          ! H^2 = H2_0 + H2_thin * coeff_thin + H2_thick*coeff_thick
          ! + H2_thin_thin * coeff_thin^2 + H2_thick_thick * coeff_thick^2
          ! + H2_thin_thick * coeff_thin * coeff_thick;
          H2_0 = - H_0_t**2 + H_0_v**2 * v2(ix1) + H_0_F**2 * F2 + 2.0d0 * &
             H_0_v * H_0_F * v_dot_F
          H2_thin = 2.0d0 * (H_0_v * H_thin_v * v2(ix1) + H_0_F * H_thin_F * &
             F2 + H_0_v * H_thin_F * v_dot_F + H_0_F * H_thin_v * v_dot_F - &
             H_0_t * H_thin_t)
          H2_thick = 2.0d0 * (H_0_v * H_thick_v * v2(ix1) + H_0_F * H_thick_F &
             * F2 + H_0_v * H_thick_F * v_dot_F + H_0_F * H_thick_v * v_dot_F &
             - H_0_t * H_thick_t)
          H2_thin_thick = 2.0d0 * (H_thin_v * H_thick_v * v2(ix1) + H_thin_F * &
             H_thick_F * F2 + H_thin_v * H_thick_F * v_dot_F + H_thin_F * &
             H_thick_v * v_dot_F - H_thin_t * H_thick_t)
          H2_thick_thick = (H_thick_v**2) * v2(ix1) + (H_thick_F**2) * F2 + &
             2.0d0 * H_thick_v * H_thick_F * v_dot_F - (H_thick_t**2) 
          H2_thin_thin = (H_thin_v**2) * v2(ix1) + (H_thin_F**2) * F2 + 2.0d0 &
             * H_thin_v * H_thin_F * v_dot_F - (H_thin_t**2)
      
          ! find zeta though root finding
      
          ! atempt Newton-Raphson method to find the root
          ! To avoid failures in Newton-Raphson solver at the boundary of
          ! the allowed domain for zeta, test the edge values first.
          zeta_m = 0.0d0
          zeta_p = 1.0d0
          if ( dabs(zeta_J2_minus_H2(zeta_m)) <= tol_rad ) then
             zeta = zeta_m
          else if ( dabs(zeta_J2_minus_H2(zeta_p)) <= tol_rad ) then
             zeta = zeta_p
          else
             call rootfinding_constrained_newton_raphson(zeta, zeta_m, zeta_p,&
                 tol_rad, 50, error_code, zeta_J2_minus_H2,&
                 dev_zeta_J2_minus_H2)
             ! just in case if newton raphson method fails, use illinois method
             if( error_code/=0 ) then
                write(*,*)&
"Warning: Newton-Raphson method failed to find zeta, now try illinois method"
                call rootfinding_illinois(zeta, zeta_m, zeta_p, tol_rad,&
                    iter_rad_max, error_code, zeta_J2_minus_H2)
             end if
             ! check the root
             select case (error_code)
             case (-1) ! nothing happened
                call mpistop&
                   ("have you ever attemp to find the root in con2prim?")
             case (1) ! failed to find the root
                call mpistop("cannot find the root in zeta_J2_minus_H2")
             case (2) ! z is NaN
                call mpistop("NaN in zeta_J2_minus_H2")
             case (3) ! z is not bracketed
                call mpistop("the root is not bracketed in con2prim")
             end select
          end if
       end if
       prim(ix1, nu_zeta) = zeta
    end do ! end loop for every points

    contains
       double precision function zeta_J2_minus_H2(zeta) result(f)
          double precision, intent(in)     :: zeta
          double precision                 :: chi
          double precision                 :: coeff_thin
          double precision                 :: coeff_thick
          double precision                 :: J
          double precision                 :: H2
   
          chi = minerbo_closure_func(zeta)
          coeff_thin = c_thin(chi)
          coeff_thick = 1.0d0 - coeff_thin
   
          J = J_0 + coeff_thin * J_Thin + coeff_thick * J_Thick
   
          H2 = H2_0 + H2_thick * coeff_thick + H2_thin * coeff_thin + &
             H2_thin_thin * (coeff_thin**2) + H2_thick_thick * &
             (coeff_thick**2) + H2_thin_thick * coeff_thin * coeff_thick
   
          f = ( ( zeta * J )**2 - H2 ) / E**2
       end function zeta_J2_minus_H2
   
       double precision function dev_zeta_J2_minus_H2(zeta) result(f)
          double precision, intent(in)     :: zeta
          double precision                 :: chi, dchi
          double precision                 :: coeff_thin, coeff_thick
          double precision                 :: dcoeff_thin, dcoeff_thick
          double precision                 :: J, H2
          double precision                 :: dJ, dH2, dH2_dcoeff_thin,&
              dH2_dcoeff_thick
   
          chi = minerbo_closure_func(zeta)
          dchi = minerbo_closure_deriv(zeta)
          coeff_thin = c_thin(chi)
          coeff_thick = 1.0d0 - coeff_thin
          dcoeff_thin = dc_thin_dchi(chi) * dchi
          dcoeff_thick = dc_thick_dchi(chi) * dchi
   
          J = J_0 + coeff_thin * J_Thin + coeff_thick * J_Thick
   
          H2 = H2_0 + H2_thick * coeff_thick + H2_thin * coeff_thin + &
             H2_thin_thin * (coeff_thin**2) + H2_thick_thick * &
             (coeff_thick**2) + H2_thin_thick * coeff_thin * coeff_thick
   
          dJ = dcoeff_thin * J_Thin + dcoeff_thick * J_Thick
   
          dH2_dcoeff_thin = H2_thin + H2_thin_thin * (2.0d0 * coeff_thin) + &
             H2_thin_thick * coeff_thick
   
          dH2_dcoeff_thick = H2_thick + H2_thick_thick * (2.0d0 * coeff_thick) &
             + H2_thin_thick * coeff_thin
   
          dH2 = dH2_dcoeff_thin * dcoeff_thin + dH2_dcoeff_thick * &
             dcoeff_thick
   
          f = ( 2.0d0 * ( zeta * J**2 + zeta**2 * J * dJ) - dH2 ) / E**2
       end function dev_zeta_J2_minus_H2

 end subroutine grrhd_get_nu_zeta

  ! --------------------------
  !  some useful functions
  ! --------------------------
  double precision function minerbo_closure_func(zeta) result(chi)
    double precision, intent(in)     :: zeta
    chi = 1.0d0/3.0d0 + zeta**2 * (0.4d0 - 2.0d0/15.0d0 * zeta + 0.4d0 * &
       zeta**2)
    !if ( (chi > 1.0d0) .or. (chi < 1.0d0/3.0d0) ) stop "wrong chi"
  end function minerbo_closure_func

  double precision function minerbo_closure_deriv(zeta) result(dchi)
    double precision, intent(in)     :: zeta
    dchi = 0.4d0 * zeta * ( 2.0d0 - zeta + 4.0d0 * zeta**2)
  end function minerbo_closure_deriv

  double precision function c_thin(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = 1.5d0 * chi - 0.5d0
  end function c_thin

  double precision function c_thick(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = 1.5d0 * ( 1.0d0 - chi )
  end function c_thick

  double precision function dc_thin_dchi(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = 1.5d0 
  end function dc_thin_dchi

  double precision function dc_thick_dchi(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = - 1.5d0
  end function dc_thick_dchi

end module mod_grrhd_phys_parameters
