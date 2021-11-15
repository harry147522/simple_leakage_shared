module mod_gremhd_phys_parameters
  use mod_physics
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: evolve_hydro = .True.
  logical                      :: evolve_EM = .True.
  logical                      :: use_GR = .false.

  !-------------------------------------------------------------------!
  ! Parameters for con2prim
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tolerance = 1.0d-14
  !> maximum iteration for the root finding
  integer                      :: iter_max = 50000
  !> maximum Lorentz factor
  double precision             :: lfac_max = 1.1d2
  double precision             :: v_max, k_max

  !-------------------------------------------------------------------!
  ! Parameters for solving E fields
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tol_im = 1.0d-14
  !> maximum iteration for the root finding
  integer                      :: iter_max_im = 500

  !-------------------------------------------------------------------!
  ! Parameters for divergence handle
  !-------------------------------------------------------------------!
  logical                      :: dive_4thorder = .False.
  logical                      :: divb_4thorder = .False.
  integer                      :: type_divb = 0
  integer, parameter           :: none = 0
  integer, parameter           :: divb_multigrid = 1
  integer, parameter           :: divb_GLM = 2

  !-------------------------------------------------------------------!
  ! Parameters for divergence with Multigrid
  !-------------------------------------------------------------------!
  ! N smooths will be applied for 1: cycling up, 2: cycling down
  integer, dimension(1:2)                   :: divB_mg_n_cycle = (/5,5/)
  logical                                   :: divB_mg_redblack = .True.
  double precision                          :: divB_mg_tol = 1.0d-6
  integer                                   :: divB_mg_it_max = 50

  !-------------------------------------------------------------------!
  ! Parameters for divergence with GLM
  !-------------------------------------------------------------------!
  double precision                          :: divB_glm_kappa = 1.0d0


  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine gremhd_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)
    ! rho and eps are given, update the rest of the primitive variables
    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_))
    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_))
    ! strictly require cs2 is physical
    if ( prim(cs2_) >= 1.0d0 ) then
       prim(cs2_) = 0.0d0
    else
       prim(cs2_) = max( prim(cs2_) , 0.0d0)
    end if
  end subroutine gremhd_update_eos_one_point

  !> get some useful variables from primitive
  subroutine gremhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,&
     ixOmax1, prim, x, h,gamma, sqrt_gamma, lfac2, lfac, v2, v_hat, E2, B2,&
      emu, bmu, qrho_e)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                     :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(in)            :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)            :: x(ixImin1:ixImax1, 1:ndim)

    double precision, intent(out), optional :: sqrt_gamma(ixImin1:ixImax1)
    double precision, intent(out), optional :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision, intent(out), optional :: v_hat(ixImin1:ixImax1,&
       1:ndir)   
    double precision, intent(out), optional :: v2(ixImin1:ixImax1)        ! 
    double precision, intent(out), optional :: lfac2(ixImin1:ixImax1) !Lorentz factor square
    double precision, intent(out), optional :: lfac(ixImin1:ixImax1) !Lorentz factor
    double precision, intent(out), optional :: h(ixImin1:ixImax1) ! enthalpy
    double precision, intent(out), optional :: E2(ixImin1:ixImax1),&
        B2(ixImin1:ixImax1) 
    double precision, intent(out), optional :: emu(ixImin1:ixImax1,0:ndir) !projection of Emu along fluid four-velocity unu
    double precision, intent(out), optional :: bmu(ixImin1:ixImax1,0:ndir) !projection of Bmu along fluid four-velocity unu
    double precision, intent(out), optional :: qrho_e(ixImin1:ixImax1) !conformal rescaled charge density

    integer                                 :: idir, jdir, kdir
    double precision                        :: W2v2(ixImin1:ixImax1)        ! 
    double precision                        :: lfac2_tmp(ixImin1:ixImax1) !Lorentz factor square
    double precision                        :: lfac_tmp(ixImin1:ixImax1) !Lorentz factor
    double precision                        :: gamma_tmp(ixImin1:ixImax1,1:3,&
       1:3)
    double precision                        :: sqrt_gamma_tmp(ixImin1:ixImax1)
    double precision                        :: qE_field(ixImin1:ixImax1,&
        1:ndir) !conformally rescaled E-field

    double precision                        :: nmu(ixImin1:ixImax1,0:ndir)
    double precision                        :: B_mu(ixImin1:ixImax1,0:ndir),&
        E_mu(ixImin1:ixImax1,0:ndir)
    double precision                        :: tmp_var(ixImin1:ixImax1)

    if ( present(h) ) then
       ! Calculate the specific enthalpy
       h(ixOmin1:ixOmax1) = 1.0d0 + prim(ixOmin1:ixOmax1,&
           eps_) + prim(ixOmin1:ixOmax1, press_) / prim(ixOmin1:ixOmax1,&
           rho_) 
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

    call get_sqrt_gamma_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1,&
        ixOmin1,ixOmax1, sqrt_gamma_tmp(ixImin1:ixImax1))
    sqrt_gamma_tmp(ixOmin1:ixOmax1) = sqrt_gamma_tmp(ixOmin1:ixOmax1) * &
       prim(ixOmin1:ixOmax1, psi_)**6
    if ( present(sqrt_gamma) ) then
       sqrt_gamma(ixOmin1:ixOmax1) = sqrt_gamma_tmp(ixOmin1:ixOmax1)
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
    lfac_tmp(ixOmin1:ixOmax1) = dsqrt( lfac2_tmp(ixOmin1:ixOmax1) )
    if ( present(lfac2) ) lfac2(ixOmin1:ixOmax1) = lfac2_tmp(ixOmin1:ixOmax1)
    if ( present(lfac) ) lfac(ixOmin1:ixOmax1) = lfac_tmp(ixOmin1:ixOmax1)
    if ( present(v2) ) v2(ixOmin1:ixOmax1) = W2v2(ixOmin1:ixOmax1) / &
       lfac2_tmp(ixOmin1:ixOmax1)

    if ( present(v_hat) ) then
       ! v_hat^i = v^i * alp - beta
       do idir = 1, ndir
          v_hat(ixOmin1:ixOmax1, idir) = prim(ixOmin1:ixOmax1,&
              alp_) * prim(ixOmin1:ixOmax1,&
              W_vel(idir)) / lfac_tmp(ixOmin1:ixOmax1) - prim(ixOmin1:ixOmax1,&
              beta(idir))
       end do
    end if

    if ( present(E2) ) then
       E2 = 0.0d0
       do idir = 1, ndir
          E2(ixOmin1:ixOmax1) = E2(ixOmin1:ixOmax1) + prim(ixOmin1:ixOmax1,&
              Evec(idir))**2 * gamma_tmp(ixOmin1:ixOmax1, idir, idir)
       end do
    end if
    if ( present(B2) ) then
       B2 = 0.0d0
       do idir = 1, ndir
          B2(ixOmin1:ixOmax1) = B2(ixOmin1:ixOmax1) + prim(ixOmin1:ixOmax1,&
              Bvec(idir))**2 * gamma_tmp(ixOmin1:ixOmax1, idir, idir)
       end do
    end if

    if ( present(emu) .or. present(bmu) ) then
       nmu(ixOmin1:ixOmax1,0) = 1.0d0 / prim(ixOmin1:ixOmax1, alp_)
       do idir = 1, ndir
          nmu(ixOmin1:ixOmax1, idir) = - prim(ixOmin1:ixOmax1,&
              beta(idir)) / prim(ixOmin1:ixOmax1, alp_)
       end do

       if (present(emu)) then
          emu = 0.0d0
          do idir = 1, ndir
             do jdir = 1, ndir
                do kdir = 1, ndir
                   emu(ixOmin1:ixOmax1, idir) = emu(ixOmin1:ixOmax1,&
                       idir) + lvc(idir,jdir,kdir) * prim(ixOmin1:ixOmax1,&
                       W_vel(jdir)) * prim(ixOmin1:ixOmax1, Bvec(kdir))
                end do
             end do
             emu(ixOmin1:ixOmax1, idir) = emu(ixOmin1:ixOmax1,&
                 idir) / sqrt_gamma_tmp(ixOmin1:ixOmax1) + &
                prim(ixOmin1:ixOmax1, Evec(idir)) * lfac_tmp(ixOmin1:ixOmax1)
          end do
          ! then get E_dot_Wv term
          tmp_var = 0.0d0
          do idir = 1, ndir
             tmp_var(ixOmin1:ixOmax1) = tmp_var(ixOmin1:ixOmax1) + &
                prim(ixOmin1:ixOmax1, Evec(idir)) * prim(ixOmin1:ixOmax1,&
                 W_vel(idir)) 
          end do
          do idir = 0, ndir
             emu(ixOmin1:ixOmax1, idir) = emu(ixOmin1:ixOmax1,&
                 idir) + tmp_var(ixOmin1:ixOmax1) * nmu(ixOmin1:ixOmax1, idir)
          end do
       end if

       if (present(bmu)) then
          bmu = 0.0d0
          do idir = 1, ndir
             do jdir = 1, ndir
                do kdir = 1, ndir
                   bmu(ixOmin1:ixOmax1, idir) = bmu(ixOmin1:ixOmax1,&
                       idir) + lvc(idir,jdir,kdir) * prim(ixOmin1:ixOmax1,&
                       W_vel(jdir)) * prim(ixOmin1:ixOmax1, Evec(kdir))
                end do
             end do
             bmu(ixOmin1:ixOmax1, idir) = - bmu(ixOmin1:ixOmax1,&
                 idir) / sqrt_gamma_tmp(ixOmin1:ixOmax1) + &
                prim(ixOmin1:ixOmax1, Bvec(idir)) * lfac_tmp(ixOmin1:ixOmax1)
          end do
          ! then get B_dot_Wv term
          tmp_var = 0.0d0
          do idir = 1, ndir
             tmp_var(ixOmin1:ixOmax1) = tmp_var(ixOmin1:ixOmax1) + &
                prim(ixOmin1:ixOmax1, Bvec(idir)) * prim(ixOmin1:ixOmax1,&
                 W_vel(idir)) 
          end do
          do idir = 0, ndir
             bmu(ixOmin1:ixOmax1, idir) = bmu(ixOmin1:ixOmax1,&
                 idir) + tmp_var(ixOmin1:ixOmax1) * nmu(ixOmin1:ixOmax1, idir)
          end do
       end if
    end if

    if (present(qrho_e)) then
       do idir = 1, ndir
          qE_field(ixImin1:ixImax1, idir) = prim(ixImin1:ixImax1,&
              Evec(idir)) * prim(ixImin1:ixImax1, psi_)**6
       end do
       call gremhd_get_div(ixImin1,ixImax1, ixOmin1,ixOmax1,&
           qE_field(ixImin1:ixImax1, 1:ndir), qrho_e(ixImin1:ixImax1),&
           dive_4thorder)
    end if
  end subroutine gremhd_get_intermediate_variables

  !> Calculate div vector within ixO
  subroutine gremhd_get_div(ixImin1,ixImax1,ixOmin1,ixOmax1, vector, div_vec,&
      fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: vector(ixImin1:ixImax1,1:ndir)
    double precision, intent(inout) :: div_vec(ixImin1:ixImax1)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: divb_corner(ixImin1:ixImax1), sign
    double precision                   :: aux_vol(ixImin1:ixImax1)
    integer                            :: ixCmin1,ixCmax1, idir, ic1, ixmin1,&
       ixmax1

    if(stagger_grid) then
      stop "fixme"
    else
      !bvector(ixI^S,:)=cons(ixI^S,Bcons(:))
      !select case(typediv)
      !case("central")
        call divvector(vector(ixImin1:ixImax1,1:ndir),ixImin1,ixImax1,ixOmin1,&
           ixOmax1,div_vec(ixImin1:ixImax1),fourthorder)
      !case("limited")
        !call divvectorS(bvec,ixI^L,ixO^L,divb)
      !end select
    end if
  end subroutine gremhd_get_div

  !> This subroutine fix the abnormal values in primitive variables !
  subroutine gremhd_handle_small_values(prim, x, ixImin1,ixImax1, ixOmin1,&
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

    ! fixme: delete small_values_fix_iw

    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       do ix1 = ixOmin1,ixOmax1 
          need_to_fix_eos = .False.
          if ( prim(ix1, rho_) < small_rho_thr ) then
             ! atmosphere handling
             prim(ix1, rho_) = small_rho
             prim(ix1, eps_) = small_eps
             prim(ix1, W_vel) = 0.0d0
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
             call gremhd_update_eos_one_point(prim(ix1,1:nprim))
          end if
       enddo
    case default
       ! nothing to do here
       return
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
    end select
  end subroutine gremhd_handle_small_values

end module mod_gremhd_phys_parameters
