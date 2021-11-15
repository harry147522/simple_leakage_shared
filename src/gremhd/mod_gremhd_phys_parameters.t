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
  subroutine gremhd_get_intermediate_variables(ixI^L, ixO^L, prim, x, h,&
             gamma, sqrt_gamma, lfac2, lfac, v2, v_hat, &
             E2, B2, emu, bmu, qrho_e)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: prim(ixI^S, 1:nprim)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(out), optional :: sqrt_gamma(ixI^S)
    double precision, intent(out), optional :: gamma(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: v_hat(ixI^S,1:ndir)   
    double precision, intent(out), optional :: v2(ixI^S)        ! 
    double precision, intent(out), optional :: lfac2(ixI^S) ! Lorentz factor square
    double precision, intent(out), optional :: lfac(ixI^S) ! Lorentz factor
    double precision, intent(out), optional :: h(ixI^S) ! enthalpy
    double precision, intent(out), optional :: E2(ixI^S), B2(ixI^S) 
    double precision, intent(out), optional :: emu(ixI^S,0:ndir)   ! projection of E^mu along fluid four-velocity u^nu
    double precision, intent(out), optional :: bmu(ixI^S,0:ndir)   ! projection of B^mu along fluid four-velocity u^nu
    double precision, intent(out), optional :: qrho_e(ixI^S) ! conformal rescaled charge density

    integer                                 :: idir, jdir, kdir
    double precision                        :: W2v2(ixI^S)        ! 
    double precision                        :: lfac2_tmp(ixI^S) ! Lorentz factor square
    double precision                        :: lfac_tmp(ixI^S) ! Lorentz factor
    double precision                        :: gamma_tmp(ixI^S,1:3,1:3)
    double precision                        :: sqrt_gamma_tmp(ixI^S)
    double precision                        :: qE_field(ixI^S, 1:ndir)   ! conformally rescaled E-field

    double precision                        :: nmu(ixI^S,0:ndir)
    double precision                        :: B_mu(ixI^S,0:ndir), E_mu(ixI^S,0:ndir)
    double precision                        :: tmp_var(ixI^S)

    if ( present(h) ) then
       ! Calculate the specific enthalpy
       h(ixO^S) = 1.0d0 + prim(ixO^S, eps_) & 
                + prim(ixO^S, press_) / prim(ixO^S, rho_) 
    end if

    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    if ( present(gamma) ) then
       gamma(ixO^S,1:3,1:3) = gamma_tmp(ixO^S,1:3,1:3) 
    end if

    call get_sqrt_gamma_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, sqrt_gamma_tmp(ixI^S))
    sqrt_gamma_tmp(ixO^S) = sqrt_gamma_tmp(ixO^S) * prim(ixO^S, psi_)**6
    if ( present(sqrt_gamma) ) then
       sqrt_gamma(ixO^S) = sqrt_gamma_tmp(ixO^S)
    end if

    W2v2 = 0.0d0
    ! calculate W^2 * v^2 first
    do idir = 1, ndir
       W2v2(ixO^S) = W2v2(ixO^S) + gamma_tmp(ixO^S,idir,idir)*prim(ixO^S, W_vel(idir))**2
    end do

    ! Calculate the Lorentz factor from velocity
    lfac2_tmp(ixO^S) = W2v2(ixO^S) + 1.0d0 
    lfac_tmp(ixO^S) = dsqrt( lfac2_tmp(ixO^S) )
    if ( present(lfac2) ) lfac2(ixO^S) = lfac2_tmp(ixO^S)
    if ( present(lfac) ) lfac(ixO^S) = lfac_tmp(ixO^S)
    if ( present(v2) ) v2(ixO^S) = W2v2(ixO^S) / lfac2_tmp(ixO^S)

    if ( present(v_hat) ) then
       ! v_hat^i = v^i * alp - beta
       do idir = 1, ndir
          v_hat(ixO^S, idir) = prim(ixO^S, alp_) * prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) - prim(ixO^S, beta(idir))
       end do
    end if

    if ( present(E2) ) then
       E2 = 0.0d0
       do idir = 1, ndir
          E2(ixO^S) = E2(ixO^S) + prim(ixO^S, Evec(idir))**2 * gamma_tmp(ixO^S, idir, idir)
       end do
    end if
    if ( present(B2) ) then
       B2 = 0.0d0
       do idir = 1, ndir
          B2(ixO^S) = B2(ixO^S) + prim(ixO^S, Bvec(idir))**2 * gamma_tmp(ixO^S, idir, idir)
       end do
    end if

    if ( present(emu) .or. present(bmu) ) then
       nmu(ixO^S,0) = 1.0d0 / prim(ixO^S, alp_)
       do idir = 1, ndir
          nmu(ixO^S, idir) = - prim(ixO^S, beta(idir)) / prim(ixO^S, alp_)
       end do

       if (present(emu)) then
          emu = 0.0d0
          do idir = 1, ndir
             do jdir = 1, ndir
                do kdir = 1, ndir
                   emu(ixO^S, idir) = emu(ixO^S, idir) &
                          + lvc(idir,jdir,kdir) * prim(ixO^S, W_vel(jdir)) * prim(ixO^S, Bvec(kdir))
                end do
             end do
             emu(ixO^S, idir) = emu(ixO^S, idir) / sqrt_gamma_tmp(ixO^S) &
                                    + prim(ixO^S, Evec(idir)) * lfac_tmp(ixO^S)
          end do
          ! then get E_dot_Wv term
          tmp_var = 0.0d0
          do idir = 1, ndir
             tmp_var(ixO^S) = tmp_var(ixO^S) + prim(ixO^S, Evec(idir)) * prim(ixO^S, W_vel(idir)) 
          end do
          do idir = 0, ndir
             emu(ixO^S, idir) = emu(ixO^S, idir) + tmp_var(ixO^S) * nmu(ixO^S, idir)
          end do
       end if

       if (present(bmu)) then
          bmu = 0.0d0
          do idir = 1, ndir
             do jdir = 1, ndir
                do kdir = 1, ndir
                   bmu(ixO^S, idir) = bmu(ixO^S, idir) &
                          + lvc(idir,jdir,kdir) * prim(ixO^S, W_vel(jdir)) * prim(ixO^S, Evec(kdir))
                end do
             end do
             bmu(ixO^S, idir) = - bmu(ixO^S, idir) / sqrt_gamma_tmp(ixO^S) &
                                    + prim(ixO^S, Bvec(idir)) * lfac_tmp(ixO^S)
          end do
          ! then get B_dot_Wv term
          tmp_var = 0.0d0
          do idir = 1, ndir
             tmp_var(ixO^S) = tmp_var(ixO^S) + prim(ixO^S, Bvec(idir)) * prim(ixO^S, W_vel(idir)) 
          end do
          do idir = 0, ndir
             bmu(ixO^S, idir) = bmu(ixO^S, idir) + tmp_var(ixO^S) * nmu(ixO^S, idir)
          end do
       end if
    end if

    if (present(qrho_e)) then
       do idir = 1, ndir
          qE_field(ixI^S, idir) = prim(ixI^S, Evec(idir)) * prim(ixI^S, psi_)**6
       end do
       call gremhd_get_div(ixI^L, ixO^L, qE_field(ixI^S, 1:ndir), qrho_e(ixI^S), dive_4thorder)
    end if
  end subroutine gremhd_get_intermediate_variables

  !> Calculate div vector within ixO
  subroutine gremhd_get_div(ixI^L,ixO^L, vector, div_vec, fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: vector(ixI^S,1:ndir)
    double precision, intent(inout) :: div_vec(ixI^S)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: divb_corner(ixI^S), sign
    double precision                   :: aux_vol(ixI^S)
    integer                            :: ixC^L, idir, ic^D, ix^L

    if(stagger_grid) then
      stop "fixme"
    else
      !bvector(ixI^S,:)=cons(ixI^S,Bcons(:))
      !select case(typediv)
      !case("central")
        call divvector(vector(ixI^S,1:ndir),ixI^L,ixO^L,div_vec(ixI^S),fourthorder)
      !case("limited")
        !call divvectorS(bvec,ixI^L,ixO^L,divb)
      !end select
    end if
  end subroutine gremhd_get_div

  !> This subroutine fix the abnormal values in primitive variables !
  subroutine gremhd_handle_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: prim(ixI^S,1:nprim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                      :: idir
    integer                      :: ix^D
    double precision             :: eps_min, eps_max
    logical                      :: need_to_fix_eos

    ! fixme: delete small_values_fix_iw

    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       {do ix^D = ixO^LIM^D \}
          need_to_fix_eos = .False.
          if ( prim(ix^D, rho_) < small_rho_thr ) then
             ! atmosphere handling
             prim(ix^D, rho_) = small_rho
             prim(ix^D, eps_) = small_eps
             prim(ix^D, W_vel) = 0.0d0
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
             call gremhd_update_eos_one_point(prim(ix^D,1:nprim))
          end if
       {enddo^D&\}
    case default
       ! nothing to do here
       return
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
    end select
  end subroutine gremhd_handle_small_values

end module mod_gremhd_phys_parameters
