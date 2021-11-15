!> General relativistic radiation hydrodynamics physics module (only in CFC now)
module mod_grrhd_phys
  use mod_physics
  use mod_grrhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grrhd_phys_init

contains

  !> Read this module's parameters from a file
  subroutine grrhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grrhd_list/ evolve_hydro, use_GR, tolerance, iter_max,&
        use_radiation

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grrhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine grrhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine grrhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!grrhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine grrhd_write_info


  !> Initialize the module
  subroutine grrhd_phys_init()
    use mod_global_parameters
    use mod_grrhd_phys_convert
    use mod_grrhd_phys_flux
    use mod_grrhd_phys_add_source

    integer :: itr, idir

    call grrhd_read_params(par_files)

    physics_type = "grrhd"
  
    ! Determine primitive variables that needed to be reconstructed
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)
    if ( use_radiation ) then
       nu_E = var_set_primvar('nu_E', var_type=hydro_var,need_bc=.false.)
       allocate(nu_F_over_E(ndir))
       do idir = 1, ndir
          nu_F_over_E(idir) = var_set_primvar('nu_F_over_E',idir,&
              var_type=hydro_var,need_bc=.false.)
       end do
    end if

    ! Determine primitive variables that needed NOT to be reconstructed
    press_ = var_set_primvar('press', var_type=hydro_var,need_bc=.false.,&
       need_rec=.false.)
    cs2_ = var_set_primvar('cs2', var_type=hydro_var,need_bc=.false.,&
       need_rec=.false.)
    if ( use_radiation ) then
       ! fixme: nu_zeta actually needed to be reconstructed, 
       ! but normally will not use the same scheme as in hydro.
       ! Most of the time PC is used(first order scheme), so I put it here
       nu_zeta = var_set_primvar('nu_zeta', var_type=hydro_var,&
          need_bc=.false.)
    end if


    ! Set index of metric variables, it has to be reconstructed so use fluxvar
    alp_ = var_set_primvar("alp", var_type=metric_var, need_bc=.False.,&
       need_rec=.false.)
    psi_ = var_set_primvar("psi", var_type=metric_var, need_bc=.False.,&
       need_rec=.false.)
    allocate(beta(ndir))
    do idir=1,ndir
       beta(idir) = var_set_primvar("beta", idir, var_type=metric_var,&
           need_bc=.False.,need_rec=.false.)
    enddo
    allocate(vecX(ndir))
    do idir=1,ndir
       vecX(idir) = var_set_primvar("X", idir, var_type=metric_var,&
           need_bc=.False.,need_rec=.false.)
    enddo

    ! Determine flux (cons) variables
    D_ = var_set_consvar('D', var_type=hydro_var)
    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = var_set_consvar('mom',idir, var_type=hydro_var)
    end do
    tau_ = var_set_consvar('tau', var_type=hydro_var)
    if ( use_radiation ) then
       nu_Econs = var_set_consvar('nu_Econs', var_type=hydro_var)
       allocate(nu_Fcons(ndir))
       do idir = 1, ndir
          nu_Fcons(idir) = var_set_consvar('nu_Fcons',idir,&
              var_type=hydro_var)
       end do
    end if

    phys_get_dt              => grrhd_get_dt
    phys_get_lfac2           => grrhd_get_lfac2
    phys_get_tilde_S         => grrhd_get_tilde_S
    phys_get_csound2         => grrhd_get_csound2
    phys_get_cmax            => grrhd_get_cmax
    phys_update_eos          => grrhd_update_eos
    phys_check_params        => grrhd_check_params
    phys_check_prim          => grrhd_check_prim
    phys_write_info          => grrhd_write_info
    phys_handle_small_values => grrhd_handle_small_values

    call grrhd_phys_convert_init()
    call grrhd_phys_flux_init()
    call grrhd_phys_add_source_init()
    if ( use_radiation ) then
       phys_modify_wLR            => grrhd_modify_wLR
    end if

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, ncons))
    else if (any(shape(flux_type) /= [ndir, ncons])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    flux_type = flux_default

    if ( .not.evolve_hydro ) then
       flux_type(:, D_) = flux_nul
       flux_type(:, tau_) = flux_nul
       flux_type(:, mom(1:ndir)) = flux_nul
    end if

    nvector      = 3 ! No. vector vars
    if ( use_radiation ) then
       nvector = nvector + 1
    end if
    allocate(iw_vector(nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1
    iw_vector(3) = vecX(1) - 1
    if ( use_radiation ) then
       iw_vector(4) = nu_F_over_E(1) - 1
    end if

  end subroutine grrhd_phys_init

  subroutine grrhd_check_params
    use mod_global_parameters
    use mod_eos
    use mod_usr_methods

    ! check whether eos module has been loaded
    call eos_check
    if ( use_radiation ) then
       if (.not.use_imex_scheme ) then
          call mpistop&
("Error: currently only IMEX schemes are supported when use radiation.")
       end if

       if(.not. associated(usr_get_opacities)) then
          call mpistop("Error: usr_get_opacities must be defined in GRRHD")
       end if
    end if

  end subroutine grrhd_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grrhd_update_eos(ixImin1,ixImax1, ixOmin1,ixOmax1, prim)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)                       :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(inout)           :: prim(ixImin1:ixImax1,&
        1:nprim)
    integer                                   :: ix1

    ! rho and eps are given, update the rest of the primitive variables
    do ix1 = ixOmin1,ixOmax1 
       call grrhd_update_eos_one_point(prim(ix1,:))
    enddo
  end subroutine grrhd_update_eos

  !> Calculate cmax_idim within ixO^L
  subroutine grrhd_get_cmax(prim, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
      cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1, idim
    double precision, intent(in)              :: prim(ixImin1:ixImax1,&
        1:nprim), x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1)
    integer                                   :: idir, ix1
    double precision                          :: lambda(ixImin1:ixImax1,1:2)
    call grrhd_get_lambda(ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lambda(ixImin1:ixImax1,1:2))
    cmax(ixOmin1:ixOmax1) = max(abs(lambda(ixOmin1:ixOmax1,1)),&
        abs(lambda(ixOmin1:ixOmax1,2)))

  end subroutine grrhd_get_cmax

  subroutine grrhd_get_csound2(ixImin1,ixImax1,ixOmin1,ixOmax1,prim,csound2)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: prim(ixImin1:ixImax1,nprim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1)
    integer                         :: idir, ix1
    ! Note: the input are the prim variables
    do ix1 = ixOmin1,ixOmax1 
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix1),prim(ix1, rho_),prim(ix1, eps_))
       ! strictly require cs2 is physical
       csound2(ix1) = max( min( 1.0d0, csound2(ix1) ) , 0.0d0)
    enddo
  end subroutine grrhd_get_csound2

  subroutine grrhd_get_lfac2(ps_in,ixImin1,ixImax1,ixOmin1,ixOmax1,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixImin1:ixImax1)
    integer                         :: idir, ix1
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:1), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma(ixImin1:ixImax1,1:3,1:3))
    do idir = 1, ndir
       gamma(ixOmin1:ixOmax1,idir,idir) = gamma(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, psi_)**4 
    end do
    
    lfac2(ixOmin1:ixOmax1) = 1.0d0 
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixOmin1:ixOmax1) = lfac2(ixOmin1:ixOmax1) + gamma(ixOmin1:ixOmax1,&
          idir,idir) * prim(ixOmin1:ixOmax1, W_vel(idir))**2
    end do
    end associate
  end subroutine grrhd_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine grrhd_get_tilde_S(ps_in,ixImin1,ixImax1,ixOmin1,ixOmax1,tilde_S)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixImin1:ixImax1)

    double precision                :: h(ixImin1:ixImax1)   
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call grrhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3), h=h(ixImin1:ixImax1) )

    ! nota that S = S_i v^i + 3 p
    tilde_S(ixOmin1:ixOmax1) = 3.0d0 * prim(ixOmin1:ixOmax1, press_)
    do idir = 1, ndir
       tilde_S(ixOmin1:ixOmax1) = tilde_S(ixOmin1:ixOmax1) + ( ( &
          prim(ixOmin1:ixOmax1, rho_) * h(ixOmin1:ixOmax1) ) * &
          prim(ixOmin1:ixOmax1, W_vel(idir)) ) * gamma(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, W_vel(idir))

    end do
    tilde_S(ixOmin1:ixOmax1) = tilde_S(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
        psi_)**6
    end associate
  end subroutine grrhd_get_tilde_S

  subroutine grrhd_get_dt(w, ixImin1,ixImax1, ixOmin1,ixOmax1, dtnew, dx1, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1, 1:1)
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nprim)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine grrhd_get_dt

  !> Returns 0 in argument flag where values are ok
  subroutine grrhd_check_prim(ixImin1,ixImax1, ixOmin1,ixOmax1, prim, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: prim(ixImin1:ixImax1, 1:nprim)
    integer, intent(inout)       :: flag(ixImin1:ixImax1)
    flag(ixOmin1:ixOmax1) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where(prim(ixOmin1:ixOmax1, rho_) < smalldouble ) 
       flag(ixOmin1:ixOmax1) = rho_
    else where(prim(ixOmin1:ixOmax1, eps_) < smalldouble ) 
       flag(ixOmin1:ixOmax1) = eps_
    end where
  end subroutine grrhd_check_prim

  subroutine grrhd_modify_wLR(ixImin1,ixImax1, ixOmin1,ixOmax1, consL, consR,&
      primL, primR, s, idir)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(inout) :: consL(ixImin1:ixImax1,1:ncons),&
        consR(ixImin1:ixImax1,1:ncons)
    double precision, intent(inout) :: primL(ixImin1:ixImax1,1:nprim),&
        primR(ixImin1:ixImax1,1:nprim)
    type(state)                     :: s

    call grrhd_get_nu_zeta(ixImin1,ixImax1, ixOmin1,ixOmax1, primL, s%x)
    call grrhd_get_nu_zeta(ixImin1,ixImax1, ixOmin1,ixOmax1, primR, s%x)
  end subroutine grrhd_modify_wLR

end module mod_grrhd_phys
