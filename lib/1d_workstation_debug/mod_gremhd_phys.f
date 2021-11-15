!> General relativistic hydrodynamics physics module (only in CFC now)
module mod_gremhd_phys
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private


  ! Public methods
  public :: gremhd_phys_init
  public :: gremhd_get_div

contains

  !> Read this module's parameters from a file
  subroutine gremhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    character(len=name_len)      :: type_divb_fix = ""

    namelist /gremhd_list/ evolve_hydro, evolve_EM, use_GR, tolerance,&
        iter_max, tol_im, iter_max_im, dive_4thorder, type_divb_fix,&
        divb_4thorder, divB_glm_kappa, divB_mg_n_cycle, divB_mg_redblack,&
        divB_mg_tol, divB_mg_it_max

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, gremhd_list, end=111)
111    close(unitpar)
    end do

    if (mype == 0) then
       if ( .not. evolve_hydro ) write(*,*) "Hydro is not evolving"
       if ( .not. evolve_EM ) write(*,*) "EM is not evolving"
    end if

    select case (type_divb_fix)
    case ('none')
       type_divb = none
    case ('multigrid')
       type_divb = divb_multigrid
    case ('glm')
       type_divb = divb_GLM
    case default
       if (mype == 0) write(*,*)&
           "type_divb_fix is not specificied, now we assume type_divb = none!"
       type_divb = none
       !call mpistop("this type_divb_fix is not supported")
    end select

  end subroutine gremhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine gremhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!gremhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine gremhd_write_info


  !> Initialize the module
  subroutine gremhd_phys_init()
    use mod_global_parameters
    use mod_gremhd_phys_convert
    use mod_gremhd_phys_flux
    use mod_gremhd_phys_add_source
    use mod_gremhd_phys_implicit_update
    use mod_gremhd_phys_divb_mg

    integer :: itr, idir

    call gremhd_read_params(par_files)

    physics_type = "gremhd"
  
    ! Determine primitive variables that needed to be reconstructed
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)

    allocate(Evec(ndir))
    do idir = 1, ndir
       Evec(idir) = var_set_primvar('Evec',idir, var_type=hydro_var,&
          need_bc=.false.)
    end do
    allocate(Bvec(ndir))
    do idir = 1, ndir
       Bvec(idir) = var_set_primvar('Bvec',idir, var_type=hydro_var,&
          need_bc=.false.)
    end do
    if (type_divb == divb_GLM) then
       ! introduce a scalar field Phi for divergence cleaning
       Bphi_ = var_set_primvar('Bphi', var_type=hydro_var,need_bc=.false.)
    end if

    ! Determine primitive variables that needed not to be reconstructed
    press_ = var_set_primvar('press', var_type=hydro_var,need_bc=.false.,&
       need_rec=.false.)
    cs2_ = var_set_primvar('cs2', var_type=hydro_var,need_bc=.false.,&
       need_rec=.false.)

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

    ! Determine hydro flux variables
    D_ = var_set_consvar('D', var_type=hydro_var)
    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = var_set_consvar('mom',idir, var_type=hydro_var)
    end do
    tau_ = var_set_consvar('tau', var_type=hydro_var)
    allocate(Econs(ndir))
    do idir = 1, ndir
       Econs(idir) = var_set_consvar('Econs',idir, var_type=hydro_var)
    end do
    allocate(Bcons(ndir))
    do idir = 1, ndir
       Bcons(idir) = var_set_consvar('Bcons',idir, var_type=hydro_var)
    end do
    if (type_divb == divb_GLM) then
       ! introduce a scalar field Phi for divergence cleaning
       Bphi_cons_ = var_set_consvar('Bphi_cons', var_type=hydro_var)
    end if

    phys_get_dt              => gremhd_get_dt
    phys_get_lfac2           => gremhd_get_lfac2
    phys_get_tilde_S         => gremhd_get_tilde_S
    phys_get_csound2         => gremhd_get_csound2
    phys_get_cmax            => gremhd_get_cmax
    phys_update_eos          => gremhd_update_eos
    phys_check_params        => gremhd_check_params
    phys_check_prim          => gremhd_check_prim
    phys_write_info          => gremhd_write_info
    phys_handle_small_values => gremhd_handle_small_values

    call gremhd_phys_convert_init()
    call gremhd_phys_flux_init()
    call gremhd_phys_add_source_init()
    ! fixme: if not evolve_EM, divB might have some problem
    if ( evolve_EM ) call gremhd_phys_implicit_update_init()
    
    select case (type_divb)
    case (divb_glm)
       ! nothing
    case (divb_multigrid)
       call gremhd_phys_divb_mg_init()
    case default
       ! nothing
    end select

    ! setup div cleaning here
    phys_step_adv_global => gremhd_step_adv_global

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
    if ( .not.evolve_EM ) then
       flux_type(:, Econs(1:ndir)) = flux_nul
       flux_type(:, Bcons(1:ndir)) = flux_nul
    end if

    nvector      = 5 ! No. vector vars
    allocate(iw_vector(1:nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1 
    iw_vector(3) = vecX(1) - 1
    iw_vector(4) = Evec(1) - 1
    iw_vector(5) = Bvec(1) - 1

  end subroutine gremhd_phys_init

  subroutine gremhd_check_params
    use mod_global_parameters
    use mod_eos
    use mod_usr_methods

    ! check whether eos module has been loaded
    call eos_check

    if (evolve_EM) then
       if (.not.use_imex_scheme ) then
          call mpistop&
             ("Error: currently only IMEX schemes are supported in gremhd.")
       end if
       
       if(.not. associated(usr_get_resistivity)) then
          call mpistop("Error: usr_get_resistivity must be defined in GRRMHD")
       end if
    end if

  end subroutine gremhd_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine gremhd_update_eos(ixImin1,ixImax1, ixOmin1,ixOmax1, prim)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1
    double precision, intent(inout)           :: prim(ixImin1:ixImax1,&
        1:nprim)
    integer                                   :: ix1
    ! rho and eps are given, update the rest of the primitive variables
    do ix1 = ixOmin1,ixOmax1 
       call gremhd_update_eos_one_point(prim(ix1, 1:nprim))
    enddo
  end subroutine gremhd_update_eos

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine gremhd_get_cmax(prim, x, ixImin1,ixImax1, ixOmin1,ixOmax1, idim,&
      cmax)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat
    integer, intent(in)                       :: ixImin1,ixImax1, ixOmin1,&
       ixOmax1, idim
    double precision, intent(in)              :: prim(ixImin1:ixImax1, nprim),&
        x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1)
    integer                                   :: idir, ix1
    double precision                          :: gamma_hat(ixImin1:ixImax1,1:3,&
       1:3)
    double precision, dimension(ixImin1:ixImax1,1:2)    :: lambda
    ! get the metric
    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma_hat(ixImin1:ixImax1,1:3,1:3))
    lambda(ixOmin1:ixOmax1,1) = dsqrt( 1.0d0 / gamma_hat(ixOmin1:ixOmax1, idim,&
        idim)) / prim(ixOmin1:ixOmax1, psi_)**2
    lambda(ixOmin1:ixOmax1,2) = - lambda(ixOmin1:ixOmax1,1)
    lambda(ixOmin1:ixOmax1,1) = prim(ixOmin1:ixOmax1,&
        alp_) * lambda(ixOmin1:ixOmax1,1) - prim(ixOmin1:ixOmax1, beta(idim))
    lambda(ixOmin1:ixOmax1,2) = prim(ixOmin1:ixOmax1,&
        alp_) * lambda(ixOmin1:ixOmax1,2) - prim(ixOmin1:ixOmax1, beta(idim))
    cmax(ixOmin1:ixOmax1) = max(abs(lambda(ixOmin1:ixOmax1,1)),&
        abs(lambda(ixOmin1:ixOmax1,2)))
  end subroutine gremhd_get_cmax

  subroutine gremhd_get_csound2(ixImin1,ixImax1,ixOmin1,ixOmax1,prim,csound2)
    use mod_global_parameters
    use mod_eos

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: prim(ixImin1:ixImax1,1:nprim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1)
    integer                         :: idir, ix1

    ! Note: the input are the prim variables
    do ix1 = ixOmin1,ixOmax1 
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix1),prim(ix1, rho_),prim(ix1, eps_))
    enddo
  end subroutine gremhd_get_csound2

  subroutine gremhd_get_lfac2(ps_in,ixImin1,ixImax1,ixOmin1,ixOmax1,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixImin1:ixImax1)
    integer                         :: idir, ix1
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma(ixImin1:ixImax1,1:3,1:3))
    do idir = 1, ndir
       gamma(ixOmin1:ixOmax1,idir,idir) = gamma(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, psi_)**4 
    end do
    
    lfac2 = 1.0d0
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixOmin1:ixOmax1) = lfac2(ixOmin1:ixOmax1) + gamma(ixOmin1:ixOmax1,&
          idir,idir)*prim(ixOmin1:ixOmax1, W_vel(idir))**2
    end do
    end associate
  end subroutine gremhd_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine gremhd_get_tilde_S(ps_in,ixImin1,ixImax1,ixOmin1,ixOmax1,tilde_S)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixImin1:ixImax1)

    double precision                :: E2(ixImin1:ixImax1),&
        B2(ixImin1:ixImax1)
    double precision                :: v2(ixImin1:ixImax1)
    double precision                :: h(ixImin1:ixImax1)    ! enthalpy
    double precision                :: lfac(ixImin1:ixImax1) ! Lorentz factor
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call gremhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3), lfac=lfac(ixImin1:ixImax1),&
        v2=v2(ixImin1:ixImax1), E2=E2(ixImin1:ixImax1), B2=B2(ixImin1:ixImax1),&
        h=h(ixImin1:ixImax1) )

    tilde_S(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1,&
        rho_) * h(ixOmin1:ixOmax1) * lfac(ixOmin1:ixOmax1) * &
       v2(ixOmin1:ixOmax1) + 3.0d0 * prim(ixOmin1:ixOmax1,&
        press_) + 0.5d0 * ( E2(ixOmin1:ixOmax1) + B2(ixOmin1:ixOmax1) )
    tilde_S(ixOmin1:ixOmax1) = tilde_S(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
        psi_)**6
    end associate
  end subroutine gremhd_get_tilde_S

  subroutine gremhd_get_dt(w, ixImin1,ixImax1, ixOmin1,ixOmax1, dtnew, dx1, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1, 1:1)
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nprim)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine gremhd_get_dt

  !> for processing after the step, but before con2prim
  subroutine gremhd_step_adv_global(iit,qt)
    use mod_global_parameters
    use mod_gremhd_phys_divb_mg
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    if (type_divb == divb_multigrid) then
       call gremhd_clean_divb_multigrid()
    end if
  end subroutine gremhd_step_adv_global

  !> Returns 0 in argument flag where values are ok
  subroutine gremhd_check_prim(ixImin1,ixImax1, ixOmin1,ixOmax1, prim, flag)
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
  end subroutine gremhd_check_prim

end module mod_gremhd_phys
