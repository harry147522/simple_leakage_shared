!> General relativistic hydrodynamics physics module (only in CFC now)
module mod_grhd_ccsn_leakage_phys
  use mod_physics
  use mod_grhd_ccsn_leakage_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grhd_ccsn_leakage_phys_init

contains

  !> Read this module's parameters from a file
  subroutine grhd_ccsn_leakage_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grhd_ccsn_leakage_list/ evolve_hydro, {^NOONED oneDcore, r_core, } &
                           use_GR, tolerance, iter_max, lfac_max

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grhd_ccsn_leakage_list, end=111)
111    close(unitpar)
    end do

  end subroutine grhd_ccsn_leakage_read_params

  !> Write this module's parameters to a snapsoht
  subroutine grhd_ccsn_leakage_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!grhd_ccsn_leakage_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine grhd_ccsn_leakage_write_info


  !> Initialize the module
  subroutine grhd_ccsn_leakage_phys_init()
    use mod_global_parameters
    use mod_grhd_ccsn_leakage_phys_convert
    use mod_grhd_ccsn_leakage_phys_flux
    use mod_grhd_ccsn_leakage_phys_add_source
    use mod_grhd_ccsn_leakage_phys_one_d_core

    integer :: itr, idir, i, neu

    call grhd_ccsn_leakage_read_params(par_files)

    physics_type = "grhd_ccsn_leakage"
  
    ! Determine primitive variables that needed to be reconstructed
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    ye_ = var_set_primvar('ye', var_type=hydro_var)
    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)

    ! Determine primitive variables that needed not to be reconstructed
    press_ = var_set_primvar('press', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    cs2_ = var_set_primvar('cs2', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    temp_ = var_set_primvar('temp', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    ent_ = var_set_primvar('ent', var_type=hydro_var,need_bc=.false.,need_rec=.false.)

    allocate(leaktau(3))
    do neu = 1, 3
       leaktau(neu) = var_set_primvar('leaktau', neu, var_type=hydro_var, need_bc=.false., need_rec=.false.)
    enddo

    allocate(chi(3))
    do neu = 1,3 
       chi(neu) = var_set_primvar('chi', neu, var_type=hydro_var, need_bc=.false., need_rec=.false.)
    enddo
    allocate(zeta(3))
    do neu = 1,3 
       zeta(neu) = var_set_primvar('zeta', neu, var_type=hydro_var, need_bc=.false., need_rec=.false.)
    enddo


    ! Set index of metric variables, it has to be reconstructed so use fluxvar
    alp_ = var_set_primvar("alp", var_type=metric_var, need_bc=.False.,need_rec=.false.)
    psi_ = var_set_primvar("psi", var_type=metric_var, need_bc=.False.,need_rec=.false.)
    allocate(beta(ndir))
    do idir=1,ndir
       beta(idir) = var_set_primvar("beta", idir, var_type=metric_var, need_bc=.False.,need_rec=.false.)
    enddo
    allocate(vecX(ndir))
    do idir=1,ndir
       vecX(idir) = var_set_primvar("X", idir, var_type=metric_var, need_bc=.False.,need_rec=.false.)
    enddo

    ! Determine hydro flux variables
    D_ = var_set_consvar('D', var_type=hydro_var)
    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = var_set_consvar('mom',idir, var_type=hydro_var)
    end do
    tau_ = var_set_consvar('tau', var_type=hydro_var)
    ye_con_ = var_set_consvar('ye_con', var_type=hydro_var)

   
!  idir is mom(idir), idir+1 is tau, idir+2 is ye _con
    allocate(coolingsource(ndir+2))
    do idir = 1, ndir+2
      coolingsource(idir) = var_set_consvar('coolingsource', idir, var_type=hydro_var)
    enddo

    phys_get_dt              => grhd_ccsn_leakage_get_dt
    phys_get_lfac2           => grhd_ccsn_leakage_get_lfac2
    phys_get_tilde_S         => grhd_ccsn_leakage_get_tilde_S
    phys_get_csound2         => grhd_ccsn_leakage_get_csound2
    phys_get_cmax            => grhd_ccsn_leakage_get_cmax
    phys_update_eos          => grhd_ccsn_leakage_update_eos
    phys_check_params        => grhd_ccsn_leakage_check_params
    phys_check_prim          => grhd_ccsn_leakage_check_prim
    phys_write_info          => grhd_ccsn_leakage_write_info
    phys_handle_small_values => grhd_ccsn_leakage_handle_small_values

    call grhd_ccsn_leakage_phys_convert_init()
    call grhd_ccsn_leakage_phys_flux_init()
    call grhd_ccsn_leakage_phys_add_source_init()
    {^NOONED call grhd_ccsn_leakage_phys_one_d_core_init() }

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
       flux_type(:, ye_con_) = flux_nul
       flux_type(:, coolingsource(:)) = flux_nul
       
    end if

    nvector      = 3 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1
    iw_vector(3) = vecX(1) - 1

  end subroutine grhd_ccsn_leakage_phys_init

  subroutine grhd_ccsn_leakage_check_params
    use mod_global_parameters
    use mod_geometry, only: coordinate, spherical
    use mod_eos

    ! check whether eos module has been loaded
    call eos_check

    {^NOONED
    if ( oneDcore .and. (mype==0) ) then
      if (.not. internalboundary) then
         write(*,*) "internalboundary must be true if use 1D core"
         write(*,*) "now set internalboundary = .true."
         internalboundary = .True.
      end if
      if (coordinate/=spherical) call mpistop("1D core support only spherical coordinate.")
      if (r_core <= 0.0d0 ) call mpistop("r_core must be larger than zero.")
    end if
    }

  end subroutine grhd_ccsn_leakage_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grhd_ccsn_leakage_update_eos(ixI^L, ixO^L, prim)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(inout)           :: prim(ixI^S, 1:nprim)
    integer                                   :: ix^D

    ! rho and eps are given, update the rest of the primitive variables
    {do ix^D = ixO^LIM^D \}
       call grhd_ccsn_leakage_update_eos_one_point(prim(ix^D,:))
    {enddo^D&\}
  end subroutine grhd_ccsn_leakage_update_eos

  !> Calculate cmax_idim within ixO^L
  subroutine grhd_ccsn_leakage_get_cmax(prim, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    integer                                   :: idir, ix^D
    double precision                          :: lambda(ixI^S,1:2)
    call grhd_ccsn_leakage_get_lambda(ixI^L, ixO^L, idim, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambda(ixI^S,1:2))
    cmax(ixO^S) = max(abs(lambda(ixO^S,1)), abs(lambda(ixO^S,2)))

    {^NOONED
    ! 1D core treatment
    where ( (idim/=1).and.(oneDcore).and.(x(ixO^S,1) <= r_core) )
      cmax(ixO^S) = 0.0d0
    end where
    }
  end subroutine grhd_ccsn_leakage_get_cmax

  subroutine grhd_ccsn_leakage_get_csound2(ixI^L,ixO^L,prim,csound2)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: prim(ixI^S,nprim)
    double precision, intent(out)   :: csound2(ixI^S)
    integer                         :: idir, ix^D
    ! Note: the input are the prim variables
    {do ix^D = ixO^LIM^D \}
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix^D),prim(ix^D, rho_),prim(ix^D, eps_), ye=prim(ix^D, ye_), temp=prim(ix^D, temp_))
       ! strictly require cs2 is physical
       csound2(ix^D) = max( min( 1.0d0, csound2(ix^D) ) , 0.0d0)
    {enddo^D&\}
  end subroutine grhd_ccsn_leakage_get_csound2

  subroutine grhd_ccsn_leakage_get_lfac2(ps_in,ixI^L,ixO^L,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixI^S)
    integer                         :: idir, ix^D
    double precision                :: gamma(ixI^S,1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:^ND), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    
    lfac2(ixO^S) = 1.0d0 
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixO^S) = lfac2(ixO^S) + gamma(ixO^S,idir,idir) * prim(ixO^S, W_vel(idir))**2
    end do
    end associate
  end subroutine grhd_ccsn_leakage_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine grhd_ccsn_leakage_get_tilde_S(ps_in,ixI^L,ixO^L,tilde_S)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixI^S)

    double precision                :: h(ixI^S)   
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call grhd_ccsn_leakage_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                h=h(ixI^S) )

    ! nota that S = S_i v^i + 3 p
    tilde_S(ixO^S) = 3.0d0 * prim(ixO^S, press_)
    do idir = 1, ndir
       tilde_S(ixO^S) = tilde_S(ixO^S) &
        + ( ( prim(ixO^S, rho_) * h(ixO^S) ) * prim(ixO^S, W_vel(idir)) ) & 
                 * gamma(ixO^S,idir,idir) * prim(ixO^S, W_vel(idir))

    end do
    tilde_S(ixO^S) = tilde_S(ixO^S) * prim(ixO^S, psi_)**6
    end associate
  end subroutine grhd_ccsn_leakage_get_tilde_S

  subroutine grhd_ccsn_leakage_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nprim)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine grhd_ccsn_leakage_get_dt

  !> Returns 0 in argument flag where values are ok
  subroutine grhd_ccsn_leakage_check_prim(ixI^L, ixO^L, prim, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: prim(ixI^S, 1:nprim)
    integer, intent(inout)       :: flag(ixI^S)
    flag(ixO^S) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where( (prim(ixO^S, rho_) < smalldouble) .or. isnan(prim(ixO^S, rho_)) ) 
       flag(ixO^S) = rho_
    else where( (prim(ixO^S, eps_) < smalldouble) .or. isnan(prim(ixO^S, eps_)) ) 
       flag(ixO^S) = eps_
    end where
  end subroutine grhd_ccsn_leakage_check_prim
end module mod_grhd_ccsn_leakage_phys
