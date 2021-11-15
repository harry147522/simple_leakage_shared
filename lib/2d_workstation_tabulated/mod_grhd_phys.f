!> General relativistic hydrodynamics physics module (only in CFC now)
module mod_grhd_phys
  use mod_physics
  use mod_grhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grhd_phys_init

contains

  !> Read this module's parameters from a file
  subroutine grhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grhd_list/ evolve_hydro,  oneDcore, r_core, use_GR, tolerance,&
        iter_max, lfac_max

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine grhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine grhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!grhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine grhd_write_info


  !> Initialize the module
  subroutine grhd_phys_init()
    use mod_global_parameters
    use mod_grhd_phys_convert
    use mod_grhd_phys_flux
    use mod_grhd_phys_add_source
    use mod_grhd_phys_one_d_core

    integer :: itr, idir

    call grhd_read_params(par_files)

    physics_type = "grhd"
  
    ! Determine primitive variables that needed to be reconstructed
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)

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

    phys_get_dt              => grhd_get_dt
    phys_get_lfac2           => grhd_get_lfac2
    phys_get_tilde_S         => grhd_get_tilde_S
    phys_get_csound2         => grhd_get_csound2
    phys_get_cmax            => grhd_get_cmax
    phys_update_eos          => grhd_update_eos
    phys_check_params        => grhd_check_params
    phys_check_prim          => grhd_check_prim
    phys_write_info          => grhd_write_info
    phys_handle_small_values => grhd_handle_small_values

    call grhd_phys_convert_init()
    call grhd_phys_flux_init()
    call grhd_phys_add_source_init()
     call grhd_phys_one_d_core_init()

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
    allocate(iw_vector(nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1
    iw_vector(3) = vecX(1) - 1

  end subroutine grhd_phys_init

  subroutine grhd_check_params
    use mod_global_parameters
    use mod_geometry, only: coordinate, spherical
    use mod_eos

    ! check whether eos module has been loaded
    call eos_check

    
    if ( oneDcore .and. (mype==0) ) then
      if (.not. internalboundary) then
         write(*,*) "internalboundary must be true if use 1D core"
         write(*,*) "now set internalboundary = .true."
         internalboundary = .True.
      end if
      if (coordinate/=spherical) call mpistop(&
         "1D core support only spherical coordinate.")
      if (r_core <= 0.0d0 ) call mpistop("r_core must be larger than zero.")
    end if
   

  end subroutine grhd_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grhd_update_eos(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, prim)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout)           :: prim(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nprim)
    integer                                   :: ix1,ix2

    ! rho and eps are given, update the rest of the primitive variables
    do ix1 = ixOmin1,ixOmax1 
    do ix2 = ixOmin2,ixOmax2 
       call grhd_update_eos_one_point(prim(ix1,ix2,:))
    enddo
    enddo
  end subroutine grhd_update_eos

  !> Calculate cmax_idim within ixO^L
  subroutine grhd_get_cmax(prim, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: prim(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                                   :: idir, ix1,ix2
    double precision                          :: lambda(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:2)
    call grhd_get_lambda(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, prim(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim),&
        x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), lambda(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:2))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(abs(lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)), abs(lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))

    
    ! 1D core treatment
    where ( (idim/=1).and.(oneDcore).and.(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1) <= r_core) )
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
    end where
   
  end subroutine grhd_get_cmax

  subroutine grhd_get_csound2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,prim,csound2)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
       nprim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                         :: idir, ix1,ix2
    ! Note: the input are the prim variables
    do ix1 = ixOmin1,ixOmax1 
    do ix2 = ixOmin2,ixOmax2 
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix1,ix2),prim(ix1,ix2, rho_),prim(ix1,&
          ix2, eps_))
       ! strictly require cs2 is physical
       csound2(ix1,ix2) = max( min( 1.0d0, csound2(ix1,ix2) ) , 0.0d0)
    enddo
    enddo
  end subroutine grhd_get_csound2

  subroutine grhd_get_lfac2(ps_in,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir, ix1,ix2
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2), ixImin1,&
       ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
        gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3))
    do idir = 1, ndir
       gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) = gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, psi_)**4 
    end do
    
    lfac2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.0d0 
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = lfac2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) + gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, W_vel(idir))**2
    end do
    end associate
  end subroutine grhd_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine grhd_get_tilde_S(ps_in,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,tilde_S)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    double precision                :: h(ixImin1:ixImax1,ixImin2:ixImax2)   
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call grhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        h=h(ixImin1:ixImax1,ixImin2:ixImax2) )

    ! nota that S = S_i v^i + 3 p
    tilde_S(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3.0d0 * prim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, press_)
    do idir = 1, ndir
       tilde_S(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tilde_S(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) + ( ( prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_) * h(ixOmin1:ixOmax1,ixOmin2:ixOmax2) ) * prim(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, W_vel(idir)) ) * gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,idir,idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           W_vel(idir))

    end do
    tilde_S(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tilde_S(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, psi_)**6
    end associate
  end subroutine grhd_get_tilde_S

  subroutine grhd_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:2)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine grhd_get_dt

  !> Returns 0 in argument flag where values are ok
  subroutine grhd_check_prim(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, prim, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    integer, intent(inout)       :: flag(ixImin1:ixImax1,ixImin2:ixImax2)
    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where( (prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_) < smalldouble) .or. isnan(prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_)) ) 
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
    else where( (prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        eps_) < smalldouble) .or. isnan(prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        eps_)) ) 
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = eps_
    end where
  end subroutine grhd_check_prim
end module mod_grhd_phys
