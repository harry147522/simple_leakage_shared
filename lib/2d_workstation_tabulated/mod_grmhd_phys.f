!> General relativistic hydrodynamics physics module (only in CFC now)
module mod_grmhd_phys
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private


  ! Public methods
  public :: grmhd_phys_init
  public :: grmhd_get_divB

contains

  !> Read this module's parameters from a file
  subroutine grmhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    character(len=name_len)      :: type_divb_fix = ""

    namelist /grmhd_list/  use_GR, tolerance, iter_max, type_divb_fix,&
        divb_4thorder, divB_glm_kappa, divB_mg_n_cycle, divB_mg_redblack,&
        divB_mg_tol, divB_mg_it_max, type_ct

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grmhd_list, end=111)
111    close(unitpar)
    end do

    select case (type_divb_fix)
    case ('none')
       type_divb = none
    case ('multigrid')
       type_divb = divb_multigrid
    case ('glm')
       type_divb = divb_GLM
    case ('ct')
      type_divb = divb_ct
    case default
       if (mype == 0) write(*,*)&
           "type_divb_fix is not specificied, now we assume type_divb = none!"
       type_divb = none
    end select

  end subroutine grmhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine grmhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!grmhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine grmhd_write_info


  !> Initialize the module
  subroutine grmhd_phys_init()
    use mod_global_parameters
    use mod_grmhd_phys_convert
    use mod_grmhd_phys_flux
    use mod_grmhd_phys_add_source
    use mod_grmhd_phys_divb_mg
    use mod_grmhd_phys_divb_ct

    integer :: itr, idir

    call grmhd_read_params(par_files)

    physics_type = "grmhd"
  
    ! Determine primitive variables that needed to be reconstructed
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)

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
    allocate(Bcons(ndir))
    do idir = 1, ndir
       Bcons(idir) = var_set_consvar('Bcons',idir, var_type=hydro_var)
    end do
    if (type_divb == divb_GLM) then
       ! introduce a scalar field Phi for divergence cleaning
       Bphi_cons_ = var_set_consvar('Bphi_cons', var_type=hydro_var)
    end if

    phys_get_dt              => grmhd_get_dt
    phys_get_lfac2           => grmhd_get_lfac2
    phys_get_tilde_S         => grmhd_get_tilde_S
    phys_get_csound2         => grmhd_get_csound2
    phys_get_cmax            => grmhd_get_cmax
    phys_update_eos          => grmhd_update_eos
    phys_check_params        => grmhd_check_params
    phys_check_prim          => grmhd_check_prim
    phys_write_info          => grmhd_write_info
    phys_handle_small_values => grmhd_handle_small_values

    call grmhd_phys_convert_init()
    call grmhd_phys_flux_init()
    call grmhd_phys_add_source_init()
    
    select case (type_divb)
    case (divb_glm)
       ! nothing
    case (divb_multigrid)
       call grmhd_phys_divb_mg_init()
       ! setup div cleaning here. fixme: maybe need to move it our later
       phys_step_adv_global => grmhd_step_adv_global
    case (divb_ct)
       call grmhd_phys_divb_ct_init()
    case default
       ! nothing
    end select

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, ncons))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, ncons])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 4 ! No. vector vars
    allocate(iw_vector(1:nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1 
    iw_vector(3) = vecX(1) - 1
    iw_vector(4) = Bvec(1) - 1

  end subroutine grmhd_phys_init

  subroutine grmhd_check_params
    use mod_global_parameters
    use mod_eos

    ! check whether eos module has been loaded
    call eos_check

  end subroutine grmhd_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grmhd_update_eos(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, prim)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout)           :: prim(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nprim)
    integer                                   :: ix1,ix2

    ! rho and eps are given, update the rest of the primitive variables
    do ix1 = ixOmin1,ixOmax1 
    do ix2 = ixOmin2,ixOmax2 
       call grmhd_update_eos_one_point(prim(ix1,ix2, :))
    enddo
    enddo
  end subroutine grmhd_update_eos

  !> Calculate cmax_idim within ixO^L
  subroutine grmhd_get_cmax(prim, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: prim(ixImin1:ixImax1,&
       ixImin2:ixImax2, nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                                   :: idir, ix1,ix2
    double precision                          :: lambda(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3)
    call grmhd_get_lambda(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, prim(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim),&
        x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), lambda(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:2))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(abs(lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)), abs(lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
  end subroutine grmhd_get_cmax

  subroutine grmhd_get_csound2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,prim,csound2)
    use mod_global_parameters
    use mod_eos

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
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
  end subroutine grmhd_get_csound2

  subroutine grmhd_get_lfac2(ps_in,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
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
    call get_gamma_ij_hat(x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim), ixImin1,&
       ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
        gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3))
    do idir = 1, ndir
       gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) = gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, psi_)**4 
    end do
    
    lfac2 = 1.0d0
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = lfac2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) + gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir)*prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, W_vel(idir))**2
    end do
    end associate
  end subroutine grmhd_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine grmhd_get_tilde_S(ps_in,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,tilde_S)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    double precision                :: B_dot_v(ixImin1:ixImax1,&
       ixImin2:ixImax2)   
    double precision                :: b2(ixImin1:ixImax1,ixImin2:ixImax2) 
    double precision                :: htot(ixImin1:ixImax1,ixImin2:ixImax2) !modified enthalpy:(h + b2/rho)
    double precision                :: Ptot(ixImin1:ixImax1,ixImin2:ixImax2) !total pressure
    double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2) !Lorentz factor
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call grmhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        lfac=lfac(ixImin1:ixImax1,ixImin2:ixImax2),&
        B_dot_v=B_dot_v(ixImin1:ixImax1,ixImin2:ixImax2),&
        b2=b2(ixImin1:ixImax1,ixImin2:ixImax2), Ptot=Ptot(ixImin1:ixImax1,&
       ixImin2:ixImax2), htot=htot(ixImin1:ixImax1,ixImin2:ixImax2) )

    ! nota that S = S_i v^i + 3 p* - b^2
    tilde_S(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3.0d0 * Ptot(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    do idir = 1, ndir
       tilde_S(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tilde_S(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) + ( ( prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_) * htot(ixOmin1:ixOmax1,ixOmin2:ixOmax2) - &
          B_dot_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2 ) * prim(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, W_vel(idir)) - B_dot_v(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir)) / lfac(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) ) * gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
          idir) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, W_vel(idir)) 

    end do
    tilde_S(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tilde_S(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, psi_)**6
    end associate
  end subroutine grmhd_get_tilde_S

  subroutine grmhd_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
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

  end subroutine grmhd_get_dt

  !> for processing after the step, but before con2prim
  subroutine grmhd_step_adv_global(iit,qt)
    use mod_global_parameters
    use mod_grmhd_phys_divb_mg
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    if (type_divb == divb_multigrid) then
       call grmhd_clean_divb_multigrid()
    end if
  end subroutine grmhd_step_adv_global

  !> Returns 0 in argument flag where values are ok
  subroutine grmhd_check_prim(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, prim, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: prim(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    integer, intent(inout)       :: flag(ixImin1:ixImax1,ixImin2:ixImax2)
    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where(prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) < smalldouble ) 
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
    else where(prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, eps_) < smalldouble ) 
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = eps_
    end where
  end subroutine grmhd_check_prim

  !> This subroutine fix the abnormal values in primitive variables !
  subroutine grmhd_handle_small_values(prim, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
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

    ! fixme: delete small_values_fix_iw

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
             prim(ix1,ix2, W_vel) = 0.0d0
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
             call grmhd_update_eos_one_point(prim(ix1,ix2,1:nprim))
          end if
       enddo
       enddo
    case default
       ! nothing to do here
       return
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
    end select
  end subroutine grmhd_handle_small_values

  subroutine grmhd_modify_wLR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)
    double precision, intent(inout) :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim), wRp(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)
    type(state)                     :: s
    double precision                :: dB(ixImin1:ixImax1,ixImin2:ixImax2),&
        dPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    if(stagger_grid) then
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Bvec(idir))=s%prims(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Bvec(idir))=s%prims(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Bvec(idir))=s%prims(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Bvec(idir))=s%prims(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
    end if

    !if(associated(usr_set_wLR)) call usr_set_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine grmhd_modify_wLR

end module mod_grmhd_phys
