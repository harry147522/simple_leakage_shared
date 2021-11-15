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

    namelist /grmhd_list/  use_GR, tolerance, iter_max, &
                           type_divb_fix, divb_4thorder, &
                           divB_glm_kappa, &
                           divB_mg_n_cycle, divB_mg_redblack, divB_mg_tol, divB_mg_it_max, &
                           type_ct

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
       if (mype == 0) write(*,*) "type_divb_fix is not specificied, now we assume type_divb = none!"
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
       Bvec(idir) = var_set_primvar('Bvec',idir, var_type=hydro_var,need_bc=.false.)
    end do
    if (type_divb == divb_GLM) then
       ! introduce a scalar field Phi for divergence cleaning
       Bphi_ = var_set_primvar('Bphi', var_type=hydro_var,need_bc=.false.)
    end if

    ! Determine primitive variables that needed not to be reconstructed
    press_ = var_set_primvar('press', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    cs2_ = var_set_primvar('cs2', var_type=hydro_var,need_bc=.false.,need_rec=.false.)

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
  subroutine grmhd_update_eos(ixI^L, ixO^L, prim)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(inout)           :: prim(ixI^S, 1:nprim)
    integer                                   :: ix^D

    ! rho and eps are given, update the rest of the primitive variables
    {do ix^D = ixO^LIM^D \}
       call grmhd_update_eos_one_point(prim(ix^D, :))
    {enddo^D&\}
  end subroutine grmhd_update_eos

  !> Calculate cmax_idim within ixO^L
  subroutine grmhd_get_cmax(prim, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: prim(ixI^S, nprim), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    integer                                   :: idir, ix^D
    double precision                          :: lambda(ixI^S,1:3)
    call grmhd_get_lambda(ixI^L, ixO^L, idim, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambda(ixI^S,1:2))
    cmax(ixO^S) = max(abs(lambda(ixO^S,1)), abs(lambda(ixO^S,2)))
  end subroutine grmhd_get_cmax

  subroutine grmhd_get_csound2(ixI^L,ixO^L,prim,csound2)
    use mod_global_parameters
    use mod_eos

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: prim(ixI^S,1:nprim)
    double precision, intent(out)   :: csound2(ixI^S)
    integer                         :: idir, ix^D

    ! Note: the input are the prim variables
    {do ix^D = ixO^LIM^D \}
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix^D),prim(ix^D, rho_),prim(ix^D, eps_))
       ! strictly require cs2 is physical
       csound2(ix^D) = max( min( 1.0d0, csound2(ix^D) ) , 0.0d0)
    {enddo^D&\}
  end subroutine grmhd_get_csound2

  subroutine grmhd_get_lfac2(ps_in,ixI^L,ixO^L,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixI^S)
    integer                         :: idir, ix^D
    double precision                :: gamma(ixI^S,1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    
    lfac2 = 1.0d0
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixO^S) = lfac2(ixO^S) + gamma(ixO^S,idir,idir)*prim(ixO^S, W_vel(idir))**2
    end do
    end associate
  end subroutine grmhd_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine grmhd_get_tilde_S(ps_in,ixI^L,ixO^L,tilde_S)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixI^S)

    double precision                :: B_dot_v(ixI^S)   
    double precision                :: b2(ixI^S) 
    double precision                :: htot(ixI^S)    ! modified enthalpy:(h + b2/rho)
    double precision                :: Ptot(ixI^S) ! total pressure
    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call grmhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), &
                B_dot_v=B_dot_v(ixI^S), b2=b2(ixI^S), &
                Ptot=Ptot(ixI^S), &
                htot=htot(ixI^S) )

    ! nota that S = S_i v^i + 3 p* - b^2
    tilde_S(ixO^S) = 3.0d0 * Ptot(ixO^S) - b2(ixO^S)
    do idir = 1, ndir
       tilde_S(ixO^S) = tilde_S(ixO^S) &
        + ( ( prim(ixO^S, rho_) * htot(ixO^S) - B_dot_v(ixO^S)**2 ) * prim(ixO^S, W_vel(idir)) &
                  - B_dot_v(ixO^S) * prim(ixO^S, Bvec(idir)) / lfac(ixO^S) ) & 
                 * gamma(ixO^S,idir,idir) * prim(ixO^S, W_vel(idir)) 

    end do
    tilde_S(ixO^S) = tilde_S(ixO^S) * prim(ixO^S, psi_)**6
    end associate
  end subroutine grmhd_get_tilde_S

  subroutine grmhd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nprim)
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
  subroutine grmhd_check_prim(ixI^L, ixO^L, prim, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: prim(ixI^S, 1:nprim)
    integer, intent(inout)       :: flag(ixI^S)
    flag(ixO^S) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where(prim(ixO^S, rho_) < smalldouble ) 
       flag(ixO^S) = rho_
    else where(prim(ixO^S, eps_) < smalldouble ) 
       flag(ixO^S) = eps_
    end where
  end subroutine grmhd_check_prim

  !> This subroutine fix the abnormal values in primitive variables !
  subroutine grmhd_handle_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
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
             call grmhd_update_eos_one_point(prim(ix^D,1:nprim))
          end if
       {enddo^D&\}
    case default
       ! nothing to do here
       return
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
    end select
  end subroutine grmhd_handle_small_values

  subroutine grmhd_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nprim), wRC(ixI^S,1:nprim)
    double precision, intent(inout) :: wLp(ixI^S,1:nprim), wRp(ixI^S,1:nprim)
    type(state)                     :: s
    double precision                :: dB(ixI^S), dPsi(ixI^S)

    if(stagger_grid) then
      wLC(ixO^S,Bvec(idir))=s%prims(ixO^S,idir)
      wRC(ixO^S,Bvec(idir))=s%prims(ixO^S,idir)
      wLp(ixO^S,Bvec(idir))=s%prims(ixO^S,idir)
      wRp(ixO^S,Bvec(idir))=s%prims(ixO^S,idir)
    end if

    !if(associated(usr_set_wLR)) call usr_set_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine grmhd_modify_wLR

end module mod_grmhd_phys
