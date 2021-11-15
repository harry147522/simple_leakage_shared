!> Module with all the methods that users can customize in GMUNU
!>
!> Each procedure pointer can be initialized in a user's mod_usr.t
module mod_usr_methods

  implicit none
  public

  !> Initialize the user's settings (after initializing gmunu)
  procedure(p_no_args), pointer       :: usr_set_parameters   => null()
  !> Initialize earch grid block data
  procedure(init_one_grid), pointer   :: usr_init_one_grid    => null()

  ! Boundary condition related
  procedure(special_bc), pointer      :: usr_special_bc       => null()
  procedure(internal_bc), pointer     :: usr_internal_bc      => null()

  ! Output related
  procedure(p_no_args), pointer       :: usr_print_log        => null()
  procedure(p_no_args), pointer       :: usr_write_analysis   => null()
  procedure(transform_w), pointer     :: usr_transform_w      => null()
  procedure(aux_output), pointer      :: usr_aux_output       => null()
  procedure(sub_modify_io), pointer   :: usr_modify_output    => null()
  procedure(special_convert), pointer :: usr_special_convert  => null()

  ! Called at the beginning of every time step (after determining dt)
  procedure(process_grid), pointer    :: usr_process_grid     => null()
  procedure(process_global), pointer  :: usr_process_global   => null()

  ! Called every time step just after advance (with w^(n+1), it^n, global_time^n)
  procedure(process_adv_grid), pointer   :: usr_process_adv_grid   => null()
  procedure(process_adv_global), pointer :: usr_process_adv_global => null()

  ! Called after initial condition before the start of the simulation
  procedure(p_no_args), pointer       :: usr_improve_initial_condition => &
     null()

  ! Called before the start of the simulation
  procedure(p_no_args), pointer       :: usr_before_main_loop => null()

  ! Source terms
  procedure(source), pointer          :: usr_source           => null()
  procedure(get_dt), pointer          :: usr_get_dt           => null()

  ! Refinement related procedures
  procedure(refine_grid), pointer     :: usr_refine_grid      => null()
  procedure(var_for_errest), pointer  :: usr_var_for_errest   => null()
  procedure(a_refine_threshold), pointer :: usr_refine_threshold => null()
  procedure(flag_grid), pointer       :: usr_flag_grid        => null()

  ! Called after the mesh has been adjuste
  procedure(after_refine), pointer      :: usr_after_refine => null()

  ! allow user to explicitly set flux at cell interfaces for finite volume scheme
  procedure(set_flux), pointer      :: usr_set_flux => null()

  ! initialize vector potential on cell edges for magnetic field
  procedure(init_vector_potential), pointer :: usr_init_vector_potential => &
     null()

  ! allow user to define conductivity or resistivity
  procedure(get_resistivity), pointer        :: usr_get_resistivity       => &
     null()
  ! allow user to define the dynamo alpha-term
  procedure(get_resistivity), pointer        :: usr_get_dynamo_coeff      => &
     null()

  ! allow user to define the frequency-integrated emissivity
  procedure(get_opacities), pointer          :: usr_get_opacities         => &
     null()

  abstract interface

    subroutine p_no_args()
    end subroutine p_no_args

    !> Initialize one grid
    subroutine init_one_grid(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
    end subroutine init_one_grid

     !> special boundary types, users must assign conservative
     !> variables in boundaries
     subroutine special_bc(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iB,w,x)
       use mod_global_parameters
       !> Shape of input arrays
       integer, intent(in)             :: ixImin1,ixImax1
       !> Region where boundary values have to be set
       integer, intent(in)             :: ixOmin1,ixOmax1
       !> Integer indicating direction of boundary
       integer, intent(in)             :: iB
       double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
       double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
     end subroutine special_bc

    !> internal boundary, user defined
    !> This subroutine can be used to artificially overwrite ALL conservative
    !> variables in a user-selected region of the mesh, and thereby act as
    !> an internal boundary region. It is called just before external (ghost cell)
    !> boundary regions will be set by the BC selection. Here, you could e.g.
    !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
    !> which can be used to identify the internal boundary region location.
    !> Its effect should always be local as it acts on the mesh.
    subroutine internal_bc(level,qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    end subroutine internal_bc

    !> this subroutine is ONLY to be used for computing auxiliary variables
    !> which happen to be non-local (like div v), and are in no way used for
    !> flux computations. As auxiliaries, they are also not advanced
    subroutine process_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,&
       x)
      use mod_global_parameters
      integer, intent(in)             :: igrid,level,ixImin1,ixImax1,ixOmin1,&
         ixOmax1
      double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
    end subroutine process_grid

    !> If defined, this routine is called before writing output, and it can
    !> set/modify the variables in the w array.
    subroutine sub_modify_io(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
    end subroutine sub_modify_io

    !> This subroutine is called at the beginning of each time step
    !> by each processor. No communication is specified, so the user
    !> has to implement MPI routines if information has to be shared
    subroutine process_global(iit,qt)
      use mod_global_parameters
      integer, intent(in)          :: iit
      double precision, intent(in) :: qt
    end subroutine process_global

    !> for processing after the advance (PIC-MHD, e.g.)
    subroutine process_adv_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,&
       w,x)
      use mod_global_parameters
      integer, intent(in)             :: igrid,level,ixImin1,ixImax1,ixOmin1,&
         ixOmax1
      double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
    end subroutine process_adv_grid

    !> for processing after the advance (PIC-MHD, e.g.)
    subroutine process_adv_global(iit,qt)
      use mod_global_parameters
      integer, intent(in)          :: iit
      double precision, intent(in) :: qt
    end subroutine process_adv_global

    !> this subroutine can be used in convert, to add auxiliary variables to the
    !> converted output file, for further analysis using tecplot, paraview, ....
    !> these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    !> the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    !> corresponding normalization values (default value 1)
    subroutine aux_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
      !double precision             :: w(ixI^S,nprim+nwauxio)
      !double precision             :: normconv(0:nprim+nwauxio)
      double precision             :: w(ixImin1:ixImax1,nprim)
      double precision             :: normconv(0:nprim)
    end subroutine aux_output

    !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
    !> iw=iwmin...iwmax.  wCT is at time qCT
    subroutine source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
       primCT,qt,cons,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
          iwmin,iwmax
      double precision, intent(in)    :: qdt, qtC, qt
      double precision, intent(in)    :: primCT(ixImin1:ixImax1,1:nprim),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: cons(ixImin1:ixImax1,1:ncons)
    end subroutine source

    !> Limit "dt" further if necessary, e.g. due to the special source terms.
    !> The getdt_courant (CFL condition) and the getdt subroutine in the GMUNUPHYS
    !> module have already been called.
    subroutine get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: w(ixImin1:ixImax1,1:nprim)
      double precision, intent(inout) :: dtnew
    end subroutine get_dt

    !> Enforce additional refinement or coarsening
    !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    !> you must set consistent values for integers refine/coarsen:
    !> refine = -1 enforce to not refine
    !> refine =  0 doesn't enforce anything
    !> refine =  1 enforce refinement
    !> coarsen = -1 enforce to not coarsen
    !> coarsen =  0 doesn't enforce anything
    !> coarsen =  1 enforce coarsen
    !> e.g. refine for negative first coordinate x < 0 as
    !> if (any(x(ix^S,1) < zero)) refine=1
    subroutine refine_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x,&
       refine,coarsen)
      use mod_global_parameters
      integer, intent(in)          :: igrid, level, ixImin1,ixImax1, ixOmin1,&
         ixOmax1
      double precision, intent(in) :: qt, w(ixImin1:ixImax1,1:nprim),&
          x(ixImin1:ixImax1,1:ndim)
      integer, intent(inout)       :: refine, coarsen
    end subroutine refine_grid

    !> this is the place to compute a local auxiliary variable to be used
    !> as refinement criterion for the Lohner error estimator only
    !>  -->it is then requiring and iflag>nw
    !> note that ixO=ixI=ixG, hence the term local (gradients need special attention!)
    subroutine var_for_errest(ixImin1,ixImax1,ixOmin1,ixOmax1,iflag,w,x,var)
      use mod_global_parameters
      integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1,iflag
      double precision, intent(in)  :: w(ixImin1:ixImax1,1:nprim),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out) :: var(ixImin1:ixImax1)
    end subroutine var_for_errest

    !> regenerate w and eqpartf arrays to output into *tf.dat
    subroutine transform_w(ixImin1,ixImax1,ixOmin1,ixOmax1,nw_in,w_in,x,w_out)
      use mod_global_parameters
      integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, nw_in
      double precision, intent(in)  :: w_in(ixImin1:ixImax1,1:nw_in)
      double precision, intent(in)  :: x(ixImin1:ixImax1, 1:ndim)
      double precision, intent(out) :: w_out(ixImin1:ixImax1,1:nprim)
    end subroutine transform_w

    !> use different threshold in special regions for AMR to
    !> reduce/increase resolution there where nothing/something interesting happens.
    subroutine a_refine_threshold(wlocal,xlocal,threshold,qt)
      use mod_global_parameters
      double precision, intent(in)    :: wlocal(1:nprim),xlocal(1:ndim),qt
      double precision, intent(inout) :: threshold
    end subroutine a_refine_threshold

    !> Allow user to use their own data-postprocessing procedures
    subroutine special_convert(qunitconvert)
      use mod_global_parameters
      integer, intent(in) :: qunitconvert
      character(len=20)   :: userconvert_type
    end subroutine special_convert

    !> flag=-1 : Treat all cells active, omit deactivation (onentry, default)
    !> flag=0  : Treat as normal domain
    !> flag=1  : Treat as passive, but reduce by safety belt
    !> flag=2  : Always treat as passive
    subroutine flag_grid(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,flag)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      integer, intent(inout)          :: flag
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    end subroutine flag_grid

    subroutine after_refine(n_coarsen, n_refine)
      integer, intent(in) :: n_coarsen
      integer, intent(in) :: n_refine
    end subroutine after_refine

    !> allow user to explicitly set flux at cell interfaces for finite volume scheme
    subroutine set_flux(ixImin1,ixImax1,ixCmin1,ixCmax1,idim,fC,x)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixCmin1,ixCmax1, idim
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      ! face-center flux
      double precision,intent(inout) :: fC(ixImin1:ixImax1,1:ncons,1:ndim)
      ! For example, to set flux at bottom boundary in a 3D box for induction equation
      ! vobs and bobs are interpolated data from original observational data for data-driven application
      !integer :: idir
      !if(idim==3) then
      !  if(block%is_physical_boundary(idim*2-1)) then
      !    do idir=1,ndir
      !       fC(ixCmin3^%3ixC^S,mag(idir),idim)=vobs(ixCmin3+1^%3ixC^S,idim)*bobs(ixCmin3+1^%3ixC^S,idir)-vobs(ixCmin3+1^%3ixC^S,idir)*bobs(ixCmin3+1^%3ixC^S,idim)
      !    end do
      !  end if
      !end if
    end subroutine set_flux

    !> initialize vector potential on cell edges for magnetic field
    subroutine init_vector_potential(ixImin1,ixImax1, ixCmin1,ixCmax1, xC, A,&
        idir)
      use mod_global_parameters

      integer, intent(in)                :: ixImin1,ixImax1, ixCmin1,ixCmax1,&
          idir
      double precision, intent(in)       :: xC(ixImin1:ixImax1,1:ndim)
      double precision, intent(out)      :: A(ixImin1:ixImax1)

    end subroutine init_vector_potential

    !> allow user to define conductivity or resistivity (eta)
    subroutine get_resistivity(ixImin1,ixImax1,ixOmin1,ixOmax1,cons,x,eta)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: cons(ixImin1:ixImax1,1:ncons)
      double precision, intent(out)   :: eta(ixImin1:ixImax1)
    end subroutine get_resistivity 

    !> allow user to define frequency-integrated emissivity, averaged absorption opacity, and scattering opacity.
    subroutine get_opacities(ixImin1,ixImax1,ixOmin1,ixOmax1,cons,prim,x,&
        var_out)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: cons(ixImin1:ixImax1,1:ncons)
      double precision, intent(in)    :: prim(ixImin1:ixImax1,1:nprim)
      double precision, intent(out)   :: var_out(ixImin1:ixImax1,1:3)
    end subroutine get_opacities 

  end interface

end module mod_usr_methods
