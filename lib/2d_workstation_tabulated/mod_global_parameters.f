!> This module contains definitions of global parameters and variables and some
!> generic functions/subroutines used in GMUNU.
!>
!> \todo Move the parameters to the relevant (physics) modules
module mod_global_parameters
  use mod_physicaldata
  use mod_connectivity
  use mpi
  use mod_constants
  use mod_variables
  use mod_basic_types

  implicit none
  public

  ! Parameters

  character(len=*), parameter :: undefined = 'undefined'

  !> @todo Move mpi related variables to e.g. mod_comm

  !> The number of MPI tasks
  integer :: npe

  !> The rank of the current MPI task
  integer :: mype

  !> The MPI communicator
  integer :: icomm

  !> A global MPI error return code
  !> @todo Make local
  integer :: ierrmpi

  !> MPI file handle for logfile
  integer :: log_fh
  !> MPI type for block including ghost cells and its size
  integer :: type_block, size_block
  !> MPI type for block coarsened by 2, and for its children blocks
  integer :: type_coarse_block, type_sub_block(2,2)
  !> MPI type for staggered block coarsened by 2, and for its children blocks
  integer :: type_coarse_block_stg(2,2,2), type_sub_block_stg(2,2,2)
  !> MPI type for IO: block excluding ghost cells
  integer :: type_block_io, size_block_io
  !> MPI type for IO of staggered variables
  integer :: type_block_io_stg, size_block_io_stg
  !> MPI type for IO: cell corner (xc) or cell center (xcc) coordinates
  integer :: type_block_xc_io,type_block_xcc_io
  !> MPI type for IO: cell corner (wc) or cell center (wcc) variables
  integer :: type_block_wc_io,type_block_wcc_io

  ! geometry and domain setups 

  !> the mesh range (within a block with ghost cells)
  integer :: ixMlo1,ixMlo2,ixMhi1,ixMhi2

  !> minimum and maximum domain boundaries for each dimension
  double precision  :: xprobmin1,xprobmin2,xprobmax1,xprobmax2

  !> Indices for cylindrical and spherical coordinates, negative value when not used:
  integer :: r_ = -1
  integer :: theta_ = -1
  integer :: phi_ = -1
  integer :: z_ = -1

  !> Number of spatial dimensions for grid variables
  integer, parameter :: ndim=2

  !> Number of spatial dimensions (components) for vector variables
  integer :: ndir=ndim

  !> Cartesian geometry or not
  logical :: slab

  !> uniform Cartesian geometry or not (stretched Cartesian)
  logical :: slab_uniform

  !> number of grid blocks in domain per dimension, in array over levels
  integer, dimension(:), allocatable :: ng1,ng2
  !> extent of grid blocks in domain per dimension, in array over levels
  double precision, dimension(:), allocatable :: dg1,dg2

  !> number of cells for each dimension in level-one mesh
  integer :: domain_nx1,domain_nx2

  !> number of cells for each dimension in grid block excluding ghostcells
  integer :: block_nx1,block_nx2

  !> Lower index of grid block arrays (always 1)
  integer, parameter :: ixGlo1 = 1, ixGlo2 = 1

  !> Upper index of grid block arrays
  integer :: ixGhi1,ixGhi2

  !> Lower index of stagger grid block arrays (always 0)
  integer, parameter :: ixGslo1 = 0, ixGslo2 = 0

  !> Upper index of stagger grid block arrays
  integer :: ixGshi1,ixGshi2

  !> Number of ghost cells surrounding a grid
  integer :: nghostcells

  integer, parameter :: stretch_none = 0 !< No stretching
  integer, parameter :: stretch_uni  = 1 !< Unidirectional stretching from a side
  integer, parameter :: stretch_symm = 2 !< Symmetric stretching around the center

  !> If true, adjust mod_geometry routines to account for grid stretching (but
  !> the flux computation will not)
  logical :: stretch_uncentered
  !> True if a dimension is stretched
  logical :: stretched_dim(ndim)
  !> What kind of stretching is used per dimension
  integer :: stretch_type(ndim)
  !> stretch factor between cells at AMR level 1, per dimension
  double precision ::  qstretch_baselevel(ndim)
  !> (even) number of (symmetrically) stretched
  !> blocks at AMR level 1, per dimension
  integer ::  nstretchedblocks_baselevel(ndim)
  !> (even) number of (symmetrically) stretched blocks per level and dimension
  integer, allocatable ::  nstretchedblocks(:,:)
  !> physical extent of stretched border in symmetric stretching
  double precision :: xstretch1,xstretch2
  !> Stretching factors and first cell size for each AMR level and dimension
  double precision, allocatable :: qstretch(:,:), dxfirst(:,:),  dxfirst_1mq(:,&
     :), dxmid(:,:)

  !> grid hierarchy info (level and grid indices)
  integer, parameter :: nodehi=2+1
  integer, parameter :: plevel_=1
  integer, parameter :: pig1_=plevel_+1,pig2_=plevel_+2

  integer, allocatable :: node(:,:)
  integer, allocatable :: node_sub(:,:)

  !> grid location info (corner coordinates and grid spacing)
  integer, parameter :: rnodehi=3*2
  integer, parameter :: rpxmin0_=0
  integer, parameter :: rpxmin1_=rpxmin0_+1,rpxmin2_=rpxmin0_+2
  integer, parameter :: rpxmax0_=2
  integer, parameter :: rpxmax1_=rpxmax0_+1,rpxmax2_=rpxmax0_+2
  integer, parameter :: rpdx1_=2*2+1,rpdx2_=2*2+2

  !> Corner coordinates
  double precision, allocatable :: rnode(:,:)
  double precision, allocatable :: rnode_sub(:,:)

  double precision, allocatable :: dx(:,:)
  double precision :: dxlevel(ndim)

  ! IO related quantities

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100

  !> Number of output methods
  integer, parameter :: nfile = 5

  !> Names of the output methods
  character(len=40), parameter  :: output_names(nfile) = ['log      ',&
      'normal   ', 'slice    ', 'collapsed', 'analysis ']

  !> User parameter file
  character(len=std_len)   :: usr_filename 

  !> If collapse(DIM) is true, generate output integrated over DIM
  logical :: collapse(ndim)

  !> Save output of type N on times tsave(:, N)
  double precision :: tsave(nsavehi,nfile)

  !> \todo Move tsavelast to gmunu.t
  double precision :: tsavelast(nfile)

  !> Repeatedly save output of type N when dtsave(N) simulation time has passed
  double precision :: dtsave(nfile)

  !> Save output of type N on iterations itsave(:, N)
  integer :: itsave(nsavehi,nfile)

  !> \todo remove itsavelast?
  integer :: itsavelast(nfile)

  !> Repeatedly save output of type N when ditsave(N) time steps have passed
  integer :: ditsave(nfile)

  !> \todo Move to gmunu.t
  integer :: isavet(nfile)

  !> \todo Move to gmunu.t
  integer :: isaveit(nfile)

  !> Start of read out (not counting specified read outs)
  double precision :: tsavestart(nfile)

  !> The level at which to produce line-integrated / collapsed output
  integer :: collapseLevel

  !> Number of saved files of each type
  !> \todo Move to mod_input_output
  integer :: n_saves(1:nfile)

  !> whether or not to save an output file
  logical :: save_file(nfile)

  !> to monitor timeintegration loop at given wall-clock time intervals
  double precision :: time_between_print

  !> accumulated wall-clock time spent on boundary conditions
  double precision :: time_bc

  !> IO: snapshot and collapsed views output numbers/labels
  integer :: snapshotnext, collapsenext

  !> Constant indicating log output
  integer, parameter :: filelog_      = 1

  !> Constant indicating regular output
  integer, parameter :: fileout_      = 2

  !> Constant indicating slice output
  integer, parameter :: fileslice_    = 3

  !> Constant indicating collapsed output
  integer, parameter :: filecollapse_ = 4

  !> Constant indicating analysis output (see @ref analysis.md)
  integer, parameter :: fileanalysis_ = 5

  !> Unit for standard input
  integer, parameter :: unitstdin=5

  !> Unit for standard output
  integer, parameter :: unitterm=6

  !> Unit for error messages
  integer, parameter :: uniterr=6

  !> \todo Move to mod_input_output
  integer, parameter :: unitpar=9
  integer, parameter :: unitconvert=10
  integer, parameter :: unitslice=11
  integer, parameter :: unitsnapshot=12
  integer, parameter :: unitcollapse=13
  integer, parameter :: unitanalysis=14

  !> IO switches for conversion
  logical          :: nocartesian
  logical, allocatable :: w_write(:)
  logical, allocatable :: writelevel(:)
  double precision :: writespshift(ndim,2)
  integer          :: level_io, level_io_min, level_io_max

  !> Which par files are used as input
  character(len=std_len), allocatable :: par_files(:)

  !> Base file name for simulation output, which will be followed by a number
  character(len=std_len) :: base_filename

  !> If not 'unavailable', resume from snapshot with this base file name
  character(len=std_len) :: restart_from_file

  !> Which type of log to write: 'normal', 'special', 'regression_test'
  character(len=std_len) :: typefilelog

  !> Resume from the snapshot with this index
  integer :: snapshotini

  !> If true, restart a previous run from the latest snapshot
  logical :: resume_previous_run

  !> If true and restart_from_file is given, convert snapshots to
  !> other file formats
  logical                :: convert

  !> If true, already convert to output format during the run
  logical                :: autoconvert

  !> Which format to use when converting
  !>
  !> Options are: tecplot, tecplotCC, vtu, vtuCC, vtuB, vtuBCC, 
  !> tecplotmpi, tecplotCCmpi, vtumpi, vtuCCmpi, vtuBmpi, vtuBCCmpi, pvtumpi, pvtuCCmpi,
  !> pvtuBmpi, pvtuBCCmpi, tecline, teclinempi, onegrid
  character(len=std_len) :: convert_type

  character(len=std_len) :: collapse_type

  !> Conversion factors the primitive variables
  double precision, allocatable :: w_convert_factor(:)

  double precision :: length_convert_factor

  !> Conversion factor for time unit
  double precision       :: time_convert_factor

  integer                :: saveigrid

  !> Stores the memory and load imbalance, used in printlog
  double precision       :: Xload, Xmemory

  !> Save a snapshot before crash a run met unphysical values
  logical :: crash=.false.

  ! Physics factors

  !> Physical scaling factor for length
  double precision :: unit_length=1.d0

  !> Physical scaling factor for time
  double precision :: unit_time=1.d0

  !> Physical scaling factor for density
  double precision :: unit_density=1.d0

  !> Physical scaling factor for velocity
  double precision :: unit_velocity=0.d0

  !> Physical scaling factor for temperature
  double precision :: unit_temperature=1.d0

  !> Physical scaling factor for pressure
  double precision :: unit_pressure=1.d0

  !> check and optionally fix unphysical small values (density, gas pressure)
  logical :: fix_small_values=.false.

  !> check unphysical small values
  logical :: check_NaN_values=.false.

  !> Whether to apply small value fixes to certain variables
  logical, allocatable :: small_values_fix_iw(:)

  !> Use SI units (.true.) or use cgs units (.false.)
  logical :: SI_unit=.false.

  !> Use multigrid (only available in 2D and 3D)
  logical :: use_multigrid = .false.

  ! AMR switches

  !> The maximum number of grid blocks in a processor
  integer :: max_blocks

  !> The maximum number of levels in the grid refinement
  integer, parameter :: nlevelshi = 20

  !> Maximal number of AMR levels
  integer :: refine_max_level

  !> Weights of variables used to calculate error for mesh refinement
  double precision, allocatable :: w_refine_weight(:)

  !> Fix the AMR grid after this time
  double precision :: tfixgrid

  !> Fix the AMR grid after this many time steps
  integer :: itfixgrid

  !> Reconstruct the AMR grid once every ditregrid iteration(s)
  integer :: ditregrid

  !> refinement: lohner estimate wavefilter setting
  double precision, allocatable :: amr_wavefilter(:)

  integer                       :: refine_criterion

  !> Error tolerance for refinement decision
  double precision, allocatable :: refine_threshold(:)
  double precision, allocatable :: derefine_ratio(:)

  !> If true, rebuild the AMR grid upon restarting
  logical :: reset_grid
  !> True for using stagger grid
  logical :: stagger_grid=.false.

  !> Number of cells as buffer zone
  !> \todo is it necessary?
  integer :: nbufferx1,nbufferx2

  integer :: levmin
  integer :: levmax
  integer :: levmax_sub

  ! Miscellaneous 

  !> problem switch allowing different setups in same usr_mod.t
  integer           :: iprob

  !> Kronecker delta tensor
  integer :: kr(3,3)

  !> Levi-Civita tensor
  integer :: lvc(3,3,3)
  !> 4-Levi-Civita tensor
  integer :: four_lvc(0:3,0:3,0:3,0:3)

  ! Time integration aspects 

  double precision :: dt
  double precision, allocatable :: dt_grid(:)

  logical :: time_advance

  !> The Courant (CFL) number used for the simulation
  double precision :: courantpar

  !> How to compute the CFL-limited time step.
  !>
  !> Options are 'maxsum': max(sum(c/dx)); 'summax': sum(max(c/dx)) and
  !> 'minimum: max(c/dx), where the summations loop over the grid dimensions and
  !> c is the velocity. The default 'maxsum' is the conventiontal way of
  !> computing CFL-limited time steps.
  character(len=std_len) :: typecourant

  !> If dtpar is positive, it sets the timestep dt, otherwise courantpar is used
  !> to limit the time step based on the Courant condition.
  double precision :: dtpar

  !> For resistive MHD, the time step is also limited by the diffusion time:
  !> \f$ dt < dtdiffpar \times dx^2/eta \f$
  double precision :: dtdiffpar

  !> The global simulation time
  double precision :: global_time

  !> Start time for the simulation
  double precision :: time_init

  !> End time for the simulation
  double precision :: time_max

  !> Ending wall time (in hours) for the simulation
  double precision :: wall_time_max

  !> Stop the simulation when the time step becomes smaller than this value
  double precision :: dtmin

  !> If true, reset iteration count and global_time to original values, and
  !> start writing snapshots at index 0
  logical :: reset_time

  !> If true, reset iteration count to 0
  logical :: reset_it

  !> If true, call initonegrid_usr upon restarting
  logical :: firstprocess

  !> If true, wall time is up, modify snapshotnext for later overwrite
  logical :: pass_wall_time

  !> Number of time steps taken
  integer :: it

  !> Stop the simulation after this many time steps have been taken
  integer :: it_max

  !> initial iteration count
  integer :: it_init

  !> If > 1, then in the first slowsteps-1 time steps dt is reduced
  !> by a factor \f$ 1 - (1- step/slowsteps)^2 \f$
  integer :: slowsteps

  ! Method switches

  !> Index of the sub-step in a multi-step time integrator
  integer :: istep

  !> How many sub-steps the time integrator takes
  integer :: nstep

  !> Which time stepper to use
  character(len=std_len) :: time_stepper

  !> Which time integrator to use
  character(len=std_len) :: time_integrator

  !> How to apply dimensional splitting to the source terms, see
  !> @ref disretization.md
  character(len=std_len) :: typesourcesplit

  !> Which spatial discretization to use (per grid level)
  character(len=std_len), allocatable :: flux_scheme(:)

  !> The spatial dicretization for the predictor step when using a two
  !> step method
  character(len=std_len), allocatable :: typepred1(:)

  !> Type of slope limiter used for reconstructing variables on cell edges
  integer, allocatable :: type_limiter(:)

  !> Type of slope limiter used for computing gradients or divergences, when
  !> typegrad or typediv are set to 'limited'
  integer, allocatable :: type_gradient_limiter(:)

  !> \todo Remove / replace with limiter
  integer :: typelimiter

  !> \todo Remove / replace with gradient_limiter
  integer :: typegradlimiter

  !> Limiter used for prolongation to refined grids and ghost cells
  character(len=std_len) :: typeprolonglimit

  !> Which type of entropy fix to use with Riemann-type solvers
  character(len=std_len), allocatable :: typeentropy(:)

  !> Which type of TVD method to use
  character(len=std_len) :: typetvd

  character(len=std_len) :: typeaverage
  character(len=std_len) :: typedimsplit
  character(len=std_len) :: geometry_name='default'
  character(len=std_len) :: typepoly

  integer                       :: nxdiffusehllc
  double precision, allocatable :: entropycoef(:)
  double precision              :: tvdlfeps
  logical, allocatable          :: loglimit(:), logflag(:)
  logical                       :: flathllc,flatcd,flatsh
  logical                       :: positivity_preserving
  !> Use split or unsplit way to add user's source terms, default: unsplit
  logical                       :: source_split_usr
  logical                       :: dimsplit

  !> RK2(alfa) method parameters from Butcher tableau
  double precision              :: rk2_alfa,rk_a21,rk_b1,rk_b2
  !> SSPRK choice of methods (both threestep and fourstep, Shu-Osher 2N* implementation)
  !> also fivestep SSPRK54
  integer                       :: ssprk_order
  double precision              :: rk_beta11,rk_beta22,rk_beta33,rk_beta44,&
     rk_c2,rk_c3,rk_c4
  double precision              :: rk_alfa21,rk_alfa22,rk_alfa31,rk_alfa33,&
     rk_alfa41,rk_alfa44
  double precision              :: rk_beta54,rk_beta55,rk_alfa53,rk_alfa54,&
     rk_alfa55,rk_c5
  !> RK3 Butcher table
  integer                       :: rk3_switch
  double precision              :: rk3_a21,rk3_a31,rk3_a32,rk3_b1,rk3_b2,&
     rk3_b3,rk3_c2,rk3_c3
  !> IMEX_ARS3 parameter ars_gamma
  double precision              :: ars_gamma
  !> IMEX_232 choice and parameters
  integer                       :: imex_switch
  double precision              :: imex_a21,imex_a31,imex_a32,imex_b1,imex_b2,&
     imex_ha21,imex_ha22
  double precision              :: imex_b3,imex_c2,imex_c3
  !> whether IMEX in use or not
  logical                       :: use_imex_scheme

  character(len=std_len) :: typediv,typegrad,typecurl

  !> need global maximal wave speed
  logical :: need_global_cmax=.false.

  ! Boundary region parameters

  !> True for dimensions with periodic boundaries
  logical :: periodB(ndim)

  !> Indicates whether there is a pole at a boundary
  logical :: poleB(2,ndim)

  !> True for dimensions with aperiodic boundaries
  logical :: aperiodB(ndim)

  !> True for save physical boundary cells in dat files
  logical :: save_physical_boundary

  !> True if a block has any physical boundary
  logical, allocatable :: phyboundblock(:)

  !> Array indicating the type of boundary condition per variable and per
  !> physical boundary
  character(len=std_len), allocatable :: typeboundary(:, :)

  character(len=std_len) :: typeghostfill='linear',prolongation_method
  logical :: internalboundary

  ! parameters for bc_phys
  integer, parameter :: ismin1=-1+2*1,ismin2=-1+2*2
  integer, parameter :: ismax1=2*1,ismax2=2*2

  !> Base file name for synthetic EUV emission output
  character(len=std_len) :: filename_euv
  !> output image
  logical :: image=.false.
  !> output spectrum
  logical :: spectrum=.false.
  ! wavelength for output
  integer :: wavelength
  !> direction of light of sight
  integer :: direction_LOS=3
  !> direction of slit (for spectrum)
  integer :: direction_slit=2
  !> location of the slit
  double precision :: location_slit=0.d0
  !> resolution of the output
  character(len=std_len) :: resolution_euv='instrument'

  !$OMP THREADPRIVATE(dxlevel)
  !$OMP THREADPRIVATE(saveigrid)
  !$OMP THREADPRIVATE(typelimiter,typegradlimiter)

contains

  !> Cross product of two vectors
  pure subroutine cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,a,b,axb)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,3)
    double precision, intent(out) :: axb(ixImin1:ixImax1,ixImin2:ixImax2,3)
    !-------------------------------------------------------------------------

    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)-a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)
    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)-a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
  end subroutine cross_product

  !> SNG function
  integer function sgn(x)
    integer, intent(in) :: x
    if (x==0) then
      sgn = 0
    else
      sgn = sign(1,x)
    end if
  end function sgn

end module mod_global_parameters
