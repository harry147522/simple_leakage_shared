module mod_physicaldata
   implicit none
   save

   type state
      !> ID of a grid block
      integer :: igrid=-1
      !> index range of block array in cell centers
      integer :: ixG^L
      !> index range of block array in cell faces
      integer :: ixGs^L
      !> location of w0-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iw0=0
      !> If it face a physical boundary
      logical, dimension(:), pointer :: is_physical_boundary(:) =>Null()
      !> Variables, cell center primitive values
      double precision, dimension(:^D&,:), allocatable :: prim
      !> Variables, cell center conservative values
      double precision, dimension(:^D&,:), allocatable :: cons
      !> Variables, cell face values
      double precision, dimension(:^D&,:), allocatable :: prims
      !> Variables, cell edge values
      double precision, dimension(:^D&,:), allocatable :: prime
      !> Variables, cell corner values
      double precision, dimension(:^D&,:), allocatable :: primc
      !> Cell-center positions
      double precision, dimension(:^D&,:), pointer :: x=>Null()
      !> Cell sizes in coordinate units
      double precision, dimension(:^D&,:), pointer :: dx=>Null()
      !> Cell sizes at cell center in length unit
      double precision, dimension(:^D&,:), pointer :: ds=>Null()
      !> Cell sizes at cell face in length unit
      double precision, dimension(:^D&,:), pointer :: dsC=>Null()
      !> Volumes of a cell
      double precision, dimension(:^D&), pointer :: dvolume=>Null()
      !> Areas of cell-center surfaces
      double precision, dimension(:^D&,:), pointer :: surface=>Null()
      !> Areas of cell-face surfaces
      double precision, dimension(:^D&,:), pointer :: surfaceC=>Null()
      !> Christoffel symbols of a cell 
      double precision, dimension(:^D&,:,:,:), pointer :: christoffel=>Null()
   end type state

{^NOONED
   type state_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> Variables, normally center
      double precision, dimension(:^DE&,:), allocatable :: prim
      !> Variables for the cornerpositions on the slice 
      double precision, dimension(:^DE&,:), allocatable :: primC
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:^DE&,:), allocatable :: primcoarse
      !> Cell-center positions
      double precision, dimension(:^DE&,:), allocatable :: x
      !> Corner positions on the slice
      double precision, dimension(:^DE&,:), allocatable :: xC
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:^DE&,:), allocatable :: xcoarse
      !> Cell sizes
      double precision, dimension(:^DE&,:), allocatable :: dx
      !> Cell sizes, one level coarser
      double precision, dimension(:^D&,:), allocatable :: dxcoarse
      !> Cell sizes in length unit
      double precision, dimension(:^D&,:), allocatable :: ds
      !> Volumes of a cell
      double precision, dimension(:^DE&), allocatable :: dvolume
      !> Volumes of a cell, one level coarser representative
      double precision, dimension(:^DE&), allocatable :: dvolumecoarse
      !> Areas of cell-center surfaces 
      double precision, dimension(:^DE&,:), allocatable :: surface
      !> Areas of cell-face surfaces
      double precision, dimension(:^DE&,:), allocatable :: surfaceC
   end type state_sub
}
{^IFONED
   type state_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> Variables, normally center
      double precision, dimension(:), allocatable :: prim
      !> Variables for the cornerpositions on the slice 
      double precision, dimension(:), allocatable :: primC
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:), allocatable :: primcoarse
      !> Cell-center positions
      double precision, dimension(:), allocatable :: x
      !> Corner positions on the slice
      double precision, dimension(:), allocatable :: xC
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:), allocatable :: xcoarse
   end type state_sub
}
   type grid_field
      !> Variables new state
      double precision, dimension(:^D&,:), allocatable :: prim
      !> Variables old state
      double precision, dimension(:^D&,:), allocatable :: primold
   end type grid_field

   !> Block pointer for using one block and its previous state
   type(state), pointer :: block, block0
   !> buffer for pole boundary
   type(state) :: pole_buf

   !> array of physical states for all blocks on my processor
   type(state), dimension(:), allocatable, target :: ps
   !> array of physical states, temp 1 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps1
   !> array of physical states, temp 2 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps2
   !> array of physical states, temp 3 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps3
   !> array of physical states, temp 4 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps4
   !> array of physical states, at the beginning of each iteration
   type(state), dimension(:), allocatable, target :: pso
   !> array of physical blocks, one level coarser representative
   type(state), dimension(:), allocatable, target :: psc

   !> array of physical blocks in reduced dimension
   type(state_sub), dimension(:), allocatable, target :: ps_sub

{^IFONED
   double precision, dimension(:), allocatable :: collapsedData
}
{^NOONED
   double precision, dimension(:^DE&,:), allocatable :: collapsedData
}
end module mod_physicaldata
