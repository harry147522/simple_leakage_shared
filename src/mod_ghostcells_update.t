!> update ghost cells of all blocks including physical boundaries 
module mod_ghostcells_update

  implicit none
  ! Special buffer for pole copy
  type wbuffer
    double precision, dimension(:^D&,:), allocatable :: prim
  end type wbuffer

  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.
  integer :: ixM^L, ixCoG^L, ixCoM^L

  !> The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  ! The first index goes from -1:2, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions, 2 when a block touched both lower and upper
  ! boundary

  ! index ranges to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(-1:2,-1:1) :: ixS_srl_^L, ixR_srl_^L

  ! index ranges of staggered variables to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(^ND,-1:1) :: ixS_srl_stg_^L, ixR_srl_stg_^L

  ! index ranges to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(-1:1,-1:1) :: ixS_r_^L

  ! index ranges of staggered variables to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(^ND,-1:1) :: ixS_r_stg_^L

  ! index ranges to receive restriced ghost cells from finer blocks 
  integer, dimension(-1:1, 0:3) :: ixR_r_^L

  ! index ranges of staggered variables to receive restriced ghost cells from finer blocks 
  integer, dimension(^ND,0:3)  :: ixR_r_stg_^L

  ! send prolongated (p) ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixS_p_^L, ixR_p_^L

  ! send prolongated (p) staggered ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(^ND,0:3)  :: ixS_p_stg_^L, ixR_p_stg_^L

  ! number of MPI receive-send pairs, srl: same refinement level; r: restrict; p: prolong
  integer :: nrecv_bc_srl, nsend_bc_srl, nrecv_bc_r, nsend_bc_r, nrecv_bc_p, nsend_bc_p

  ! total size of buffer arrays
  integer :: nbuff_bc_recv_srl, nbuff_bc_send_srl, nbuff_bc_recv_r, nbuff_bc_send_r, nbuff_bc_recv_p, nbuff_bc_send_p

  ! record index position of buffer arrays
  integer :: ibuf_send_srl, ibuf_recv_srl, ibuf_send_r, ibuf_recv_r, ibuf_send_p, ibuf_recv_p

  ! count of times of send and receive
  integer :: isend_srl, irecv_srl, isend_r, irecv_r, isend_p, irecv_p

  ! count of times of send and receive for cell center ghost cells
  integer :: isend_c, irecv_c

  ! tag of MPI send and recv
  integer, private :: itag

  ! total sizes = cell-center normal flux + stagger-grid flux of send and receive
  integer, dimension(-1:1^D&) :: sizes_srl_send_total, sizes_srl_recv_total

  ! sizes of buffer arrays for center-grid variable for siblings and restrict
  integer, dimension(:), allocatable :: recvrequest_c_sr, sendrequest_c_sr
  integer, dimension(:,:), allocatable :: recvstatus_c_sr, sendstatus_c_sr

  ! sizes of buffer arrays for center-grid variable for prolongation
  integer, dimension(:), allocatable :: recvrequest_c_p, sendrequest_c_p
  integer, dimension(:,:), allocatable :: recvstatus_c_p, sendstatus_c_p

  ! sizes of buffer arrays for stagger-grid variable
  integer, dimension(^ND,-1:1^D&) :: sizes_srl_send_stg, sizes_srl_recv_stg

  integer, dimension(:), allocatable :: recvrequest_srl, sendrequest_srl
  integer, dimension(:,:), allocatable :: recvstatus_srl, sendstatus_srl
 
  ! buffer arrays for send and receive of siblings, allocate in build_connectivity
  double precision, dimension(:), allocatable :: recvbuffer_srl, sendbuffer_srl

  integer, dimension(:), allocatable :: recvrequest_r, sendrequest_r
  integer, dimension(:,:), allocatable :: recvstatus_r, sendstatus_r

  ! buffer arrays for send and receive in restriction
  double precision, dimension(:), allocatable :: recvbuffer_r, sendbuffer_r

  integer, dimension(:), allocatable :: recvrequest_p, sendrequest_p
  integer, dimension(:,:), allocatable :: recvstatus_p, sendstatus_p

  ! buffer arrays for send and receive in prolongation
  double precision, dimension(:), allocatable :: recvbuffer_p, sendbuffer_p

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(-1:1^D&)     :: sizes_r_send_total
  integer, dimension(0:3^D&)      :: sizes_r_recv_total
  integer, dimension(^ND,-1:1^D&) :: sizes_r_send_stg
  integer, dimension(^ND,0:3^D&)  :: sizes_r_recv_stg

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(0:3^D&)      :: sizes_p_send_total, sizes_p_recv_total
  integer, dimension(^ND,0:3^D&)  :: sizes_p_send_stg, sizes_p_recv_stg

  ! There are two variants, _f indicates that all flux variables are filled,
  ! whereas _p means that part of the variables is filled 
  ! Furthermore _r_ stands for restrict, _p_ for prolongation.
  integer, dimension(-1:2^D&,-1:1^D&), target :: type_send_srl_f, type_recv_srl_f
  integer, dimension(-1:1^D&,-1:1^D&), target :: type_send_r_f
  integer, dimension(-1:1^D&, 0:3^D&), target :: type_recv_r_f, type_send_p_f, type_recv_p_f
  integer, dimension(-1:2^D&,-1:1^D&), target :: type_send_srl_p1, type_recv_srl_p1
  integer, dimension(-1:1^D&,-1:1^D&), target :: type_send_r_p1
  integer, dimension(-1:1^D&, 0:3^D&), target :: type_recv_r_p1, type_send_p_p1, type_recv_p_p1
  integer, dimension(-1:2^D&,-1:1^D&), target :: type_send_srl_p2, type_recv_srl_p2
  integer, dimension(-1:1^D&,-1:1^D&), target :: type_send_r_p2
  integer, dimension(-1:1^D&, 0:3^D&), target :: type_recv_r_p2, type_send_p_p2, type_recv_p_p2
  integer, dimension(:^D&,:^D&), pointer :: type_send_srl, type_recv_srl, type_send_r
  integer, dimension(:^D&,:^D&), pointer :: type_recv_r, type_send_p, type_recv_p

contains

  subroutine init_bc()
    use mod_global_parameters 
    use mod_physics, only: phys_req_diagonal

    integer :: nghostcellsCo, interpolation_order
    integer :: nx^D, nxCo^D, ixG^L, i^D, ic^D, inc^D, idir

    ixG^L=ixG^LL;
    ixM^L=ixG^L^LSUBnghostcells;
    ixCoGmin^D=1;
    !ixCoGmax^D=ixGmax^D/2+nghostcells;
    ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;

    ixCoM^L=ixCoG^L^LSUBnghostcells;

    nx^D=ixMmax^D-ixMmin^D+1;
    nxCo^D=nx^D/2;
    
    select case (typeghostfill)
    case ("copy")
       interpolation_order=1
    case ("linear")
       interpolation_order=2
    case default
       write (unitterm,*) "Undefined typeghostfill ",typeghostfill
       call mpistop("Undefined typeghostfill")
    end select
    nghostcellsCo=int((nghostcells+1)/2)
    
    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc too high")
    end if

    ! (iib,i) index has following meanings: iib = 0 means it is not at any physical boundary
    ! iib=-1 means it is at the minimum side of a physical boundary  
    ! iib= 1 means it is at the maximum side of a physical boundary  
    ! i=-1 means subregion prepared for the neighbor at its minimum side 
    ! i= 1 means subregion prepared for the neighbor at its maximum side 
    {
    ixS_srl_min^D(:,-1)=ixMmin^D
    ixS_srl_min^D(:, 0)=ixMmin^D
    ixS_srl_min^D(:, 1)=ixMmax^D+1-nghostcells
    ixS_srl_max^D(:,-1)=ixMmin^D-1+nghostcells
    ixS_srl_max^D(:, 0)=ixMmax^D
    ixS_srl_max^D(:, 1)=ixMmax^D
     
    ixR_srl_min^D(:,-1)=1
    ixR_srl_min^D(:, 0)=ixMmin^D
    ixR_srl_min^D(:, 1)=ixMmax^D+1
    ixR_srl_max^D(:,-1)=nghostcells
    ixR_srl_max^D(:, 0)=ixMmax^D
    ixR_srl_max^D(:, 1)=ixGmax^D
    
    ixS_r_min^D(:,-1)=ixCoMmin^D
    ixS_r_min^D(:, 0)=ixCoMmin^D
    ixS_r_min^D(:, 1)=ixCoMmax^D+1-nghostcells
    ixS_r_max^D(:,-1)=ixCoMmin^D-1+nghostcells
    ixS_r_max^D(:, 0)=ixCoMmax^D
    ixS_r_max^D(:, 1)=ixCoMmax^D
    
    ixR_r_min^D(:, 0)=1
    ixR_r_min^D(:, 1)=ixMmin^D
    ixR_r_min^D(:, 2)=ixMmin^D+nxCo^D
    ixR_r_min^D(:, 3)=ixMmax^D+1
    ixR_r_max^D(:, 0)=nghostcells
    ixR_r_max^D(:, 1)=ixMmin^D-1+nxCo^D
    ixR_r_max^D(:, 2)=ixMmax^D
    ixR_r_max^D(:, 3)=ixGmax^D

    ixS_p_min^D(:, 0)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(:, 1)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(:, 2)=ixMmin^D+nxCo^D-nghostcellsCo-(interpolation_order-1)
    ixS_p_min^D(:, 3)=ixMmax^D+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max^D(:, 0)=ixMmin^D-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(:, 1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(:, 2)=ixMmax^D+(interpolation_order-1)
    ixS_p_max^D(:, 3)=ixMmax^D+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
      ! exclude ghost-cell region when diagonal cells are unknown
      ixS_p_min^D(:, 0)=ixMmin^D
      ixS_p_max^D(:, 3)=ixMmax^D
      ixS_p_max^D(:, 1)=ixMmin^D-1+nxCo^D+(interpolation_order-1)
      ixS_p_min^D(:, 2)=ixMmin^D+nxCo^D-(interpolation_order-1)
    end if

    ixR_p_min^D(:, 0)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(:, 1)=ixCoMmin^D-(interpolation_order-1)
    ixR_p_min^D(:, 2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(:, 3)=ixCoMmax^D+1-(interpolation_order-1)
    ixR_p_max^D(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max^D(:, 1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
    ixR_p_max^D(:, 2)=ixCoMmax^D+(interpolation_order-1)
    ixR_p_max^D(:, 3)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
      ixR_p_max^D(:, 0)=nghostcells
      ixR_p_min^D(:, 3)=ixCoMmax^D+1
      ixR_p_max^D(:, 1)=ixCoMmax^D+(interpolation_order-1)
      ixR_p_min^D(:, 2)=ixCoMmin^D-(interpolation_order-1)
    end if

    \}

    if (stagger_grid) then
      allocate(pole_buf%prims(ixGs^T,nprims))
      ! Staggered (face-allocated) variables
      do idir=1,ndim
      { ixS_srl_stg_min^D(idir,-1)=ixMmin^D-kr(idir,^D)
        ixS_srl_stg_max^D(idir,-1)=ixMmin^D-1+nghostcells
        ixS_srl_stg_min^D(idir,0) =ixMmin^D-kr(idir,^D)
        ixS_srl_stg_max^D(idir,0) =ixMmax^D
        ixS_srl_stg_min^D(idir,1) =ixMmax^D-nghostcells+1-kr(idir,^D)
        ixS_srl_stg_max^D(idir,1) =ixMmax^D
        
        ixR_srl_stg_min^D(idir,-1)=1-kr(idir,^D)
        ixR_srl_stg_max^D(idir,-1)=nghostcells
        ixR_srl_stg_min^D(idir,0) =ixMmin^D-kr(idir,^D)
        ixR_srl_stg_max^D(idir,0) =ixMmax^D
        ixR_srl_stg_min^D(idir,1) =ixMmax^D+1-kr(idir,^D)
        ixR_srl_stg_max^D(idir,1) =ixGmax^D

        ixS_r_stg_min^D(idir,-1)=ixCoMmin^D-kr(idir,^D)
        ixS_r_stg_max^D(idir,-1)=ixCoMmin^D-1+nghostcells
        ixS_r_stg_min^D(idir,0) =ixCoMmin^D-kr(idir,^D)
        ixS_r_stg_max^D(idir,0) =ixCoMmax^D
        ixS_r_stg_min^D(idir,1) =ixCoMmax^D+1-nghostcells-kr(idir,^D)
        ixS_r_stg_max^D(idir,1) =ixCoMmax^D
 
        ixR_r_stg_min^D(idir,0)=1-kr(idir,^D)
        ixR_r_stg_max^D(idir,0)=nghostcells
        ixR_r_stg_min^D(idir,1)=ixMmin^D-kr(idir,^D)
        ixR_r_stg_max^D(idir,1)=ixMmin^D-1+nxCo^D
        ixR_r_stg_min^D(idir,2)=ixMmin^D+nxCo^D-kr(idir,^D)
        ixR_r_stg_max^D(idir,2)=ixMmax^D
        ixR_r_stg_min^D(idir,3)=ixMmax^D+1-kr(idir,^D)
        ixR_r_stg_max^D(idir,3)=ixGmax^D
        \}
       {if (idir==^D) then
          ! Parallel components
          {
          ixS_p_stg_min^D(idir,0)=ixMmin^D-1 ! -1 to make redundant 
          ixS_p_stg_max^D(idir,0)=ixMmin^D-1+nghostcellsCo
          ixS_p_stg_min^D(idir,1)=ixMmin^D-1 ! -1 to make redundant 
          ixS_p_stg_max^D(idir,1)=ixMmin^D-1+nxCo^D+nghostcellsCo
          ixS_p_stg_min^D(idir,2)=ixMmax^D-nxCo^D-nghostcellsCo
          ixS_p_stg_max^D(idir,2)=ixMmax^D
          ixS_p_stg_min^D(idir,3)=ixMmax^D-nghostcellsCo
          ixS_p_stg_max^D(idir,3)=ixMmax^D

          ixR_p_stg_min^D(idir,0)=ixCoMmin^D-1-nghostcellsCo
          ixR_p_stg_max^D(idir,0)=ixCoMmin^D-1
          ixR_p_stg_min^D(idir,1)=ixCoMmin^D-1 ! -1 to make redundant 
          ixR_p_stg_max^D(idir,1)=ixCoMmax^D+nghostcellsCo
          ixR_p_stg_min^D(idir,2)=ixCoMmin^D-1-nghostcellsCo
          ixR_p_stg_max^D(idir,2)=ixCoMmax^D
          ixR_p_stg_min^D(idir,3)=ixCoMmax^D+1-1 ! -1 to make redundant 
          ixR_p_stg_max^D(idir,3)=ixCoMmax^D+nghostcellsCo
          \}
        else
          {
          ! Perpendicular component
          ixS_p_stg_min^D(idir,0)=ixMmin^D
          ixS_p_stg_max^D(idir,0)=ixMmin^D-1+nghostcellsCo+(interpolation_order-1)
          ixS_p_stg_min^D(idir,1)=ixMmin^D
          ixS_p_stg_max^D(idir,1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
          ixS_p_stg_min^D(idir,2)=ixMmax^D+1-nxCo^D-nghostcellsCo-(interpolation_order-1)
          ixS_p_stg_max^D(idir,2)=ixMmax^D
          ixS_p_stg_min^D(idir,3)=ixMmax^D+1-nghostcellsCo-(interpolation_order-1)
          ixS_p_stg_max^D(idir,3)=ixMmax^D
 
          ixR_p_stg_min^D(idir,0)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
          ixR_p_stg_max^D(idir,0)=ixCoMmin^D-1
          ixR_p_stg_min^D(idir,1)=ixCoMmin^D
          ixR_p_stg_max^D(idir,1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
          ixR_p_stg_min^D(idir,2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
          ixR_p_stg_max^D(idir,2)=ixCoMmax^D
          ixR_p_stg_min^D(idir,3)=ixCoMmax^D+1
          ixR_p_stg_max^D(idir,3)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
          \}
        end if
        }
      end do
      ! calculate sizes for buffer arrays for siblings
      {do i^DB=-1,1\}
         ! Staggered (face-allocated) variables
         do idir=1,ndim
           sizes_srl_send_stg(idir,i^D)={(ixS_srl_stg_max^D(idir,i^D)-ixS_srl_stg_min^D(idir,i^D)+1)|*}
           sizes_srl_recv_stg(idir,i^D)={(ixR_srl_stg_max^D(idir,i^D)-ixR_srl_stg_min^D(idir,i^D)+1)|*}
           sizes_r_send_stg(idir,i^D)={(ixS_r_stg_max^D(idir,i^D)-ixS_r_stg_min^D(idir,i^D)+1)|*}
         end do
         sizes_srl_send_total(i^D)=sum(sizes_srl_send_stg(:,i^D))
         sizes_srl_recv_total(i^D)=sum(sizes_srl_recv_stg(:,i^D))
         sizes_r_send_total(i^D)=sum(sizes_r_send_stg(:,i^D))
      {end do\}

      {do i^DB=0,3\}
         ! Staggered (face-allocated) variables
           do idir=1,ndim
             sizes_r_recv_stg(idir,i^D)={(ixR_r_stg_max^D(idir,i^D)-ixR_r_stg_min^D(idir,i^D)+1)|*}
             sizes_p_send_stg(idir,i^D)={(ixS_p_stg_max^D(idir,i^D)-ixS_p_stg_min^D(idir,i^D)+1)|*}
             sizes_p_recv_stg(idir,i^D)={(ixR_p_stg_max^D(idir,i^D)-ixR_p_stg_min^D(idir,i^D)+1)|*}
           end do
           sizes_r_recv_total(i^D)=sum(sizes_r_recv_stg(:,i^D))
           sizes_p_send_total(i^D)=sum(sizes_p_send_stg(:,i^D))
           sizes_p_recv_total(i^D)=sum(sizes_p_recv_stg(:,i^D))
      {end do\}
    else
      ! extend index range to physical boundary
      {
      ixS_srl_min^D(-1,0)=1
      ixS_srl_min^D( 1,0)=ixMmin^D
      ixS_srl_min^D( 2,0)=1
      ixS_srl_max^D(-1,0)=ixMmax^D
      ixS_srl_max^D( 1,0)=ixGmax^D
      ixS_srl_max^D( 2,0)=ixGmax^D
       
      ixR_srl_min^D(-1,0)=1
      ixR_srl_min^D( 1,0)=ixMmin^D
      ixR_srl_min^D( 2,0)=1
      ixR_srl_max^D(-1,0)=ixMmax^D
      ixR_srl_max^D( 1,0)=ixGmax^D
      ixR_srl_max^D( 2,0)=ixGmax^D
      
      ixS_r_min^D(-1,0)=1
      ixS_r_min^D( 1,0)=ixCoMmin^D
      ixS_r_max^D(-1,0)=ixCoMmax^D
      ixS_r_max^D( 1,0)=ixCoGmax^D
      
      ixR_r_min^D(-1,1)=1
      ixR_r_max^D(-1,1)=ixMmin^D-1+nxCo^D
      ixR_r_min^D( 1,2)=ixMmin^D+nxCo^D
      ixR_r_max^D( 1,2)=ixGmax^D

      ixS_p_min^D(-1,1)=1
      ixS_p_max^D( 1,2)=ixGmax^D

      ixR_p_min^D(-1,1)=1
      ixR_p_max^D( 1,2)=ixCoGmax^D
      \}
    end if

  end subroutine init_bc

  subroutine create_bc_mpi_datatype(nwstart,nwbc) 
    use mod_global_parameters 

    integer, intent(in) :: nwstart, nwbc
    integer :: i^D, ic^D, inc^D, iib^D

    {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle
      {do iib^DB=-1,2\}
         call get_bc_comm_type(type_send_srl(iib^D,i^D),ixS_srl_^L(iib^D,i^D),ixG^LL,nwstart,nwbc)
         call get_bc_comm_type(type_recv_srl(iib^D,i^D),ixR_srl_^L(iib^D,i^D),ixG^LL,nwstart,nwbc)
         if (iib^D==2|.or.) cycle
         call get_bc_comm_type(type_send_r(iib^D,i^D),  ixS_r_^L(iib^D,i^D),ixCoG^L,nwstart,nwbc)
         {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
            inc^DB=2*i^DB+ic^DB\}
            call get_bc_comm_type(type_recv_r(iib^D,inc^D),ixR_r_^L(iib^D,inc^D), ixG^LL,nwstart,nwbc)
            call get_bc_comm_type(type_send_p(iib^D,inc^D),ixS_p_^L(iib^D,inc^D), ixG^LL,nwstart,nwbc)
            call get_bc_comm_type(type_recv_p(iib^D,inc^D),ixR_p_^L(iib^D,inc^D),ixCoG^L,nwstart,nwbc)
         {end do\}
      {end do\}
    {end do\}
  
  end subroutine create_bc_mpi_datatype

  subroutine get_bc_comm_type(comm_type,ix^L,ixG^L,nwstart,nwbc)
    use mod_global_parameters 
  
    integer, intent(inout) :: comm_type
    integer, intent(in) :: ix^L, ixG^L, nwstart, nwbc
    
    integer, dimension(ndim+1) :: fullsize, subsize, start

    ^D&fullsize(^D)=ixGmax^D;
    fullsize(ndim+1)=nprim
    ^D&subsize(^D)=ixmax^D-ixmin^D+1;
    subsize(ndim+1)=nwbc
    ^D&start(^D)=ixmin^D-1;
    start(ndim+1)=nwstart-1
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,fullsize,subsize,start,MPI_ORDER_FORTRAN, &
                                  MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
    call MPI_TYPE_COMMIT(comm_type,ierrmpi)
    
  end subroutine get_bc_comm_type

  subroutine put_bc_comm_types()
    use mod_global_parameters 
 
    integer :: i^D, ic^D, inc^D, iib^D

    {do i^DB=-1,1\}
       if (i^D==0|.and.) cycle
       {do iib^DB=-1,2\}
           call MPI_TYPE_FREE(type_send_srl(iib^D,i^D),ierrmpi)
           call MPI_TYPE_FREE(type_recv_srl(iib^D,i^D),ierrmpi)
           if (levmin==levmax) cycle
           if (iib^D==2|.or.) cycle
           call MPI_TYPE_FREE(type_send_r(iib^D,i^D),ierrmpi)
           {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
              inc^DB=2*i^DB+ic^DB\}
              call MPI_TYPE_FREE(type_recv_r(iib^D,inc^D),ierrmpi)
              call MPI_TYPE_FREE(type_send_p(iib^D,inc^D),ierrmpi)
              call MPI_TYPE_FREE(type_recv_p(iib^D,inc^D),ierrmpi)
           {end do\}
       {end do\}
    {end do\}
  
  end subroutine put_bc_comm_types

  !> do update ghost cells of all blocks including physical boundaries
  subroutine getbc(time,qdt,psb,nwstart,nwbc,req_diag)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in)      :: time, qdt
    type(state), target               :: psb(max_blocks)
    integer, intent(in)               :: nwstart ! Fill from nwstart
    integer, intent(in)               :: nwbc    ! Number of variables to fill
    logical, intent(in), optional     :: req_diag ! If false, skip diagonal ghost cells

    double precision :: time_bcin
    integer :: my_neighbor_type, ipole, idims, iside, nwhead, nwtail
    integer :: iigrid, igrid, ineighbor, ipe_neighbor
    integer :: nrecvs, nsends, isizes
    integer :: ixG^L, ixR^L, ixS^L, ixB^L, ixI^L, k^L
    integer :: i^D, n_i^D, ic^D, inc^D, n_inc^D, iib^D, idir
    ! store physical boundary indicating index
    integer :: idphyb(ndim,max_blocks),bindex(ndim)
    integer :: isend_buf(npwbuf), ipwbuf, nghostcellsco,iB
    ! index pointer for buffer arrays as a start for a segment
    integer :: ibuf_start, ibuf_next
    ! shapes of reshape
    integer, dimension(1) :: shapes
    logical  :: req_diagonal, NeedProlong(-1:1^D&)
    type(wbuffer) :: pwbuf(npwbuf)

    ! Stretching grid parameters for coarsened block of the current block

    nwhead=nwstart
    nwtail=nwstart+nwbc-1

    req_diagonal = .true.
    if (present(req_diag)) req_diagonal = req_diag

    time_bcin=MPI_WTIME()
    ixG^L=ixG^LL;
    
    if (internalboundary) then 
       call getintbc(time,ixG^L)
    end if
    ! fill ghost cells in physical boundaries
    if(bcphys.and. .not.stagger_grid) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>psb(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idims=1,ndim
          ! to avoid using as yet unknown corner info in more than 1D, we
          ! fill only interior mesh ranges of the ghost cell ranges at first,
          ! and progressively enlarge the ranges to include corners later
          {
           kmin^D=merge(0, 1, idims==^D)
           kmax^D=merge(0, 1, idims==^D)
           ixBmin^D=ixGmin^D+kmin^D*nghostcells
           ixBmax^D=ixGmax^D-kmax^D*nghostcells
          \}
          {^IFTWOD
          if(idims > 1 .and. neighbor_type(-1,0,igrid)==neighbor_boundary) ixBmin1=ixGmin1
          if(idims > 1 .and. neighbor_type( 1,0,igrid)==neighbor_boundary) ixBmax1=ixGmax1}
          {^IFTHREED
          if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==neighbor_boundary) ixBmin1=ixGmin1
          if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==neighbor_boundary) ixBmax1=ixGmax1
          if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==neighbor_boundary) ixBmin2=ixGmin2
          if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==neighbor_boundary) ixBmax2=ixGmax2}
          do iside=1,2
            i^D=kr(^D,idims)*(2*iside-3);
            if (aperiodB(idims)) then
              if (neighbor_type(i^D,igrid) /= neighbor_boundary .and. &
                   .not. psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
            else
              if (neighbor_type(i^D,igrid) /= neighbor_boundary) cycle
            end if
            call bc_phys(iside,idims,time,qdt,psb(igrid),ixG^L,ixB^L)
          end do
        end do
      end do
    end if

    ! default : no singular axis
    ipole=0
    
    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    if(stagger_grid) then
      ibuf_recv_srl=1
      ibuf_recv_r=1
      ibuf_recv_p=1
      ibuf_send_srl=1
      ibuf_send_r=1
      ibuf_send_p=1
      irecv_srl=0
      irecv_r=0
      irecv_p=0
      isend_srl=0
      isend_r=0
      isend_p=0
    end if
    ! receiving ghost-cell values from sibling blocks and finer neighbors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       call identifyphysbound(ps(igrid),iib^D)   
       ^D&idphyb(^D,igrid)=iib^D;
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          my_neighbor_type=neighbor_type(i^D,igrid)
          select case (my_neighbor_type)
          case (neighbor_sibling)
             call bc_recv_srl
          case (neighbor_fine)
             call bc_recv_restrict
          end select
       {end do\}
    end do
    
    ! sending ghost-cell values to sibling blocks and coarser neighbors
    nghostcellsco=ceiling(nghostcells*0.5d0)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>psb(igrid)

       ! Used stored data to identify physical boundaries
       ^D&iib^D=idphyb(^D,igrid);

       if (any(neighbor_type(:^D&,igrid)==neighbor_coarse)) then
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    {#IFDEF EVOLVINGBOUNDARY
          if(phyboundblock(igrid)) then
            ! coarsen finer ghost cells at physical boundaries
            ixCoMmin^D=ixCoGmin^D+nghostcellsco;
            ixCoMmax^D=ixCoGmax^D-nghostcellsco;
            ixMmin^D=ixGmin^D+(nghostcellsco-1);
            ixMmax^D=ixGmax^D-(nghostcellsco-1);
          else
            ixCoM^L=ixCoG^L^LSUBnghostcells;
            ixM^L=ixG^L^LSUBnghostcells;
          end if
    }
          call coarsen_grid(psb(igrid),ixG^L,ixM^L,psc(igrid),ixCoG^L,ixCoM^L)
       end if
    
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if (phi_ > 0) ipole=neighbor_pole(i^D,igrid)
          my_neighbor_type=neighbor_type(i^D,igrid)
          select case (my_neighbor_type)
          case (neighbor_sibling)
             call bc_send_srl
          case (neighbor_coarse)
             call bc_send_restrict
          end select
       {end do\}
    end do
    
    call MPI_WAITALL(irecv_c,recvrequest_c_sr,recvstatus_c_sr,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_sr,sendstatus_c_sr,ierrmpi)

    if(stagger_grid) then
      call MPI_WAITALL(nrecv_bc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
      call MPI_WAITALL(nsend_bc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
      call MPI_WAITALL(nrecv_bc_r,recvrequest_r,recvstatus_r,ierrmpi)
      call MPI_WAITALL(nsend_bc_r,sendrequest_r,sendstatus_r,ierrmpi)
      ! unpack the received data to fill ghost cells
      ibuf_recv_srl=1
      ibuf_recv_r=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          select case (neighbor_type(i^D,igrid))
          case (neighbor_sibling)
             call bc_fill_srl
          case (neighbor_fine)
             call bc_fill_r
          end select
       {end do\}
      end do
    end if

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%prim)
    end do

    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    ! receiving ghost-cell values from coarser neighbors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          my_neighbor_type=neighbor_type(i^D,igrid)
          if (my_neighbor_type==neighbor_coarse) call bc_recv_prolong
       {end do\}
    end do
    ! sending ghost-cell values to finer neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>psb(igrid)
       ^D&iib^D=idphyb(^D,igrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       if (any(neighbor_type(:^D&,igrid)==neighbor_fine)) then
          {do i^DB=-1,1\}
             if (skip_direction([ i^D ])) cycle
             if (phi_ > 0) ipole=neighbor_pole(i^D,igrid)
             my_neighbor_type=neighbor_type(i^D,igrid)
             if (my_neighbor_type==neighbor_fine) call bc_send_prolong
          {end do\}
       end if
    end do

    call MPI_WAITALL(irecv_c,recvrequest_c_p,recvstatus_c_p,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_p,sendstatus_c_p,ierrmpi)

    if(stagger_grid) then
      call MPI_WAITALL(nrecv_bc_p,recvrequest_p,recvstatus_p,ierrmpi)
      call MPI_WAITALL(nsend_bc_p,sendrequest_p,sendstatus_p,ierrmpi)

      ! fill coarser representative after receipt
      ibuf_recv_p=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         ^D&iib^D=idphyb(^D,igrid);
         {do i^DB=-1,1\}
            if (skip_direction([ i^D ])) cycle
            if(neighbor_type(i^D,igrid)==neighbor_coarse) call bc_fill_p
         {end do\}
      end do
    end if
    ! do prolongation on the ghost-cell values received from coarser neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&iib^D=idphyb(^D,igrid);
       if (any(neighbor_type(:^D&,igrid)==neighbor_coarse)) then
         NeedProlong=.false.
         {do i^DB=-1,1\}
            if (skip_direction([ i^D ])) cycle
            my_neighbor_type=neighbor_type(i^D,igrid)
            if (my_neighbor_type==neighbor_coarse) then
              call bc_prolong
              NeedProlong(i^D)=.true.
            end if
         {end do\}
         if(stagger_grid) then
           ! Ghost cell prolongation for staggered variables
           ! must be done in a specific order.
           ! First the first neighbours, which have 2 indices=0 in 3D
           ! or one index=0 in 2D
           do idims=1,ndim
             i^D=0;
             select case(idims)
            {case(^D)
               do i^D=-1,1,2
                 if (NeedProlong(i^DD)) call bc_prolong_stg(NeedProlong)
               end do
             \}
             end select
           end do
           ! Then the second neighbours which have 1 index=0 in 3D
           ! (Only in 3D)
           {^IFTHREED
           i1=0;
           do i2=-1,1,2
             do i3=-1,1,2
               if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
             end do
           end do
           i2=0;
           do i3=-1,1,2
             do i1=-1,1,2
               if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
             end do
           end do
           i3=0;
           do i1=-1,1,2
             do i2=-1,1,2
               if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
             end do
           end do
           }
           ! Finally, the corners, that have no index=0
          {do i^D=-1,1,2\}
             if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
          {end do\}
         end if
       end if
    end do
    
    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%prim)
    end do

    if(bcphys.and.stagger_grid) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>psb(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idims=1,ndim
          ! to avoid using as yet unknown corner info in more than 1D, we
          ! fill only interior mesh ranges of the ghost cell ranges at first,
          ! and progressively enlarge the ranges to include corners later
          kmin1=0; kmax1=0;
          {^IFTWOD
           kmin2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0,-1,igrid)==1)
           kmax2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0, 1,igrid)==1)}
          {^IFTHREED
           kmin2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0,-1,0,igrid)==1)
           kmax2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0, 1,0,igrid)==1)
           kmin3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0,-1,igrid)==1)
           kmax3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0, 1,igrid)==1)}
          ixBmin^D=ixGlo^D+kmin^D*nghostcells;
          ixBmax^D=ixGhi^D-kmax^D*nghostcells;
          do iside=1,2
            i^D=kr(^D,idims)*(2*iside-3);
            if (aperiodB(idims)) then 
              if (neighbor_type(i^D,igrid) /= neighbor_boundary .and. &
                 .not. psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
            else 
              if (neighbor_type(i^D,igrid) /= neighbor_boundary) cycle
            end if
            call bc_phys(iside,idims,time,qdt,psb(igrid),ixG^L,ixB^L)
          end do
        end do
      end do
    end if

     ! modify normal component of magnetic field to fix divB=0 
    if(bcphys.and.associated(phys_boundary_adjust)) call phys_boundary_adjust()
    
    time_bc=time_bc+(MPI_WTIME()-time_bcin)
    
    contains

      logical function skip_direction(dir)
        integer, intent(in) :: dir(^ND)

        if (all(dir == 0)) then
           skip_direction = .true.
        else if (.not. req_diagonal .and. count(dir /= 0) > 1) then
           skip_direction = .true.
        else
           skip_direction = .false.
        end if
      end function skip_direction

      !> Receive from sibling at same refinement level
      subroutine bc_recv_srl

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           irecv_c=irecv_c+1
           itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
           call MPI_IRECV(psb(igrid)%prim,1,type_recv_srl(iib^D,i^D), &
                          ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
           if(stagger_grid) then
             irecv_srl=irecv_srl+1
             call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),sizes_srl_recv_total(i^D),MPI_DOUBLE_PRECISION, &
                            ipe_neighbor,itag,icomm,recvrequest_srl(irecv_srl),ierrmpi)
             ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i^D)
           end if
        end if

      end subroutine bc_recv_srl

      !> Receive from fine neighbor
      subroutine bc_recv_restrict

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
           if (ipe_neighbor/=mype) then
              irecv_c=irecv_c+1
              itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
              call MPI_IRECV(psb(igrid)%prim,1,type_recv_r(iib^D,inc^D), &
                             ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
              if(stagger_grid) then
                irecv_r=irecv_r+1
                call MPI_IRECV(recvbuffer_r(ibuf_recv_r),sizes_r_recv_total(inc^D), &
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                               icomm,recvrequest_r(irecv_r),ierrmpi)
                ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc^D)
              end if
           end if
        {end do\}

      end subroutine bc_recv_restrict

      !> Send to sibling at same refinement level
      subroutine bc_send_srl

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)

        if (ipole==0) then
           n_i^D=-i^D;
           if (ipe_neighbor==mype) then
              ixS^L=ixS_srl_^L(iib^D,i^D);
              ixR^L=ixR_srl_^L(iib^D,n_i^D);
              psb(ineighbor)%prim(ixR^S,nwhead:nwtail)=&
                  psb(igrid)%prim(ixS^S,nwhead:nwtail)
           else
              isend_c=isend_c+1
              itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
              call MPI_ISEND(psb(igrid)%prim,1,type_send_srl(iib^D,i^D), &
                             ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
              if(stagger_grid) then
                ibuf_start=ibuf_send_srl
                do idir=1,ndim
                  ixS^L=ixS_srl_stg_^L(idir,i^D);
                  ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i^D)
                  shapes=(/sizes_srl_send_stg(idir,i^D)/)
                  sendbuffer_srl(ibuf_start:ibuf_next-1)=&
                    reshape(psb(igrid)%prims(ixS^S,idir),shapes)
                  ibuf_start=ibuf_next
                end do
                isend_srl=isend_srl+1
                call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),sizes_srl_send_total(i^D),&
                  MPI_DOUBLE_PRECISION, ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
                ibuf_send_srl=ibuf_next
              end if
           end if
        else
           ixS^L=ixS_srl_^L(iib^D,i^D);
           select case (ipole)
           {case (^D)
              n_i^D=i^D^D%n_i^DD=-i^DD;\}
           end select
           if (ipe_neighbor==mype) then
              ixR^L=ixR_srl_^L(iib^D,n_i^D);
              call pole_copy(psb(ineighbor)%prim,ixG^L,ixR^L,psb(igrid)%prim,ixG^L,ixS^L)
           else
              if (isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)), &
                               sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%prim)
              end if
              allocate(pwbuf(ipwbuf)%prim(ixS^S,nwhead:nwtail))
              call pole_buffer(pwbuf(ipwbuf)%prim,ixS^L,ixS^L,psb(igrid)%prim,ixG^L,ixS^L)
              isend_c=isend_c+1
              isend_buf(ipwbuf)=isend_c
              itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
              isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
              call MPI_ISEND(pwbuf(ipwbuf)%prim,isizes,MPI_DOUBLE_PRECISION, &
                             ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
              ipwbuf=1+modulo(ipwbuf,npwbuf)
              if(stagger_grid) then
                ibuf_start=ibuf_send_srl
                do idir=1,ndim
                  ixS^L=ixS_srl_stg_^L(idir,i^D);
                  ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i^D)
                  shapes=(/sizes_srl_send_stg(idir,i^D)/)
                  sendbuffer_srl(ibuf_start:ibuf_next-1)=&
                    reshape(psb(igrid)%prims(ixS^S,idir),shapes)
                  ibuf_start=ibuf_next
                end do
                isend_srl=isend_srl+1
                call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),sizes_srl_send_total(i^D),&
                  MPI_DOUBLE_PRECISION, ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
                ibuf_send_srl=ibuf_next
              end if
           end if
        end if

      end subroutine bc_send_srl

      !> Send to coarser neighbor
      subroutine bc_send_restrict
        integer :: ii^D

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return
        if(phyboundblock(igrid).and..not.stagger_grid.and.bcphys) then
          ! to use block in physical boundary setup for coarse representative
          block=>psc(igrid)
          ! filling physical boundary ghost cells of a coarser representative block for
          ! sending swap region with width of nghostcells to its coarser neighbor
          do idims=1,ndim
             ! to avoid using as yet unknown corner info in more than 1D, we
             ! fill only interior mesh ranges of the ghost cell ranges at first,
             ! and progressively enlarge the ranges to include corners later
             {kmin^D=merge(0, 1, idims==^D)
             kmax^D=merge(0, 1, idims==^D)
             ixBmin^D=ixCoGmin^D+kmin^D*nghostcells
             ixBmax^D=ixCoGmax^D-kmax^D*nghostcells\}
             {^IFTWOD
             if(idims > 1 .and. neighbor_type(-1,0,igrid)==neighbor_boundary) ixBmin1=ixCoGmin1
             if(idims > 1 .and. neighbor_type( 1,0,igrid)==neighbor_boundary) ixBmax1=ixCoGmax1}
             {^IFTHREED
             if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==neighbor_boundary) ixBmin1=ixCoGmin1
             if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==neighbor_boundary) ixBmax1=ixCoGmax1
             if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==neighbor_boundary) ixBmin2=ixCoGmin2
             if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==neighbor_boundary) ixBmax2=ixCoGmax2}
             {if(i^D==-1) then
               ixBmin^D=ixCoGmin^D+nghostcells
               ixBmax^D=ixCoGmin^D+2*nghostcells-1
             else if(i^D==1) then
               ixBmin^D=ixCoGmax^D-2*nghostcells+1
               ixBmax^D=ixCoGmax^D-nghostcells
             end if\}
             do iside=1,2
                ii^D=kr(^D,idims)*(2*iside-3);
                if ({abs(i^D)==1.and.abs(ii^D)==1|.or.}) cycle
                if (neighbor_type(ii^D,igrid)/=neighbor_boundary) cycle
                call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoG^L,ixB^L)
             end do
          end do
        end if

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)

        if (ipole==0) then
           n_inc^D=-2*i^D+ic^D;
           if (ipe_neighbor==mype) then
              ixS^L=ixS_r_^L(iib^D,i^D);
              ixR^L=ixR_r_^L(iib^D,n_inc^D);
              psb(ineighbor)%prim(ixR^S,nwhead:nwtail)=&
               psc(igrid)%prim(ixS^S,nwhead:nwtail)
           else
              isend_c=isend_c+1
              itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
              call MPI_ISEND(psc(igrid)%prim,1,type_send_r(iib^D,i^D), &
                             ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
              if(stagger_grid) then 
                ibuf_start=ibuf_send_r
                do idir=1,ndim
                  ixS^L=ixS_r_stg_^L(idir,i^D);
                  ibuf_next=ibuf_start+sizes_r_send_stg(idir,i^D)
                  shapes=(/sizes_r_send_stg(idir,i^D)/)
                  sendbuffer_r(ibuf_start:ibuf_next-1)=&
                    reshape(psc(igrid)%prims(ixS^S,idir),shapes)
                  ibuf_start=ibuf_next
                end do
                isend_r=isend_r+1
                call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i^D),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                               icomm,sendrequest_r(isend_r),ierrmpi)
                ibuf_send_r=ibuf_next
              end if
           end if
        else
           ixS^L=ixS_r_^L(iib^D,i^D);
           select case (ipole)
           {case (^D)
              n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
           end select
           if (ipe_neighbor==mype) then
              ixR^L=ixR_r_^L(iib^D,n_inc^D);
              call pole_copy(psb(ineighbor)%prim,ixG^L,ixR^L,psc(igrid)%prim,ixCoG^L,ixS^L)
           else
              if (isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)), &
                               sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%prim)
              end if
              allocate(pwbuf(ipwbuf)%prim(ixS^S,nwhead:nwtail))
              call pole_buffer(pwbuf(ipwbuf)%prim,ixS^L,ixS^L,psc(igrid)%prim,ixCoG^L,ixS^L)
              isend_c=isend_c+1
              isend_buf(ipwbuf)=isend_c
              itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
              isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
              call MPI_ISEND(pwbuf(ipwbuf)%prim,isizes,MPI_DOUBLE_PRECISION, &
                             ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
              ipwbuf=1+modulo(ipwbuf,npwbuf)
              if(stagger_grid) then 
                ibuf_start=ibuf_send_r
                do idir=1,ndim
                  ixS^L=ixS_r_stg_^L(idir,i^D);
                  ibuf_next=ibuf_start+sizes_r_send_stg(idir,i^D)
                  shapes=(/sizes_r_send_stg(idir,i^D)/)
                  sendbuffer_r(ibuf_start:ibuf_next-1)=&
                    reshape(psc(igrid)%prims(ixS^S,idir),shapes)
                  ibuf_start=ibuf_next
                end do
                isend_r=isend_r+1
                call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i^D),&
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                               icomm,sendrequest_r(isend_r),ierrmpi)
                ibuf_send_r=ibuf_next
              end if
           end if
        end if

      end subroutine bc_send_restrict


      !> fill siblings ghost cells with received data
      subroutine bc_fill_srl
        double precision :: tmp(ixGs^T)
        integer :: ixS^L,ixR^L,n_i^D,ixSsync^L,ixRsync^L
        integer :: idir, idirect

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)
        ipole=neighbor_pole(i^D,igrid)
        idirect={abs(i^D)|+}

        !! Now the special treatment of the pole is done here, at the receive step
        if (ipole==0) then    
          ixR^L=ixR_srl_^L(iib^D,i^D);
          if (ipe_neighbor==mype) then
            n_i^D=-i^D;
            do idir=1,ndim
              ixS^L=ixS_srl_stg_^L(idir,n_i^D);
              ixR^L=ixR_srl_stg_^L(idir,i^D);
              !if (idirect == 1 .and. qdt>0.d0) then
              !  ! use the same value at the face shared by two neighors
              !  call indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L)
              !   if({idir==^D .and. i^D /=0  | .or.}) write(*,*) it,'difference',maxval(abs(psb(igrid)%prims(ixRsync^S,idir)-psb(ineighbor)%prims(ixSsync^S,idir)))
              !  !psb(igrid)%prims(ixRsync^S,idir) = half*(psb(igrid)%prims(ixRsync^S,idir)+psb(ineighbor)%prims(ixSsync^S,idir))
              !end if
              psb(igrid)%prims(ixR^S,idir) = psb(ineighbor)%prims(ixS^S,idir)
            end do
          else
            !! Unpack the buffer and fill the ghost cells
            n_i^D=-i^D;
            do idir=1,ndim
              ixS^L=ixS_srl_stg_^L(idir,n_i^D);
              ixR^L=ixR_srl_stg_^L(idir,i^D);
              ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
              tmp(ixS^S) = reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),shape=shape(psb(igrid)%prims(ixS^S,idir)))       
              !if (idirect == 1) then
              !   ! ixR ixS maybe changed
              !   call indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L) ! Overwrites ixR, ixS
              !   if(qdt==0) psb(igrid)%prims(ixRsync^S,idir) = half*(tmp(ixSsync^S) + psb(igrid)%prims(ixRsync^S,idir))
              !   if({idir==^D .and. i^D /=0  | .or.}) write(*,*) it,'betweenpe',maxval(abs(psb(igrid)%prims(ixRsync^S,idir)-tmp(ixSsync^S)))
              !end if
              psb(igrid)%prims(ixR^S,idir) = tmp(ixS^S)
              ibuf_recv_srl=ibuf_next
            end do
          end if

        else ! There is a pole
          select case (ipole)
          {case (^D)
             n_i^D=i^D^D%n_i^DD=-i^DD;\}
          end select
          if (ipe_neighbor==mype) then
            !! Fill ghost cells
            do idir=1,ndim
              ixR^L=ixR_srl_stg_^L(idir,i^D);
              ixS^L=ixS_srl_stg_^L(idir,n_i^D);
              !! Fill ghost cells
            end do
          else
             pole_buf%prims=zero
             do idir=1,ndim
              ixR^L=ixR_srl_stg_^L(idir,i^D);
              ixS^L=ixS_srl_stg_^L(idir,n_i^D);
              ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
              pole_buf%prims(ixS^S,idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
                shape=shape(psb(igrid)%prims(ixS^S,idir)))
              ibuf_recv_srl=ibuf_next
             end do
          end if
        end if

      end subroutine bc_fill_srl

      subroutine indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L)
        integer, intent(in)       :: i^D,idir
        integer, intent(inout)    :: ixR^L,ixS^L
        integer, intent(out)      :: ixRsync^L,ixSsync^L
      
        ixRsync^L=ixR^L;
        ixSsync^L=ixS^L;
        
        {
        if (i^D == -1 .and. idir == ^D) then
           ixRsyncmin^D = ixRmax^D
           ixRsyncmax^D = ixRmax^D
           ixSsyncmin^D = ixSmax^D
           ixSsyncmax^D = ixSmax^D
           ixRmax^D = ixRmax^D - 1
           ixSmax^D = ixSmax^D - 1
        else if (i^D == 1 .and. idir == ^D) then
           ixRsyncmin^D = ixRmin^D
           ixRsyncmax^D = ixRmin^D
           ixSsyncmin^D = ixSmin^D
           ixSsyncmax^D = ixSmin^D
           ixRmin^D = ixRmin^D + 1
           ixSmin^D = ixSmin^D + 1
        end if
        \}

      end subroutine indices_for_syncing

      !> fill restricted ghost cells after receipt
      subroutine bc_fill_r

        ipole=neighbor_pole(i^D,igrid)
        if (ipole==0) then
          ! Loop over the children ic^D to and their neighbors inc^D
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             n_i^D=-i^D;
             ineighbor=neighbor_child(1,inc^D,igrid)
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             if (ipe_neighbor==mype) then ! Same processor
               do idir=1,ndim
                  ixS^L=ixS_r_stg_^L(idir,n_i^D);
                  ixR^L=ixR_r_stg_^L(idir,inc^D);
                  psb(igrid)%prims(ixR^S,idir)=psc(ineighbor)%prims(ixS^S,idir)
               end do
             else ! Different processor
               !! Unpack the buffer and fill the ghost cells
               do idir=1,ndim
                 ixR^L=ixR_r_stg_^L(idir,inc^D);
                 ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
                 psb(igrid)%prims(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                       shape=shape(psb(igrid)%prims(ixR^S,idir)))
                 ibuf_recv_r=ibuf_next
               end do
             end if
          {end do\}
        
        else !! There is a pole
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             select case(ipole)
            {case (^D)
               n_i^D=i^D^D%n_i^DD=-i^DD;\}
             end select
             ineighbor=neighbor_child(1,inc^D,igrid)
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             ixR^L=ixR_r_^L(iib^D,inc^D);
             if (ipe_neighbor==mype) then ! Same processor
               do idir=1,ndim
                 ixS^L=ixS_r_stg_^L(idir,n_i^D);
                 ixR^L=ixR_r_stg_^L(idir,inc^D);
                 !! Fill ghost cells
               end do
             else ! Different processor
               !! Unpack the buffer and fill an auxiliary array
               pole_buf%prims=zero
               do idir=1,ndim
                 ixS^L=ixS_r_stg_^L(idir,n_i^D);
                 ixR^L=ixR_r_stg_^L(idir,inc^D);
                 ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
                 pole_buf%prims(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                   shape=shape(psb(igrid)%prims(ixR^S,idir)))
                 ibuf_recv_r=ibuf_next
               end do
             end if
          {end do\}
        
        end if

      end subroutine bc_fill_r

      !> Receive from coarse neighbor
      subroutine bc_recv_prolong

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           irecv_c=irecv_c+1
           inc^D=ic^D+i^D;
           itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
           call MPI_IRECV(psc(igrid)%prim,1,type_recv_p(iib^D,inc^D), &
                          ipe_neighbor,itag,icomm,recvrequest_c_p(irecv_c),ierrmpi)  
           if(stagger_grid) then
             irecv_p=irecv_p+1
             call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc^D),& 
                            MPI_DOUBLE_PRECISION,ipe_neighbor,itag,&
                            icomm,recvrequest_p(irecv_p),ierrmpi)
             ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc^D)
           end if
        end if

      end subroutine bc_recv_prolong

      !> Send to finer neighbor
      subroutine bc_send_prolong
        integer :: ii^D

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ixS^L=ixS_p_^L(iib^D,inc^D);

           ineighbor=neighbor_child(1,inc^D,igrid)
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
        
           if (ipole==0) then
              n_i^D=-i^D;
              n_inc^D=ic^D+n_i^D;
              if (ipe_neighbor==mype) then
                ixR^L=ixR_p_^L(iib^D,n_inc^D);
                psc(ineighbor)%prim(ixR^S,nwhead:nwtail) &
                   =psb(igrid)%prim(ixS^S,nwhead:nwtail)
              else
                isend_c=isend_c+1
                itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
                call MPI_ISEND(psb(igrid)%prim,1,type_send_p(iib^D,inc^D), &
                               ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
                if(stagger_grid) then 
                  ibuf_start=ibuf_send_p
                  do idir=1,ndim
                    ixS^L=ixS_p_stg_^L(idir,inc^D);
                    ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc^D)
                    shapes=(/sizes_p_send_stg(idir,inc^D)/)
                    sendbuffer_p(ibuf_start:ibuf_next-1)=&
                      reshape(psb(igrid)%prims(ixS^S,idir),shapes)   
                    ibuf_start=ibuf_next
                  end do
                  isend_p=isend_p+1
                  call MPI_ISEND(sendbuffer_p(ibuf_send_p),sizes_p_send_total(inc^D),&
                                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                                 icomm,sendrequest_p(isend_p),ierrmpi)
                  ibuf_send_p=ibuf_next
                end if
              end if
           else
              select case (ipole)
              {case (^D)
                 n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;\}
              end select
              if (ipe_neighbor==mype) then
                 ixR^L=ixR_p_^L(iib^D,n_inc^D);
                 call pole_copy(psc(ineighbor)%prim,ixCoG^L,ixR^L,psb(igrid)%prim,ixG^L,ixS^L)
                if(stagger_grid) then
                  do idir=1,ndim
                    ixS^L=ixS_p_stg_^L(idir,inc^D);
                    ixR^L=ixR_p_stg_^L(idir,n_inc^D);
                  end do
                end if
              else
                if (isend_buf(ipwbuf)/=0) then
                   call MPI_WAIT(sendrequest_c_p(isend_buf(ipwbuf)), &
                                 sendstatus_c_p(:,isend_buf(ipwbuf)),ierrmpi)
                   deallocate(pwbuf(ipwbuf)%prim)
                end if
                allocate(pwbuf(ipwbuf)%prim(ixS^S,nwhead:nwtail))
                call pole_buffer(pwbuf(ipwbuf)%prim,ixS^L,ixS^L,psb(igrid)%prim,ixG^L,ixS^L)
                isend_c=isend_c+1
                isend_buf(ipwbuf)=isend_c
                itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
                isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
                call MPI_ISEND(pwbuf(ipwbuf)%prim,isizes,MPI_DOUBLE_PRECISION, &
                               ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
                ipwbuf=1+modulo(ipwbuf,npwbuf)
                if(stagger_grid) then 
                  ibuf_start=ibuf_send_p
                  do idir=1,ndim
                    ixS^L=ixS_p_stg_^L(idir,inc^D);
                    ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc^D)
                    shapes=(/sizes_p_send_stg(idir,inc^D)/)
                    sendbuffer_p(ibuf_start:ibuf_next-1)=&
                      reshape(psb(igrid)%prims(ixS^S,idir),shapes)   
                    ibuf_start=ibuf_next
                  end do
                  isend_p=isend_p+1
                  call MPI_ISEND(sendbuffer_p(ibuf_send_p),sizes_p_send_total(inc^D),&
                                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                                 icomm,sendrequest_p(isend_p),ierrmpi)
                  ibuf_send_p=ibuf_next
                end if
              end if
           end if
        {end do\}

      end subroutine bc_send_prolong

      !> fill coarser representative with data from coarser neighbors
      subroutine bc_fill_p
        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)
        ipole=neighbor_pole(i^D,igrid)

        if (ipole==0) then   !! There is no pole 
          inc^D=ic^D+i^D;
          ixR^L=ixR_p_^L(iib^D,inc^D);
          if(ipe_neighbor==mype) then !! Same processor
            n_inc^D=-2*i^D+ic^D;
            do idir=1,ndim
              ixS^L=ixS_p_stg_^L(idir,n_inc^D);
              ixR^L=ixR_p_stg_^L(idir,inc^D);
              psc(igrid)%prims(ixR^S,idir)=psb(ineighbor)%prims(ixS^S,idir)
            end do
          else !! Different processor
            do idir=1,ndim
              ixR^L=ixR_p_stg_^L(idir,inc^D);
              ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
              psc(igrid)%prims(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                    shape=shape(psc(igrid)%prims(ixR^S,idir)))
              ibuf_recv_p=ibuf_next
            end do
          end if

        else !! There is a pole
          inc^D=ic^D+i^D;
          select case (ipole)
          {case (^D)
             n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
          end select
          if (ipe_neighbor==mype) then
            do idir=1,ndim
              ixS^L=ixS_p_stg_^L(idir,n_inc^D);
              ixR^L=ixR_p_stg_^L(idir,inc^D);
            end do
          else
            !! Unpack the buffer and fill an auxiliary array
            pole_buf%prims=zero
            do idir=1,ndim
              ixR^L=ixR_p_stg_^L(idir,inc^D);
              ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
              pole_buf%prims(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                shape=shape(psc(igrid)%prims(ixR^S,idir)))
              ibuf_recv_p=ibuf_next
            end do
          end if

        end if

      end subroutine bc_fill_p

      !> do prolongation for fine blocks after receipt data from coarse neighbors
      subroutine bc_prolong

        integer :: ixFi^L,ixCo^L,ii^D
        double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D

        ixFi^L=ixR_srl_^L(iib^D,i^D);
        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        ! compute the enlarged grid lower left corner coordinates
        ! these are true coordinates for an equidistant grid, 
        ! but we can temporarily also use them for getting indices 
        ! in stretched grids
        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;

        if(stagger_grid.and.phyboundblock(igrid).and.bcphys) then
          block=>psc(igrid)
          do idims=1,ndim
            ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
            ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;
            {^IFTHREED
            ! avoid using undetermined ghost cells at physical boundary edges
            if(idims == 1) then
              if(neighbor_type(-1,0,0,igrid)==neighbor_boundary .or. &
                 neighbor_type(1,0,0,igrid)==neighbor_boundary) then
                if(neighbor_type(0,-1,0,igrid)==neighbor_boundary) ixComin2=ixCoMmin2
                if(neighbor_type(0,0,-1,igrid)==neighbor_boundary) ixComin3=ixCoMmin3
                if(neighbor_type(0,1,0,igrid)==neighbor_boundary) ixComax2=ixCoMmax2
                if(neighbor_type(0,0,1,igrid)==neighbor_boundary) ixComax3=ixCoMmax3
              end if
            else if(idims == 2) then
              if(neighbor_type(0,-1,0,igrid)==neighbor_boundary .or. &
                 neighbor_type(0,1,0,igrid)==neighbor_boundary) then
                if(neighbor_type(0,0,-1,igrid)==neighbor_boundary) ixComin3=ixCoMmin3
                if(neighbor_type(0,0,1,igrid)==neighbor_boundary) ixComax3=ixCoMmax3
              end if
            end if
            }
            do iside=1,2
              ii^D=kr(^D,idims)*(2*iside-3);
              if(neighbor_type(ii^D,igrid)/=neighbor_boundary) cycle
              if(( {(iside==1.and.idims==^D.and.ixComin^D<ixCoGmin^D+nghostcells)|.or.} ) &
               .or.( {(iside==2.and.idims==^D.and.ixComax^D>ixCoGmax^D-nghostcells)|.or. })) then
                {ixBmin^D=merge(ixCoGmin^D,ixComin^D,idims==^D);}
                {ixBmax^D=merge(ixCoGmax^D,ixComax^D,idims==^D);}
                call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoG^L,ixB^L)
              end if
            end do
          end do
        end if

!        if(prolongprimitive) then
!           ! following line again assumes equidistant grid, but 
!           ! just computes indices, so also ok for stretched case
!           ! reason for +1-1 and +1+1: the coarse representation has 
!           ! also nghostcells at each side. During
!           ! prolongation, we need cells to left and right, hence -1/+1
!           ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
!           ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;
!           call phys_to_primitive(ixCoG^L,ixCo^L,&
!             psc(igrid)%prim,psc(igrid)%x)
!        endif

        select case (typeghostfill)
        case ("linear")
           call interpolation_linear(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case ("copy")
           call interpolation_copy(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case default
           write (unitterm,*) "Undefined typeghostfill ",typeghostfill
           call mpistop("Undefined typeghostfill")
        end select

      end subroutine bc_prolong

      subroutine bc_prolong_stg(NeedProlong)
        use mod_amr_fct
        logical,dimension(-1:1^D&) :: NeedProlong
        logical                    :: fine_^Lin
        integer                    :: ixFi^L,ixCo^L
        double precision           :: dxFi^D,dxCo^D,xFimin^D,xComin^D,invdxCo^D
        ! Check what is already at the desired level
        fine_^Lin=.false.;
        {
        if(i^D>-1) fine_min^Din=(.not.NeedProlong(i^DD-kr(^D,^DD)).and.neighbor_type(i^DD-kr(^D,^DD),igrid)/=1)
        if(i^D<1)  fine_max^Din=(.not.NeedProlong(i^DD+kr(^D,^DD)).and.neighbor_type(i^DD+kr(^D,^DD),igrid)/=1)
        \}

        ixFi^L=ixR_srl_^L(iib^D,i^D);

        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;

        ! moved the physical boundary filling here, to only fill the
        ! part needed

        ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
        ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

        call prolong_2nd_stg(psc(igrid),psb(igrid),ixCo^L,ixFi^L,dxCo^D,xComin^D,dxFi^D,xFimin^D,.true.,fine_^Lin)

        ! The current region has already been refined, so it does not need to be prolonged again
        NeedProlong(i^D)=.false. 

      end subroutine bc_prolong_stg

      subroutine interpolation_linear(ixFi^L,dxFi^D,xFimin^D, &
                                      dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims, nwmin,nwmax
        double precision :: xCo^D, xFi^D, eta^D
        double precision :: slopeL, slopeR, slopeC, signC, signR
        double precision :: slope(1:nprim,ndim)
        !!double precision :: local_invdxCo^D
        double precision :: signedfactorhalf^D
        !integer :: ixshift^D, icase

        !icase=mod(nghostcells,2)

        nwmin=1
        nwmax=nprim

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           ! here we temporarily use an equidistant grid
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB
        
           ! indices of coarse cell which contains the fine cell
           ! since we computed lower left corner earlier 
           ! in equidistant fashion: also ok for stretched case
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1
        
           ! cell-centered coordinates of coarse grid point
           ! here we temporarily use an equidistant grid
           xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB \}

           !if(.not.slab) then
           !   ^D&local_invdxCo^D=1.d0/psc(igrid)%dx({ixCo^DD},^D)\
           !endif

           if(slab_uniform) then
             ! actual cell-centered coordinates of fine grid point
             !!^D&xFi^D=block%x({ixFi^DD},^D)\
             ! actual cell-centered coordinates of coarse grid point
             !!^D&xCo^D=psc(igrid)%x({ixCo^DD},^D)\
             ! normalized distance between fine/coarse cell center
             ! in coarse cell: ranges from -0.5 to 0.5 in each direction
             ! (origin is coarse cell center)
             ! this is essentially +1/4 or -1/4 on cartesian mesh
             eta^D=(xFi^D-xCo^D)*invdxCo^D;
           else
             !select case(icase)
             ! case(0)
             !{! here we assume an even number of ghostcells!!!
             !ixshift^D=2*(mod(ixFi^D,2)-1)+1
             !if(ixshift^D>0.0d0)then
             !   ! oneven fine grid points
             !   eta^D=-0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D:ixFi^D+1^D%ixFi^DD))) 
             !else
             !   ! even fine grid points
             !   eta^D=+0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D-1:ixFi^D^D%ixFi^DD))) 
             !endif\}
             ! case(1)
             !{! here we assume an odd number of ghostcells!!!
             !ixshift^D=2*(mod(ixFi^D,2)-1)+1
             !if(ixshift^D>0.0d0)then
             !   ! oneven fine grid points
             !   eta^D=+0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D-1:ixFi^D^D%ixFi^DD))) 
             !else
             !   ! even fine grid points
             !   eta^D=-0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D:ixFi^D+1^D%ixFi^DD))) 
             !endif\}
             ! case default
             !  call mpistop("no such case")
             !end select
             ! the different cases for even/uneven number of ghost cells 
             ! are automatically handled using the relative index to ixMlo
             ! as well as the pseudo-coordinates xFi and xCo 
             ! these latter differ from actual cell centers when stretching is used
             ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
             {if(xFi^D>xCo^D) then
                signedfactorhalf^D=0.5d0
              else
                signedfactorhalf^D=-0.5d0
              end if
              eta^D=signedfactorhalf^D*(one-psb(igrid)%dvolume(ixFi^DD) &
                   /sum(psb(igrid)%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
             !{eta^D=(xFi^D-xCo^D)*invdxCo^D &
             !      *two*(one-block%dvolume(ixFi^DD) &
             !      /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
           end if
        
           do idims=1,ndim
              hxCo^D=ixCo^D-kr(^D,idims)\
              jxCo^D=ixCo^D+kr(^D,idims)\
        
              do iw=nwmin,nwmax
                 slopeL=psc(igrid)%prim(ixCo^D,iw)-psc(igrid)%prim(hxCo^D,iw)
                 slopeR=psc(igrid)%prim(jxCo^D,iw)-psc(igrid)%prim(ixCo^D,iw)
                 slopeC=half*(slopeR+slopeL)
        
                 ! get limited slope
                 signR=sign(one,slopeR)
                 signC=sign(one,slopeC)
                 select case(typeprolonglimit)
                 case('unlimit')
                   slope(iw,idims)=slopeC
                 case('minmod')
                   slope(iw,idims)=signR*max(zero,min(dabs(slopeR), &
                                                     signR*slopeL))
                 case('woodward')
                   slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), &
                                      signR*slopeL,signR*half*slopeC))
                 case('koren')
                   slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, &
                    (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
                 case default
                   slope(iw,idims)=signC*max(zero,min(dabs(slopeC), &
                                     signC*slopeL,signC*slopeR))
                 end select
              end do
           end do
        
           ! Interpolate from coarse cell using limited slopes
           psb(igrid)%prim(ixFi^D,nwmin:nwmax)=psc(igrid)%prim(ixCo^D,nwmin:nwmax)+&
             {(slope(nwmin:nwmax,^D)*eta^D)+}
        
        {end do\}
        
      end subroutine interpolation_linear

      subroutine interpolation_copy(ixFi^L,dxFi^D,xFimin^D, &
                                    dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        integer :: ixCo^D, ixFi^D, nwmin,nwmax
        double precision :: xFi^D

        nwmin=1
        nwmax=nprim

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB
        
           ! indices of coarse cell which contains the fine cell
           ! note: this also works for stretched grids
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1\}
        
           ! Copy from coarse cell
           psb(igrid)%prim(ixFi^D,nwmin:nwmax)=psc(igrid)%prim(ixCo^D,nwmin:nwmax)
        
        {end do\}
        
      end subroutine interpolation_copy

      subroutine pole_copy(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L
        double precision :: wrecv(ixIR^S,1:nprim), wsend(ixIS^S,1:nprim)

        integer :: iw, iB

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case ("symm")
               wrecv(ixR^S,iw) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case ("asymm")
               wrecv(ixR^S,iw) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do \}
        end select
      
      end subroutine pole_copy

      subroutine pole_buffer(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L
        double precision :: wrecv(ixIR^S,nwhead:nwtail), wsend(ixIS^S,1:nprim)

        integer :: iw, iB

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case ("symm")
               wrecv(ixR^S,iw) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case ("asymm")
               wrecv(ixR^S,iw) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do \}
        end select
      
      end subroutine pole_buffer

  end subroutine getbc

  subroutine identifyphysbound(s,iib^D)
    use mod_global_parameters

    type(state)          :: s
    integer, intent(out) :: iib^D

    {
    if(s%is_physical_boundary(2*^D) .and. &
       s%is_physical_boundary(2*^D-1)) then
      iib^D=2
    else if(s%is_physical_boundary(2*^D-1)) then
      iib^D=-1
    else if(s%is_physical_boundary(2*^D)) then
      iib^D=1
    else
      iib^D=0
    end if
    \}

  end subroutine identifyphysbound

end module mod_ghostcells_update
