!> Get first available igrid on processor ipe
integer function getnode(ipe)
  use mod_forest, only: igrid_inuse
  use mod_global_parameters

  integer, intent(in) :: ipe
  integer :: igrid, igrid_available

  igrid_available=0

  do igrid=1,max_blocks
     if (igrid_inuse(igrid,ipe)) cycle

     igrid_available=igrid
     exit
  end do

  if (igrid_available == 0) then
     getnode = -1
     print *, "Current maximum number of grid blocks:", max_blocks
     call mpistop("Insufficient grid blocks; increase max_blocks in meshlist")
  else
     getnode=igrid_available
     igrid_inuse(igrid,ipe)=.true.
  end if

  if (ipe==mype) then
     ! initialize nodal block
     node(1:nodehi,getnode) = 0
     rnode(1:rnodehi,getnode) = zero
  end if

end function getnode

! put igrid on processor ipe to be not in use
subroutine putnode(igrid,ipe)
  use mod_forest
  implicit none

  integer, intent(in) :: igrid, ipe

  igrid_inuse(igrid,ipe)=.false.

end subroutine putnode

!> allocate arrays on igrid node
subroutine alloc_node(igrid)
  use mod_forest
  use mod_global_parameters
  use mod_geometry

  integer, intent(in) :: igrid

  integer :: level, ig1, ign1, ixCoGmin1,ixCoGmax1, ix, i1
  integer :: imin, imax, index, igCo1, ixshift, offset, ifirst
  integer:: icase, ixGextmin1,ixGextmax1
  double precision :: dx1, summeddx, sizeuniformpart1
  double precision :: xext(ixGlo1-1:ixGhi1+1,1:ndim)

  ixCoGmin1=1;
  ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells;

  icase=mod(nghostcells,2)
  if(stagger_grid) icase=1
  select case(icase)
    case(0)
      ixGextmin1=ixGlo1;ixGextmax1=ixGhi1;
    case(1)
      ! for ghost cell related prolongations, we need
      ! an extra layer with known volumes and dx-intervals
      ! in case the number of ghost cells is odd
      ixGextmin1=ixGlo1-1;ixGextmax1=ixGhi1+1;
    case default
      call mpistop("no such case")
  end select

  if(.not. allocated(ps(igrid)%prim)) then

    ! allocate arrays for solution and space
    call alloc_state(igrid, ps(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
        .true.)
    ! allocate arrays for one level coarser solution
    call alloc_state(igrid, psc(igrid), ixCoGmin1,ixCoGmax1, ixCoGmin1,&
       ixCoGmax1, .true.)
    ! allocate arrays for old solution
    call alloc_state(igrid, pso(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
        .false.)
    ! allocate arrays for temp solution 1
    call alloc_state(igrid, ps1(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
        .false.)

    ! allocate temperary solution space
    select case (time_integrator)
    case("ssprk3","ssprk4","jameson","IMEX_Midpoint","IMEX_Trapezoidal")
      call alloc_state(igrid, ps2(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
          .false.)
    case("RK3_BT","rk4","ssprk5")
      call alloc_state(igrid, ps2(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
          .false.)
      call alloc_state(igrid, ps3(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
          .false.)
    case("IMEX_ARS3","IMEX_232")
      call alloc_state(igrid, ps2(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
          .false.)
      call alloc_state(igrid, ps3(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
          .false.)
      call alloc_state(igrid, ps4(igrid), ixGlo1,ixGhi1, ixGextmin1,ixGextmax1,&
          .false.)
    end select

  end if

  ps(igrid)%prim=0.d0
  ps1(igrid)%prim=0.d0
  psc(igrid)%prim=0.d0
  ps(igrid)%cons=0.d0
  ps1(igrid)%cons=0.d0
  psc(igrid)%cons=0.d0
  ps(igrid)%igrid=igrid

  if(stagger_grid) then
    ps(igrid)%prims=0.d0
    ps1(igrid)%prims=0.d0
    psc(igrid)%prims=0.d0
  end if


  ! block pointer to current block
  block=>ps(igrid)

  ! set level information
  level=igrid_to_node(igrid,mype)%node%level
  ig1=igrid_to_node(igrid,mype)%node%ig1;

  node(plevel_,igrid)=level
  node(pig1_,igrid)=ig1

  ! set dx information
  rnode(rpdx1_,igrid)=dx(1,level)
  dxlevel(:)=dx(:,level)

  ! uniform cartesian case as well as all unstretched coordinates
  ! determine the minimal and maximal corners
  rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
  rnode(rpxmax1_,igrid)=xprobmin1+dble(ig1)*dg1(level)
!  ^D&rnode(rpxmax^D_,igrid)=xprobmax^D-dble(ng^D(level)-ig^D)*dg^D(level)\

  dx1=rnode(rpdx1_,igrid)
 do ix=ixGlo1,ixMhi1-nghostcells
    ps(igrid)%x(ix,1)=rnode(rpxmin1_,igrid)+(dble(ix-nghostcells)-0.5d0)*dx1
  end do
 ! update overlap cells of neighboring blocks in the same way to get the same values
 do ix=ixMhi1-nghostcells+1,ixGhi1
    ps(igrid)%x(ix,1)=rnode(rpxmax1_,igrid)+(dble(ix-ixMhi1)-0.5d0)*dx1
  end do

  dx1=2.0d0*rnode(rpdx1_,igrid)
 do ix=ixCoGmin1,ixCoGmax1
    psc(igrid)%x(ix,1)=rnode(rpxmin1_,igrid)+(dble(ix-nghostcells)-0.5d0)*dx1
  end do

  ps(igrid)%dx(ixGextmin1:ixGextmax1,1)=rnode(rpdx1_,igrid);
  psc(igrid)%dx(ixCoGmin1:ixCoGmax1,1)=2.0d0*rnode(rpdx1_,igrid);
  dx1=rnode(rpdx1_,igrid)
 do ix=ixGextmin1,ixGextmax1
    xext(ix,1)=rnode(rpxmin1_,igrid)+(dble(ix-nghostcells)-0.5d0)*dx1
  end do

  if(any(stretched_dim)) then
   if(stretch_type(1) == stretch_uni)then
      imin=(ig1-1)*block_nx1
      imax=ig1*block_nx1
      rnode(rpxmin1_,igrid)=xprobmin1+dxfirst_1mq(level,&
         1) *(1.0d0-qstretch(level,1)**imin)
      rnode(rpxmax1_,igrid)=xprobmin1+dxfirst_1mq(level,&
         1) *(1.0d0-qstretch(level,1)**imax)
      ! fix possible out of bound due to precision
      if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
      ixshift=(ig1-1)*block_nx1-nghostcells
      do ix=ixGextmin1,ixGextmax1
        index=ixshift+ix
        ps(igrid)%dx(ix,1)=dxfirst(level,1)*qstretch(level,1)**(index-1)
      enddo
      igCo1=(ig1-1)/2
      ixshift=igCo1*block_nx1+(1-modulo(ig1,2))*block_nx1/2-nghostcells
      do ix=ixCoGmin1,ixCoGmax1
        index=ixshift+ix
        psc(igrid)%dx(ix,1)=dxfirst(level-1,1)*qstretch(level-1,1)**(index-1)
        psc(igrid)%x(ix,1)=xprobmin1+dxfirst_1mq(level-1,&
           1)*(1.0d0-qstretch(level-1,1)**(index-1))+ 0.5d0*dxfirst(level-1,&
           1)*qstretch(level-1,1)**(index-1)
      end do
      ! now that dx and grid boundaries are known: fill cell centers
      ifirst=nghostcells+1
      ! first fill the mesh
      summeddx=0.0d0
      do ix=ixMlo1,ixMhi1
        ps(igrid)%x(ix,1)=rnode(rpxmin1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
           1)
        summeddx=summeddx+ps(igrid)%dx(ix,1)
      enddo
      ! then ghost cells to left
      summeddx=0.0d0
      do ix=nghostcells,1,-1
        ps(igrid)%x(ix,1)=rnode(rpxmin1_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,&
           1)
        summeddx=summeddx+ps(igrid)%dx(ix,1)
      enddo
      ! then ghost cells to right
      summeddx=0.0d0
      do ix=ixGhi1-nghostcells+1,ixGhi1
        ps(igrid)%x(ix,1)=rnode(rpxmax1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
           1)
        summeddx=summeddx+ps(igrid)%dx(ix,1)
      enddo
      select case(icase)
        case(0)
          ! if even number of ghost cells: xext is just copy of local x
          xext(ixGextmin1:ixGextmax1,1)=ps(igrid)%x(ixGextmin1:ixGextmax1,1)
        case(1)
          ! if uneven number of ghost cells: extra layer left/right
          summeddx=0.0d0
          do ix=ixMlo1,ixMhi1
            xext(ix,1)=rnode(rpxmin1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,1)
           summeddx=summeddx+ps(igrid)%dx(ix,1)
          enddo
          ! then ghost cells to left
          summeddx=0.0d0
          do ix=nghostcells,ixGextmin1,-1
            xext(ix,1)=rnode(rpxmin1_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,1)
            summeddx=summeddx+ps(igrid)%dx(ix,1)
          enddo
          ! then ghost cells to right
          summeddx=0.0d0
          do ix=ixGhi1-nghostcells+1,ixGextmax1
             xext(ix,1)=rnode(rpxmax1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                1)
             summeddx=summeddx+ps(igrid)%dx(ix,1)
          enddo
        case default
          call mpistop("no such case")
      end select
     endif
    if(stretch_type(1) == stretch_symm)then
       ! here we distinguish three kinds of grid blocks
       ! depending on their ig-index, set per level 
       !      the first n_stretchedblocks/2  will stretch to the left
       !      the middle ntotal-n_stretchedblocks will be uniform
       !      the last  n_stretchedblocks/2  will stretch to the right
       if(ig1<=nstretchedblocks(level,1)/2)then
         ! stretch to the left
         offset=block_nx1*nstretchedblocks(level,1)/2
         imin=(ig1-1)*block_nx1
         imax=ig1*block_nx1
         rnode(rpxmin1_,igrid)=xprobmin1+xstretch1-dxfirst_1mq(level,&
            1) *(1.0d0-qstretch(level,1)**(offset-imin))
         rnode(rpxmax1_,igrid)=xprobmin1+xstretch1-dxfirst_1mq(level,&
            1) *(1.0d0-qstretch(level,1)**(offset-imax))
         ! fix possible out of bound due to precision
         if(rnode(rpxmin1_,igrid)<xprobmin1) rnode(rpxmin1_,igrid)=xprobmin1
         ixshift=(ig1-1)*block_nx1-nghostcells
         do ix=ixGextmin1,ixGextmax1
           index=ixshift+ix
           ps(igrid)%dx(ix,1)=dxfirst(level,1)*qstretch(level,&
              1)**(offset-index)
         enddo
         ixshift=(nstretchedblocks(level,&
            1)/2-ig1)*(block_nx1/2)+block_nx1/2+nghostcells
         do ix=ixCoGmin1,ixCoGmax1
           index=ixshift-ix
           psc(igrid)%dx(ix,1)=dxfirst(level-1,1)*qstretch(level-1,1)**index
         enddo
         ! last block: to modify ghost cells!!!
         if(ig1==nstretchedblocks(level,1)/2)then
           if(ng1(level)==nstretchedblocks(level,1))then
             ! if middle blocks do not exist then use symmetry
             do ix=ixGhi1-nghostcells+1,ixGextmax1
                ps(igrid)%dx(ix,1)= ps(igrid)%dx(2*(ixGhi1-nghostcells)+1-ix,&
                   1)
             enddo
             do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
                psc(igrid)%dx(ix,1)= psc(igrid)%dx(2*(ixCoGmax1-nghostcells)+&
                   1-ix,1)
             enddo
           else
             ! if middle blocks exist then use same as middle blocks: 
             do ix=ixGhi1-nghostcells+1,ixGextmax1
                ps(igrid)%dx(ix,1)=dxmid(level,1)
             enddo
             do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
                psc(igrid)%dx(ix,1)=dxmid(level-1,1)
             enddo
           endif
         endif
         ! first block: make ghost cells symmetric (to allow periodicity)
         if(ig1==1)then
           do ix=ixGextmin1,nghostcells
             ps(igrid)%dx(ix,1)=ps(igrid)%dx(2*nghostcells+1-ix,1)
           enddo
           do ix=1,nghostcells
             psc(igrid)%dx(ix,1)=psc(igrid)%dx(2*nghostcells+1-ix,1)
           enddo
         endif
       else 
         if(ig1<=ng1(level)-nstretchedblocks(level,1)/2) then
           ! keep uniform
           ps(igrid)%dx(ixGextmin1:ixGextmax1,1)=dxmid(level,1)
           psc(igrid)%dx(ixCoGmin1:ixCoGmax1,1)=dxmid(level-1,1)
           rnode(rpxmin1_,igrid)=xprobmin1+xstretch1+&
              (ig1-nstretchedblocks(level,1)/2-1)*block_nx1*dxmid(level,1)
           rnode(rpxmax1_,igrid)=xprobmin1+xstretch1+&
              (ig1-nstretchedblocks(level,1)/2)  *block_nx1*dxmid(level,1)
           ! first and last block: to modify the ghost cells!!!
           if(ig1==nstretchedblocks(level,1)/2+1)then
             do ix=ixGextmin1,nghostcells
               ps(igrid)%dx(ix,1)=dxfirst(level,1)*qstretch(level,&
                  1)**(nghostcells-ix)
             enddo
             do ix=1,nghostcells
               psc(igrid)%dx(ix,1)=dxfirst(level-1,1)*qstretch(level-1,&
                  1)**(nghostcells-ix)
             enddo
           endif
           if(ig1==ng1(level)-nstretchedblocks(level,1))then
             do ix=ixGhi1-nghostcells+1,ixGextmax1
               ps(igrid)%dx(ix,1)=dxfirst(level,1)*qstretch(level,&
                  1)**(ix-block_nx1-nghostcells-1)
             enddo
             do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
               psc(igrid)%dx(ix,1)=dxfirst(level-1,1)*qstretch(level-1,&
                  1)**(ix-ixCoGmax1+nghostcells-1)
             enddo
           endif
         else
           ! stretch to the right
           offset=block_nx1*(ng1(level)-nstretchedblocks(level,1)/2)
           sizeuniformpart1=dxmid(1,1)*(domain_nx1-&
              nstretchedblocks_baselevel(1)*block_nx1)
           imin=(ig1-1)*block_nx1-offset
           imax=ig1*block_nx1-offset
           rnode(rpxmin1_,igrid)=xprobmin1+xstretch1+sizeuniformpart1+&
              dxfirst_1mq(level,1) *(1.0d0-qstretch(level,1)**imin)
           rnode(rpxmax1_,igrid)=xprobmin1+xstretch1+sizeuniformpart1+&
              dxfirst_1mq(level,1) *(1.0d0-qstretch(level,1)**imax)
           ! fix possible out of bound due to precision
           if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
           ixshift=(ig1-1)*block_nx1-nghostcells-offset
           do ix=ixGextmin1,ixGextmax1
             index=ixshift+ix
             ps(igrid)%dx(ix,1)=dxfirst(level,1)*qstretch(level,1)**(index-1)
           enddo
           ixshift=(ig1+nstretchedblocks(level,&
              1)/2-ng1(level)-1)*(block_nx1/2)-nghostcells
           do ix=ixCoGmin1,ixCoGmax1
             index=ixshift+ix
             psc(igrid)%dx(ix,1)=dxfirst(level-1,1)*qstretch(level-1,&
                1)**(index-1)
           enddo
           ! first block: modify ghost cells!!!
           if(ig1==ng1(level)-nstretchedblocks(level,1)+1)then
             if(ng1(level)==nstretchedblocks(level,1))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGextmin1,nghostcells
                 ps(igrid)%dx(ix,1)=ps(igrid)%dx(2*nghostcells+1-ix,1)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ix,1)=psc(igrid)%dx(2*nghostcells+1-ix,1)
               enddo
             else
               ! if middle blocks exist then use same as middle blocks: 
               do ix=ixGextmin1,nghostcells
                 ps(igrid)%dx(ix,1)=dxmid(level,1)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ix,1)=dxmid(level-1,1)
               enddo
             endif
           endif
           ! last block: make ghost cells symmetric (to allow periodicity)
           if(ig1==ng1(level))then
             do ix=ixGhi1-nghostcells+1,ixGextmax1
               ps(igrid)%dx(ix,1)=ps(igrid)%dx(2*(ixGhi1-nghostcells)+1-ix,1)
             enddo
             do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
               psc(igrid)%dx(ix,1)=psc(igrid)%dx(2*(ixCoGmax1-nghostcells)+&
                  1-ix,1)
             enddo
           endif
         endif
       endif
       ! now that dx and grid boundaries are known: fill cell centers
       ifirst=nghostcells+1
       ! first fill the mesh
       summeddx=0.0d0
       do ix=ixMlo1,ixMhi1
         ps(igrid)%x(ix,1)=rnode(rpxmin1_,&
            igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,1)
         summeddx=summeddx+ps(igrid)%dx(ix,1)
       enddo
       ! then ghost cells to left
       summeddx=0.0d0
       do ix=nghostcells,1,-1
         ps(igrid)%x(ix,1)=rnode(rpxmin1_,&
            igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,1)
         summeddx=summeddx+ps(igrid)%dx(ix,1)
       enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi1-nghostcells+1,ixGhi1
         ps(igrid)%x(ix,1)=rnode(rpxmax1_,&
            igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,1)
         summeddx=summeddx+ps(igrid)%dx(ix,1)
       enddo
       ! and next for the coarse representation
       ! first fill the mesh
       summeddx=0.0d0
       do ix=nghostcells+1,ixCoGmax1-nghostcells
         psc(igrid)%x(ix,1)=rnode(rpxmin1_,&
            igrid)+summeddx+0.5d0*psc(igrid)%dx(ix,1)
         summeddx=summeddx+psc(igrid)%dx(ix,1)
       enddo
       ! then ghost cells to left
       summeddx=0.0d0
       do ix=nghostcells,1,-1
         psc(igrid)%x(ix,1)=rnode(rpxmin1_,&
            igrid)-summeddx-0.5d0*psc(igrid)%dx(ix,1)
         summeddx=summeddx+psc(igrid)%dx(ix,1)
       enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
         psc(igrid)%x(ix,1)=rnode(rpxmax1_,&
            igrid)+summeddx+0.5d0*psc(igrid)%dx(ix,1)
         summeddx=summeddx+psc(igrid)%dx(ix,1)
       enddo
       select case(icase)
         case(0)
           ! if even number of ghost cells: xext is just copy of local x
           xext(ixGextmin1:ixGextmax1,1)=ps(igrid)%x(ixGextmin1:ixGextmax1,1)
         case(1)
           ! if uneven number of ghost cells: extra layer left/right
           summeddx=0.0d0
           do ix=ixMlo1,ixMhi1
             xext(ix,1)=rnode(rpxmin1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                1)
            summeddx=summeddx+ps(igrid)%dx(ix,1)
           enddo
           ! then ghost cells to left
           summeddx=0.0d0
           do ix=nghostcells,ixGextmin1,-1
             xext(ix,1)=rnode(rpxmin1_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,&
                1)
             summeddx=summeddx+ps(igrid)%dx(ix,1)
           enddo
          ! then ghost cells to right
          summeddx=0.0d0
          do ix=ixGhi1-nghostcells+1,ixGextmax1
             xext(ix,1)=rnode(rpxmax1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                1)
             summeddx=summeddx+ps(igrid)%dx(ix,1)
          enddo
         case default
           call mpistop("no such case")
       end select
     endif
  endif

  ! calculate area of cell surfaces for standard block
  call get_surface_area(ps(igrid),ixGlo1,ixGhi1)
  ! calculate area of cell surfaces for coarser representative block
  call get_surface_area(psc(igrid),ixCoGmin1,ixCoGmax1)
  ! calculate volume and distance of cells
  ps(igrid)%dsC=1.d0
  select case (coordinate)
    case (Cartesian)
      ps(igrid)%dvolume(ixGextmin1:ixGextmax1)= rnode(rpdx1_,igrid)
      ps(igrid)%ds(ixGextmin1:ixGextmax1,&
         1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,1:ndim)
      ps(igrid)%dsC(ixGextmin1:ixGextmax1,&
         1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,1:ndim)
      psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1)= 2.d0*rnode(rpdx1_,igrid)
      psc(igrid)%ds(ixCoGmin1:ixCoGmax1,1:ndim)=psc(igrid)%dx(&
         ixCoGmin1:ixCoGmax1,1:ndim)
    case (Cartesian_stretched)
      ps(igrid)%dvolume(ixGextmin1:ixGextmax1)= &
         ps(igrid)%dx(ixGextmin1:ixGextmax1,1)
      ps(igrid)%ds(ixGextmin1:ixGextmax1,&
         1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,1:ndim)
      ps(igrid)%dsC(ixGextmin1:ixGextmax1,&
         1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,1:ndim)
      psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1)= &
         psc(igrid)%dx(ixCoGmin1:ixCoGmax1,1)
      psc(igrid)%ds(ixCoGmin1:ixCoGmax1,1:ndim)=psc(igrid)%dx(&
         ixCoGmin1:ixCoGmax1,1:ndim)
    case (spherical)
      ps(igrid)%dvolume(ixGextmin1:ixGextmax1)=(xext(ixGextmin1:ixGextmax1,&
         1)**2 +ps(igrid)%dx(ixGextmin1:ixGextmax1,&
         1)**2/12.0d0)*ps(igrid)%dx(ixGextmin1:ixGextmax1,1)
      psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1)=(psc(igrid)%x(&
         ixCoGmin1:ixCoGmax1,1)**2 +psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
         1)**2/12.0d0)*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,1)
      ps(igrid)%ds(ixGextmin1:ixGextmax1,1)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
         1)
      
      
      ps(igrid)%dsC(ixGextmin1:ixGextmax1,&
         1)=ps(igrid)%dx(ixGextmin1:ixGextmax1,1)
      
      
    case (cylindrical)
      ps(igrid)%dvolume(ixGextmin1:ixGextmax1)=dabs(xext(ixGextmin1:ixGextmax1,&
         1)) *ps(igrid)%dx(ixGextmin1:ixGextmax1,1)
      psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1)=dabs(psc(igrid)%x(&
         ixCoGmin1:ixCoGmax1,1)) *psc(igrid)%dx(ixCoGmin1:ixCoGmax1,1)
      ps(igrid)%ds(ixGextmin1:ixGextmax1,&
         r_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,r_)
      ps(igrid)%dsC(ixGextmin1:ixGextmax1,&
         r_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,r_)
      if(z_>0.and.z_<=ndim) then
        ps(igrid)%ds(ixGextmin1:ixGextmax1,&
           z_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,z_)
        ps(igrid)%dsC(ixGextmin1:ixGextmax1,&
           z_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,z_)
        if(phi_>z_.and.ndir>ndim) then
          ps(igrid)%dsC(ixGextmin1:ixGextmax1,phi_)=xext(ixGextmin1:ixGextmax1,&
             1)+0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,1)
        end if
      end if
      if(phi_>0.and.phi_<=ndim) then
        ps(igrid)%ds(ixGextmin1:ixGextmax1,phi_)=xext(ixGextmin1:ixGextmax1,&
           1)*ps(igrid)%dx(ixGextmin1:ixGextmax1,phi_)
        ps(igrid)%dsC(ixGextmin1:ixGextmax1,phi_)=(xext(ixGextmin1:ixGextmax1,&
           1)+0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           1))*ps(igrid)%dx(ixGextmin1:ixGextmax1,phi_)
        if(z_>phi_.and.ndir>ndim) ps(igrid)%dsC(ixGextmin1:ixGextmax1,z_)=1.d0
      end if
    case default
      call mpistop("Sorry, coordinate unknown")
  end select

  call get_christoffel(ps(igrid),ixGlo1,ixGhi1)

  ! find the blocks on the boundaries
  ps(igrid)%is_physical_boundary=.false.
  
  do i1=-1,1
    if(i1==0) cycle
    ign1=ig1+i1
    ! blocks at periodic boundary have neighbors in the physical domain
    ! thus threated at internal blocks with no physical boundary
    if (periodB(1)) ign1=1+modulo(ign1-1,ng1(level))
    if (ign1 > ng1(level)) then
       if(phi_ > 0 .and. poleB(2,1)) then
         ! if at a pole, the boundary is not physical boundary
         ps(igrid)%is_physical_boundary(2*1)=.false.
       else
         ps(igrid)%is_physical_boundary(2*1)=.true.
       end if
    else if (ign1 < 1) then
       if(phi_ > 0 .and. poleB(1,1)) then
         ! if at a pole, the boundary is not physical boundary
         ps(igrid)%is_physical_boundary(2*1-1)=.false.
       else
         ps(igrid)%is_physical_boundary(2*1-1)=.true.
       end if
    end if
  end do
  
  if(any(ps(igrid)%is_physical_boundary)) then
    phyboundblock(igrid)=.true.
  else
    phyboundblock(igrid)=.false.
  end if

end subroutine alloc_node

!> allocate memory to physical state of igrid node
subroutine alloc_state(igrid, s, ixGmin1,ixGmax1, ixGextmin1,ixGextmax1,&
    alloc_x)
  use mod_global_parameters
  use mod_geometry
  type(state) :: s
  integer, intent(in) :: igrid, ixGmin1,ixGmax1, ixGextmin1,ixGextmax1
  logical, intent(in) :: alloc_x
  integer             :: ixGsmin1,ixGsmax1

  allocate(s%prim(ixGmin1:ixGmax1,1:nprim))
  allocate(s%cons(ixGmin1:ixGmax1,1:ncons))
  s%ixGmin1=ixGmin1;s%ixGmax1=ixGmax1;
   ixGsmin1 = ixGmin1-1; ixGsmax1 = ixGmax1
  if(stagger_grid) then
    allocate(s%prims(ixGsmin1:ixGsmax1,1:nprims))
    s%ixGsmin1=ixGsmin1;s%ixGsmax1=ixGsmax1;
  end if
  if(alloc_x) then
    ! allocate coordinates
    allocate(s%x(ixGmin1:ixGmax1,1:ndim))
    allocate(s%dx(ixGextmin1:ixGextmax1,1:ndim), s%ds(ixGextmin1:ixGextmax1,&
       1:ndim),s%dsC(ixGextmin1:ixGextmax1,1:3))
    allocate(s%dvolume(ixGextmin1:ixGextmax1))
    allocate(s%surfaceC(ixGsmin1:ixGsmax1,1:ndim), s%surface(ixGmin1:ixGmax1,&
       1:ndim))
    allocate(s%christoffel(ixGextmin1:ixGextmax1,1:3,1:3,1:3))
    ! allocate physical boundary flag
    allocate(s%is_physical_boundary(2*ndim))
  else
    ! use spatial info on ps states to save memory
    s%x=>ps(igrid)%x
    s%dx=>ps(igrid)%dx
    s%ds=>ps(igrid)%ds
    s%dsC=>ps(igrid)%dsC
    s%dvolume=>ps(igrid)%dvolume
    s%surfaceC=>ps(igrid)%surfaceC
    s%surface=>ps(igrid)%surface
    s%christoffel=>ps(igrid)%christoffel
    s%is_physical_boundary=>ps(igrid)%is_physical_boundary
  end if
end subroutine alloc_state

subroutine dealloc_state(igrid, s,dealloc_x)
  use mod_global_parameters
  integer, intent(in) :: igrid
  type(state) :: s
  logical, intent(in) :: dealloc_x

  deallocate(s%prim)
  deallocate(s%cons)
  if(stagger_grid) then
    deallocate(s%prims)
  end if
  if(dealloc_x) then
    ! deallocate coordinates
    deallocate(s%x)
    deallocate(s%dx,s%ds)
    deallocate(s%dvolume)
    deallocate(s%surfaceC,s%surface)
    deallocate(s%christoffel)
  else
    nullify(s%x,s%dx,s%ds,s%dvolume,s%surfaceC,s%surface,s%christoffel)
  end if
end subroutine dealloc_state

subroutine dealloc_node(igrid)
  use mod_global_parameters
  use mod_geometry

  integer, intent(in) :: igrid

  if (igrid==0) then
     call mpistop("trying to delete a non-existing grid in dealloc_node")
  end if

  call dealloc_state(igrid, ps(igrid),.true.)
  call dealloc_state(igrid, psc(igrid),.true.)
  call dealloc_state(igrid, ps1(igrid),.false.)
  call dealloc_state(igrid, pso(igrid),.false.)
  ! deallocate temperary solution space
  select case (time_integrator)
  case("ssprk3","ssprk4","jameson","IMEX_Midpoint","IMEX_Trapezoidal")
    call dealloc_state(igrid, ps2(igrid),.false.)
  case("RK3_BT","rk4","ssprk5")
    call dealloc_state(igrid, ps2(igrid),.false.)
    call dealloc_state(igrid, ps3(igrid),.false.)
  case("IMEX_ARS3","IMEX_232")
    call dealloc_state(igrid, ps2(igrid),.false.)
    call dealloc_state(igrid, ps3(igrid),.false.)
    call dealloc_state(igrid, ps4(igrid),.false.)
  end select

end subroutine dealloc_node
