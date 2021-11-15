!> refine one block to its children blocks
subroutine refine_grids(child_igrid,child_ipe,igrid,ipe,active)
  use mod_global_parameters

  integer, dimension(2), intent(in) :: child_igrid, child_ipe
  integer, intent(in) :: igrid, ipe
  logical, intent(in) :: active

  integer :: ic1

  ! allocate solution space for new children
  do ic1=1,2
     call alloc_node(child_igrid(ic1))
  end do

  if ((time_advance .and. active).or.convert.or.firstprocess) then
     ! prolong igrid to new children
     call prolong_grid(child_igrid,child_ipe,igrid,ipe)
  else
     ! Fill new created children with initial condition
     do ic1=1,2
        call initial_condition(child_igrid(ic1))
     end do
  end if

  ! remove solution space of igrid
  !call dealloc_node(igrid)
end subroutine refine_grids

!> prolong one block
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)
  use mod_global_parameters
  use mod_amr_fct, only: old_neighbors

  integer, dimension(2), intent(in) :: child_igrid, child_ipe
  integer, intent(in) :: igrid, ipe

  integer :: ixmin1,ixmax1, ichild, ixComin1,ixComax1, ic1
  double precision :: dxCo1, xComin1, dxFi1, xFimin1

  if (prolongation_method=="linear") then
     dxlevel(1)=rnode(rpdx1_,igrid);
     
     ixmin1=ixMlo1-1;ixmax1=ixMhi1+1;
     


     xComin1=rnode(rpxmin1_,igrid)
     dxCo1=rnode(rpdx1_,igrid)
  end if

  if(stagger_grid) call old_neighbors(child_igrid,child_ipe,igrid,ipe)

  do ic1=1,2
    ichild=child_igrid(ic1)

    ixComin1=ixMlo1+(ic1-1)*block_nx1/2
    ixComax1=ixMhi1+(ic1-2)*block_nx1/2

    if (prolongation_method=="linear") then
       xFimin1=rnode(rpxmin1_,ichild)
       dxFi1=rnode(rpdx1_,ichild)
       call prolong_2nd(ps(igrid),ixComin1,ixComax1,ps(ichild), dxCo1,xComin1,&
          dxFi1,xFimin1,igrid,ichild)
    else
       call prolong_1st(ps(igrid)%prim,ixComin1,ixComax1,ps(ichild)%prim,&
          ps(ichild)%x)
    end if
  end do

end subroutine prolong_grid

!> do 2nd order prolongation
subroutine prolong_2nd(sCo,ixComin1,ixComax1,sFi,dxCo1,xComin1,dxFi1,xFimin1,&
   igridCo,igridFi)
  use mod_physics, only: phys_to_conserved
  use mod_global_parameters
  use mod_amr_fct, only: already_fine, prolong_2nd_stg

  integer, intent(in) :: ixComin1,ixComax1, igridFi, igridCo
  double precision, intent(in) :: dxCo1, xComin1, dxFi1, xFimin1
  type(state), intent(in)      :: sCo
  type(state), intent(inout)   :: sFi

  integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, idim, iw, ixCgmin1,ixCgmax1, el
  double precision :: slopeL, slopeR, slopeC, signC, signR
  double precision :: slope(nprim,ndim)
  double precision :: eta1
  logical :: fine_min1,fine_max1

  associate(wCo=>sCo%prim, wFi=>sFi%prim)
  ixCgmin1=ixComin1;ixCgmax1=ixComax1;
  
  do ixCo1 = ixCgmin1,ixCgmax1
     ! lower left grid index in finer child block
     ixFi1=2*(ixCo1-ixComin1)+ixMlo1

     do idim=1,ndim
        hxCo1=ixCo1-kr(1,idim)
        jxCo1=ixCo1+kr(1,idim)

        do iw=1,nprim
           slopeL=wCo(ixCo1,iw)-wCo(hxCo1,iw)
           slopeR=wCo(jxCo1,iw)-wCo(ixCo1,iw)
           slopeC=half*(slopeR+slopeL)

           ! get limited slope
           signR=sign(one,slopeR)
           signC=sign(one,slopeC)
           select case(typeprolonglimit)
           case('unlimit')
             slope(iw,idim)=slopeC
           case('minmod')
             slope(iw,idim)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
           case('woodward')
             slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), signR*slopeL,&
                signR*half*slopeC))
           case('koren')
             slope(iw,idim)=signR*max(zero,min(two*signR*slopeL,&
                 (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
           case default
             slope(iw,idim)=signC*max(zero,min(dabs(slopeC), signC*slopeL,&
                signC*slopeR))
           end select
        end do
     end do
     ! cell-centered coordinates of coarse grid point
     !^D&xCo^D=xCo({ixCo^DD},^D)\
     do ix1=ixFi1,ixFi1+1 
        ! cell-centered coordinates of fine grid point
        !^D&xFi^D=xFi({ix^DD},^D)\
        if(slab_uniform) then
          ! normalized distance between fine/coarse cell center
          ! in coarse cell: ranges from -0.5 to 0.5 in each direction
          ! (origin is coarse cell center)
          ! hence this is +1/4 or -1/4 on cartesian mesh
          !eta^D=(xFi^D-xCo^D)*invdxCo^D;
          eta1=0.5d0*(dble(ix1-ixFi1)-0.5d0);
        else
          ! forefactor is -0.5d0 when ix=ixFi and +0.5d0 for ixFi+1
          eta1=(dble(ix1-ixFi1)-0.5d0)*(one-sFi%dvolume(ix1) &
             /sum(sFi%dvolume(ixFi1:ixFi1+1)))  
        end if
        wFi(ix1,1:nprim) = wCo(ixCo1,1:nprim) + (slope(1:nprim,1)*eta1)
     end do
  end do
  if(stagger_grid) then
    call already_fine(sFi,igridFi,fine_min1,fine_max1)
    call prolong_2nd_stg(sCo,sFi,ixComin1,ixComax1,ixMlo1,ixMhi1,dxCo1,xComin1,&
       dxFi1,xFimin1,.false.,fine_min1,fine_max1)
  end if

  end associate

end subroutine prolong_2nd

!> do 1st order prolongation
subroutine prolong_1st(wCo,ixComin1,ixComax1,wFi,xFi)
  use mod_global_parameters

  integer, intent(in) :: ixComin1,ixComax1
  double precision, intent(in) :: wCo(ixGlo1:ixGhi1,nprim), xFi(ixGlo1:ixGhi1,&
     1:ndim)
  double precision, intent(out) :: wFi(ixGlo1:ixGhi1,nprim)

  integer :: ixCo1, ixFi1, iw
  integer :: ixFimin1,ixFimax1

  do ixCo1 = ixComin1,ixComax1
     ixFi1=2*(ixCo1-ixComin1)+ixMlo1
     forall(iw=1:nprim) wFi(ixFi1:ixFi1+1,iw)=wCo(ixCo1,iw)
  end do

end subroutine prolong_1st
