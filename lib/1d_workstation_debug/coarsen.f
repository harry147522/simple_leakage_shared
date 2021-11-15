!> coarsen one grid to its coarser representative
subroutine coarsen_grid(sFi,ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,sCo,&
   ixCoGmin1,ixCoGmax1,ixComin1,ixComax1)
  use mod_global_parameters
  use mod_physics

  type(state), intent(inout)      :: sFi, sCo
  integer, intent(in) :: ixFiGmin1,ixFiGmax1, ixFimin1,ixFimax1, ixCoGmin1,&
     ixCoGmax1, ixComin1,ixComax1

  integer :: ixCo1, ixFi1, iw
  double precision :: CoFiratio

  associate(wFi=>sFi%prim(ixFiGmin1:ixFiGmax1,1:nprim),&
      wCo=>sCo%prim(ixCoGmin1:ixCoGmax1,1:nprim))
  staggered: associate(wFis=>sFi%prims,wCos=>sCo%prims)
  ! coarsen by 2 in every direction - conservatively

  if(slab_uniform) then
    CoFiratio=one/dble(2**ndim)
    do iw=1,nprim
       do ixCo1 = ixComin1,ixComax1
          ixFi1=2*(ixCo1-ixComin1)+ixFimin1
          wCo(ixCo1,iw)=sum(wFi(ixFi1:ixFi1+1,iw))*CoFiratio
       end do
    end do
  else
    do iw=1,nprim
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,iw)= sum(sFi%dvolume(ixFi1:ixFi1+1)*wFi(ixFi1:ixFi1+1,&
            iw)) /sCo%dvolume(ixCo1)
      end do
    end do
  end if

  if(stagger_grid) then
    do iw=1,nprims
      ! Start one layer before
      do ixCo1 = ixComin1-kr(1,iw),ixComax1
         ixFi1=2*(ixCo1-ixComin1+kr(1,iw))+ixFimin1-kr(1,iw)
         ! This if statement catches the axis where surface is zero:
         if (sCo%surfaceC(ixCo1,iw)>1.0d-9*sCo%dvolume(ixCo1)) then !Normal case
           wCos(ixCo1,iw)=sum(sFi%surfaceC(ixFi1:ixFi1+1-kr(iw,1),&
              iw)*wFis(ixFi1:ixFi1+1-kr(iw,1),iw)) /sCo%surfaceC(ixCo1,iw)
         else ! On axis
           wCos(ixCo1,iw)=zero
         end if
      end do
    end do
    ! average to fill cell-centred values
    call phys_face_to_center(ixComin1,ixComax1,sCo)
  end if

  end associate staggered
  end associate
end subroutine coarsen_grid
