!> coarsen one grid to its coarser representative
subroutine coarsen_grid(sFi,ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2,ixFimin1,&
   ixFimin2,ixFimax1,ixFimax2,sCo,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
   ixComin1,ixComin2,ixComax1,ixComax2)
  use mod_global_parameters
  use mod_physics

  type(state), intent(inout)      :: sFi, sCo
  integer, intent(in) :: ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2, ixFimin1,&
     ixFimin2,ixFimax1,ixFimax2, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
      ixComin1,ixComin2,ixComax1,ixComax2

  integer :: ixCo1,ixCo2, ixFi1,ixFi2, iw
  double precision :: CoFiratio

  associate(wFi=>sFi%prim(ixFiGmin1:ixFiGmax1,ixFiGmin2:ixFiGmax2,1:nprim),&
      wCo=>sCo%prim(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nprim))
  staggered: associate(wFis=>sFi%prims,wCos=>sCo%prims)
  ! coarsen by 2 in every direction - conservatively

  if(slab_uniform) then
    CoFiratio=one/dble(2**ndim)
    do iw=1,nprim
       do ixCo2 = ixComin2,ixComax2
          ixFi2=2*(ixCo2-ixComin2)+ixFimin2
       do ixCo1 = ixComin1,ixComax1
          ixFi1=2*(ixCo1-ixComin1)+ixFimin1
          wCo(ixCo1,ixCo2,iw)=sum(wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
             iw))*CoFiratio
       end do
       end do
    end do
  else
    do iw=1,nprim
      do ixCo2 = ixComin2,ixComax2
         ixFi2=2*(ixCo2-ixComin2)+ixFimin2
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,ixCo2,iw)= sum(sFi%dvolume(ixFi1:ixFi1+1,&
            ixFi2:ixFi2+1)*wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
            iw)) /sCo%dvolume(ixCo1,ixCo2)
      end do
      end do
    end do
  end if

  if(stagger_grid) then
    do iw=1,nprims
      ! Start one layer before
      do ixCo2 = ixComin2-kr(2,iw),ixComax2
         ixFi2=2*(ixCo2-ixComin2+kr(2,iw))+ixFimin2-kr(2,iw)
      do ixCo1 = ixComin1-kr(1,iw),ixComax1
         ixFi1=2*(ixCo1-ixComin1+kr(1,iw))+ixFimin1-kr(1,iw)
         ! This if statement catches the axis where surface is zero:
         if (sCo%surfaceC(ixCo1,ixCo2,iw)>1.0d-9*sCo%dvolume(ixCo1,&
            ixCo2)) then !Normal case
           wCos(ixCo1,ixCo2,iw)=sum(sFi%surfaceC(ixFi1:ixFi1+1-kr(iw,1),&
              ixFi2:ixFi2+1-kr(iw,2),iw)*wFis(ixFi1:ixFi1+1-kr(iw,1),&
              ixFi2:ixFi2+1-kr(iw,2),iw)) /sCo%surfaceC(ixCo1,ixCo2,iw)
         else ! On axis
           wCos(ixCo1,ixCo2,iw)=zero
         end if
      end do
      end do
    end do
    ! average to fill cell-centred values
    call phys_face_to_center(ixComin1,ixComin2,ixComax1,ixComax2,sCo)
  end if

  end associate staggered
  end associate
end subroutine coarsen_grid
