!> fill ghost cells at a physical boundary
subroutine bc_phys(iside,idims,time,qdt,s,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
   ixBmin1,ixBmin2,ixBmax1,ixBmax2)
  use mod_usr_methods, only: usr_special_bc
  use mod_bc_data, only: bc_data_set
  use mod_global_parameters
  use mod_physics

  integer, intent(in) :: iside, idims, ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2
  double precision, intent(in) :: time,qdt
  type(state), intent(inout) :: s
  double precision :: wtmp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nprim)

  integer :: idir, is
  integer :: ixOsmin1,ixOsmin2,ixOsmax1,ixOsmax2,hxOmin1,hxOmin2,hxOmax1,&
     hxOmax2,jxOmin1,jxOmin2,jxOmax1,jxOmax2
  double precision :: Q(ixGmin1:ixGmax1,ixGmin2:ixGmax2),Qp(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2) 
  integer :: iw, iB, ix1,ix2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixMmin1,ixMmin2,&
     ixMmax1,ixMmax2, nghostcellsi,iib1,iib2
  logical  :: isphysbound
  integer :: vel_index

  associate(x=>s%x,w=>s%prim,ws=>s%prims)
  select case (idims)
  case (1)
     if (iside==2) then
        ! maximal boundary
        iB=ismax1
        ixOmin1=ixBmax1+1-nghostcells;ixOmin2=ixBmin2;
        ixOmax1=ixBmax1;ixOmax2=ixBmax2;
        ! cont/symm/asymm types
        do iw=1,nprim
           select case (typeboundary(iw,iB))
           case ("symm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) = w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,iw)
           case ("asymm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) =-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,iw)
           case ("cont")
              do ix1=ixOmin1,ixOmax1
                 w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmin1-1,ixOmin2:ixOmax2,iw)
              end do
           case("noinflow")
 !noinflow applies only on velocities, for other vectors, use cont
 !Note that although 1 = ndim is not ndir, in case of ndir > ndim, 
              ! the lastest velocity has no corresponding boundary
              ! so apply cont should be enough.
              if ( allocated(veloc) ) then
                 vel_index = veloc(1)
              else if ( allocated(W_vel) ) then
                 vel_index = W_vel(1)
              else
                 call mpistop("something got wrong, no veloc is allocated")
              end if

              if ( iw==vel_index )then
                do ix1=ixOmin1,ixOmax1
                    w(ix1,ixOmin2:ixOmax2,iw) = max(w(ixOmin1-1,&
                       ixOmin2:ixOmax2,iw),zero)
                end do
              else
                do ix1=ixOmin1,ixOmax1
                    w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmin1-1,ixOmin2:ixOmax2,&
                       iw)
                end do
              end if
           case ("special", "bc_data")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("character")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case ("periodic")
  !            call mpistop("periodic bc info should come from neighbors")
           case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB),&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
     else
        ! minimal boundary
        iB=ismin1
        ixOmin1=ixBmin1;ixOmin2=ixBmin2;
        ixOmax1=ixBmin1-1+nghostcells;ixOmax2=ixBmax2;
        ! cont/symm/asymm types
        do iw=1,nprim
           select case (typeboundary(iw,iB))
           case ("symm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) = w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,iw)
           case ("asymm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) =-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,iw)
           case ("cont")
              do ix1=ixOmin1,ixOmax1
                 w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmax1+1,ixOmin2:ixOmax2,iw)
              end do
           case("noinflow")
              if ( allocated(veloc) ) then
                 vel_index = veloc(1)
              else if ( allocated(W_vel) ) then
                 vel_index = W_vel(1)
              else
                 call mpistop("something got wrong, no veloc is allocated")
              end if

              if ( iw==vel_index )then
                 do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,iw) = min(w(ixOmax1+1,ixOmin2:ixOmax2,&
                      iw),zero)
                 end do
              else
                 do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmax1+1,ixOmin2:ixOmax2,iw)
                 end do
              end if
           case ("special", "bc_data")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("character")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case ("periodic")
  !            call mpistop("periodic bc info should come from neighbors")
           case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB),&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
     end if 
  case (2)
     if (iside==2) then
        ! maximal boundary
        iB=ismax2
        ixOmin1=ixBmin1;ixOmin2=ixBmax2+1-nghostcells;
        ixOmax1=ixBmax1;ixOmax2=ixBmax2;
        ! cont/symm/asymm types
        do iw=1,nprim
           select case (typeboundary(iw,iB))
           case ("symm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = w(ixOmin1:ixOmax1,&
                 ixOmin2-1:ixOmin2-nghostcells:-1,iw)
           case ("asymm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =-w(ixOmin1:ixOmax1,&
                 ixOmin2-1:ixOmin2-nghostcells:-1,iw)
           case ("cont")
              do ix2=ixOmin2,ixOmax2
                 w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmin2-1,iw)
              end do
           case("noinflow")
 !noinflow applies only on velocities, for other vectors, use cont
 !Note that although 2 = ndim is not ndir, in case of ndir > ndim, 
              ! the lastest velocity has no corresponding boundary
              ! so apply cont should be enough.
              if ( allocated(veloc) ) then
                 vel_index = veloc(2)
              else if ( allocated(W_vel) ) then
                 vel_index = W_vel(2)
              else
                 call mpistop("something got wrong, no veloc is allocated")
              end if

              if ( iw==vel_index )then
                do ix2=ixOmin2,ixOmax2
                    w(ixOmin1:ixOmax1,ix2,iw) = max(w(ixOmin1:ixOmax1,&
                       ixOmin2-1,iw),zero)
                end do
              else
                do ix2=ixOmin2,ixOmax2
                    w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmin2-1,&
                       iw)
                end do
              end if
           case ("special", "bc_data")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("character")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case ("periodic")
  !            call mpistop("periodic bc info should come from neighbors")
           case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB),&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
     else
        ! minimal boundary
        iB=ismin2
        ixOmin1=ixBmin1;ixOmin2=ixBmin2;
        ixOmax1=ixBmax1;ixOmax2=ixBmin2-1+nghostcells;
        ! cont/symm/asymm types
        do iw=1,nprim
           select case (typeboundary(iw,iB))
           case ("symm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = w(ixOmin1:ixOmax1,&
                 ixOmax2+nghostcells:ixOmax2+1:-1,iw)
           case ("asymm")
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =-w(ixOmin1:ixOmax1,&
                 ixOmax2+nghostcells:ixOmax2+1:-1,iw)
           case ("cont")
              do ix2=ixOmin2,ixOmax2
                 w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmax2+1,iw)
              end do
           case("noinflow")
              if ( allocated(veloc) ) then
                 vel_index = veloc(2)
              else if ( allocated(W_vel) ) then
                 vel_index = W_vel(2)
              else
                 call mpistop("something got wrong, no veloc is allocated")
              end if

              if ( iw==vel_index )then
                 do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,iw) = min(w(ixOmin1:ixOmax1,ixOmax2+1,&
                      iw),zero)
                 end do
              else
                 do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmax2+1,iw)
                 end do
              end if
           case ("special", "bc_data")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("character")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case ("periodic")
  !            call mpistop("periodic bc info should come from neighbors")
           case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB),&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
     end if 
  end select

  ! do user defined special boundary conditions
  if (any(typeboundary(1:nprim,iB)=="special")) then
     if (.not. associated(usr_special_bc)) call &
        mpistop("usr_special_bc not defined")
     call usr_special_bc(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,iB,w,x)
  end if

  ! fill boundary conditions from external data vtk files
  if (any(typeboundary(1:nprim,iB)=="bc_data")) then
     call bc_data_set(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,iB,w,x)
  end if

  
  !end do

  end associate
end subroutine bc_phys

!> fill inner boundary values
subroutine getintbc(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2)
  use mod_usr_methods, only: usr_internal_bc
  use mod_global_parameters

  double precision, intent(in)   :: time
  integer, intent(in)            :: ixGmin1,ixGmin2,ixGmax1,ixGmax2

  integer :: iigrid, igrid, ixOmin1,ixOmin2,ixOmax1,ixOmax2,level

  ixOmin1=ixGmin1+nghostcells;ixOmin2=ixGmin2+nghostcells
  ixOmax1=ixGmax1-nghostcells;ixOmax2=ixGmax2-nghostcells;

  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
  !do iigrid=1,igridstail; igrid=igrids(iigrid);
     dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
     block=>ps(igrid)
     typelimiter=type_limiter(node(plevel_,igrid))
     typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
     level=node(plevel_,igrid)
     saveigrid=igrid

     if (associated(usr_internal_bc)) then
        call usr_internal_bc(level,time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2,ps(igrid)%prim,ps(igrid)%x)
     end if
  end do

end subroutine getintbc
