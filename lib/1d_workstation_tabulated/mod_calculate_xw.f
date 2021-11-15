!> Handles computations for coordinates and variables in output
module mod_calculate_xw
  implicit none
 
contains 

  !> Compute both corner as well as cell-centered values for output
  subroutine calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
     normconv,ixCmin1,ixCmax1,ixCCmin1,ixCCmax1,first)
    ! this subroutine computes both corner as well as cell-centered values
    ! it handles how we do the center to corner averaging, as well as 
    ! whether we switch to cartesian or want primitive or conservative output,
    ! handling the addition of B0 in B0+B1 cases, ...
    !
    ! the normconv is passed on to usr_aux_output for extending with
    ! possible normalization values for the nw+1:nw+nwauxio entries

    use mod_usr_methods, only: usr_aux_output
    use mod_global_parameters
    use mod_limiter
    use mod_physics, only: physics_type, phys_to_primitive

    integer, intent(in) :: qunit, igrid
    double precision, intent(in), dimension(ixMlo1-1:ixMhi1,ndim) :: xC
    double precision, intent(in), dimension(ixMlo1:ixMhi1,ndim)   :: xCC
    integer :: ixCmin1,ixCmax1,ixCCmin1,ixCCmax1
    logical, intent(in) :: first

    double precision, dimension(ixMlo1-1:ixMhi1,ndim) :: xC_TMP
    double precision, dimension(ixMlo1:ixMhi1,ndim)   :: xCC_TMP
    double precision, dimension(ixMlo1-1:ixMhi1,nprim+nwauxio)   :: wC_TMP
    double precision, dimension(ixMlo1:ixMhi1,nprim+nwauxio)     :: wCC_TMP
    double precision,dimension(0:nprim+nwauxio),intent(out)       :: normconv 

    double precision :: ldw(ixGlo1:ixGhi1), dwC(ixGlo1:ixGhi1)
    double precision, dimension(ixMlo1-1:ixMhi1,nprim+nwauxio)   :: wC
    double precision, dimension(ixMlo1:ixMhi1,nprim+nwauxio)     :: wCC
    double precision, dimension(ixGlo1:ixGhi1,1:nprim+nwauxio)   :: w
    integer :: nxCC1,idims,jxCmin1,jxCmax1,iwe
    integer :: nx1, nxC1, ix1, ix, iw, level, idir
    logical, save :: subfirst=.true.

    ixCmin1=ixMlo1-1; ixCmax1=ixMhi1; ! Corner indices
    ixCCmin1=ixMlo1; ixCCmax1=ixMhi1; ! Center indices

    saveigrid=igrid
    nx1=ixMhi1-ixMlo1+1;
    level=node(plevel_,igrid)

    normconv(0) = length_convert_factor
    normconv(1:nprim) = w_convert_factor
    block=>ps(igrid)
    w(ixGlo1:ixGhi1,1:nprim)=ps(igrid)%prim(ixGlo1:ixGhi1,1:nprim)

!    if (nwextra>0) then
!     ! here we actually fill the ghost layers for the nwextra variables using 
!     ! continuous extrapolation (as these values do not exist normally in ghost cells)
!     do idims=1,ndim
!      select case(idims)
!       {case(^D)
!         jxCmin^DD=ixGhi^D+1-nghostcells^D%jxCmin^DD=ixGlo^DD;
!         jxCmax^DD=ixGhi^DD;
!         do ix^D=jxCmin^D,jxCmax^D
!             w(ix^D^D%jxC^S,nprim-nwextra+1:nprim) = w(jxCmin^D-1^D%jxC^S,nprim-nwextra+1:nprim)
!         end do 
!         jxCmin^DD=ixGlo^DD;
!         jxCmax^DD=ixGlo^D-1+nghostcells^D%jxCmax^DD=ixGhi^DD;
!         do ix^D=jxCmin^D,jxCmax^D
!             w(ix^D^D%jxC^S,nprim-nwextra+1:nprim) = w(jxCmax^D+1^D%jxC^S,nprim-nwextra+1:nprim)
!         end do \}
!      end select
!     end do
!    end if

    ! next lines needed when usr_aux_output uses gradients
    ! and later on when dwlimiter2 is used 
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
    dxlevel(1)=rnode(rpdx1_,igrid);
    block=>ps(igrid)
    if(nwauxio>0)then
      ! auxiliary io variables can be computed and added by user
      ! next few lines ensure correct usage of routines like divvector etc
      ! default (no) normalization for auxiliary variables
      normconv(nprim+1:nprim+nwauxio)=one

      if (.not. associated(usr_aux_output)) then
         call mpistop("usr_aux_output not defined")
      else
         call usr_aux_output(ixGlo1,ixGhi1,ixMlo1-1,ixMhi1+1,w,ps(igrid)%x,&
            normconv)
      end if
    endif

    ! In case primitives to be saved: use primitive subroutine
    !  extra layer around mesh only needed when storing corner values and averaging

    ! compute the cell-center values for w first
    ! cell center values obtained from mere copy
    wCC(ixCCmin1:ixCCmax1,:)=w(ixCCmin1:ixCCmax1,:)

    ! compute the corner values for w now by averaging

    if(slab) then
       ! for slab symmetry: no geometrical info required
       do iw=1,nprim+nwauxio
         do ix1=ixCmin1,ixCmax1
            wC(ix1,iw)=sum(w(ix1:ix1+1,iw))/dble(2**ndim)
         end do
       end do
    else
       do iw=1,nprim+nwauxio
         do ix1=ixCmin1,ixCmax1
           wC(ix1,iw)=sum(w(ix1:ix1+1,iw)*ps(igrid)%dvolume(ix1:ix1+1)) &
              /sum(ps(igrid)%dvolume(ix1:ix1+1))
         end do
       end do
    endif

    if(nocartesian) then
      ! keep the coordinate and vector components
      xC_TMP(ixCmin1:ixCmax1,1:ndim)          = xC(ixCmin1:ixCmax1,1:ndim)
      wC_TMP(ixCmin1:ixCmax1,1:nprim+nwauxio)    = wC(ixCmin1:ixCmax1,&
         1:nprim+nwauxio)
      xCC_TMP(ixCCmin1:ixCCmax1,1:ndim)        = xCC(ixCCmin1:ixCCmax1,1:ndim)
      wCC_TMP(ixCCmin1:ixCCmax1,1:nprim+nwauxio)  = wCC(ixCCmin1:ixCCmax1,&
         1:nprim+nwauxio)
    else
      ! do all conversions to cartesian coordinates and vector components
      ! start for the corner values
      call to_cartesian(xC_TMP,wC_TMP,ixCmin1,ixCmax1,xC,wC)
      ! then cell center values
      call to_cartesian(xCC_TMP,wCC_TMP,ixCCmin1,ixCCmax1,xCC,wCC)
    endif

    ! Warning: differentiate between idl/idlCC/tecplot/tecplotCC/vtu(B)/vtu(B)CC
    !if(nwaux>0 .and. mype==0 .and. first.and.subfirst) then
      ! when corner values are computed and auxiliaries present: warn!
    !  if(convert_type=='idl'.or.convert_type=='tecplot' &
    !   .or.convert_type=='vtu'.or.convert_type=='vtuB') &
    !      write(*,*) 'Warning: also averaged auxiliaries within calc_grid'
    !  subfirst=.false.
    !endif

  end subroutine calc_grid

  !> convert to cartesian coordinates and vector components
  subroutine to_cartesian(x_TMP,w_TMP,ixmin1,ixmax1,xC,wC)
    ! conversion of coordinate and vector components from cylindrical/spherical
    ! to cartesian coordinates and components done here
    ! Also: nullifying values lower than smalldouble
    use mod_global_parameters
    use mod_geometry

    integer :: ixmin1,ixmax1, ix1, idim, iw, ivector, iw0
    integer, dimension(nprim) :: vectoriw
    double precision :: x_TEC(ndim), w_TEC(nprim+nwauxio)
    double precision, dimension(ndim,ndim) :: normal

    double precision, dimension(ixmin1:ixmax1,ndim) :: xC
    double precision, dimension(ixmin1:ixmax1,nprim+nwauxio)   :: wC

    double precision, dimension(ixmin1:ixmax1,ndim) :: x_TMP
    double precision, dimension(ixmin1:ixmax1,nprim+nwauxio)   :: w_TMP

    iw0=0
    vectoriw=-1
    if(nvector>0) then
      do ivector=1,nvector
         do idim=1,ndim
            vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
         end do
      end do
    endif
    x_TEC=0.d0
    do ix1=ixmin1,ixmax1
       select case (coordinate)
       case (Cartesian,Cartesian_stretched)
          x_TEC(1:ndim)=xC(ix1,1:ndim)
          w_TEC(1:nprim+nwauxio)=wC(ix1,1:nprim+nwauxio)
       case (cylindrical)
          
          x_TEC(1)=xC(ix1,1)
          
          

       case (spherical)
          x_TEC(1)=xC(ix1,1)
          
          
       case default
          write(*,*) "No converter for coordinate=",coordinate
       end select
       ! there is no need to convert the vectors to other coordinates
       w_TEC(1:nprim+nwauxio)=wC(ix1,1:nprim+nwauxio)
       x_TMP(ix1,1:ndim)=x_TEC(1:ndim)
       w_TMP(ix1,1:nprim+nwauxio)=w_TEC(1:nprim+nwauxio)
       ! Be aware that small values are nullified here!!!
       where(dabs(w_TMP(ix1,1:nprim+nwauxio))<smalldouble)
             w_TMP(ix1,1:nprim+nwauxio)=zero
       endwhere
    end do

  end subroutine to_cartesian

  !> get all variables names
  subroutine getheadernames(wnamei,xandwnamei,outfilehead)
    ! this collects all variables names in the wnamei character array, getting the info from
    ! the prim_wnames/cons_wnames strings (depending on saveprim). It combines this info with names
    ! for the dimensional directions in the xandwnamei array. In the outfilehead, it collects
    ! the dimensional names, and only those names from the nprim variables for output (through w_write)
    ! together with the added names for nwauxio variables
    use mod_global_parameters
    use mod_geometry

    character(len=name_len)   :: wnamei(1:nprim+nwauxio),&
       xandwnamei(1:ndim+nprim+nwauxio)
    character(len=1024) :: outfilehead

    integer::  space_position,iw
    character(len=name_len)::  wname
    character(len=std_len):: aux_variable_names
    character(len=std_len)::  scanstring

    logical, save:: first=.true.

    ! --- part to provide variable names from prim_wnames/cons_wnames strings
    scanstring=TRIM(aux_variable_names)
    wnamei(1:nprim)=prim_names(1:nprim)

    space_position=index(scanstring,' ')
    do iw=nprim+1,nprim+nwauxio
       do while (space_position==1)
         scanstring=scanstring(2:)
         space_position=index(scanstring,' ')
       enddo
       wname=scanstring(:space_position-1)
       scanstring=scanstring(space_position+1:)
       space_position=index(scanstring,' ')

       ! fill all names, even those that we will not write away from the first nprim
       wnamei(iw)=TRIM(wname)
    enddo
    ! --- end of part to provide variable names 

    select case (coordinate)
       case( spherical )
          xandwnamei(1)="r";;
       case( cylindrical )
          xandwnamei(1)="R";
          
          
       case default
          xandwnamei(1)="X";;
    end select

    xandwnamei(ndim+1:ndim+nprim+nwauxio)=wnamei(1:nprim+nwauxio)

    ! in outfilehead, collect the dimensional names, and all output variable names
    ! first all dimensions
    write(outfilehead,'(a)') TRIM(xandwnamei(1))
    
    ! then all nw variables, with w_write control for inclusion
    do iw=ndim+1,ndim+nprim
       wname=xandwnamei(iw)
       if(w_write(iw-ndim)) then
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(&
       wname)
       endif
    enddo
    ! then all nwauxio variables
    if(nwauxio>0) then
      do iw=ndim+nprim+1,ndim+nprim+nwauxio
         wname=xandwnamei(iw)
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(&
       wname)
      enddo
    endif

    if(first.and.mype==0)then
      print*,'-------------------------------------------------------------------------------'
      write(unitterm,*&
         )'Saving visual data. Coordinate directions and variable names are:'
      do iw=1,ndim
        print *,iw,xandwnamei(iw)
      enddo
      do iw=ndim+1,ndim+nprim+nwauxio
        print *,iw,wnamei(iw-ndim),xandwnamei(iw)
      enddo
      write(unitterm,*)'time =', global_time
      print*,'-------------------------------------------------------------------------------'
      first=.false.
    endif

  end subroutine getheadernames

  !> computes cell corner (xC) and cell center (xCC) coordinates
  subroutine calc_x(igrid,xC,xCC)
    use mod_global_parameters

    integer, intent(in)               :: igrid
    double precision, intent(out)     :: xC(ixMlo1-1:ixMhi1,ndim)
    double precision, intent(out)     :: xCC(ixMlo1:ixMhi1,ndim)

    integer                           :: ixCmin1,ixCmax1, idims, level, ix

    level=node(plevel_,igrid)

    ! coordinates of cell centers
    xCC(ixMlo1:ixMhi1,:)=ps(igrid)%x(ixMlo1:ixMhi1,:)

    ! coordinates of cell corners
    ixCmin1=ixMlo1-1; ixCmax1=ixMhi1;
    if(slab)then
       do idims=1,ndim
         xC(ixCmin1:ixCmax1,idims)=ps(igrid)%x(ixCmin1:ixCmax1,&
            idims)+0.5d0*dx(idims,level)
       end do
    else
       ! for any non-cartesian or stretched coordinate (allow multiple stretched directions)
       do ix=ixCmin1,ixCmax1
         xC(ix,1)=ps(igrid)%x(ix,1)+0.5d0*ps(igrid)%dx(ix,1)
       end do
    endif

  end subroutine calc_x

end module mod_calculate_xw
