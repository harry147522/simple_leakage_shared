module mod_grmhd_phys_divb_ct
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_divb_ct_init

contains

  !> Initialize the module
  subroutine grmhd_phys_divb_ct_init()
    use mod_global_parameters
    stagger_grid = .true.
    nprims=ndim
    !phys_update_faces   => grmhd_update_faces
    !phys_face_to_center => grmhd_face_to_center
    !phys_modify_wLR     => grmhd_modify_wLR
  end subroutine grmhd_phys_divb_ct_init

  subroutine grmhd_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,s)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nprim)
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,qt,qdt,fC,fE,sCT,s)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine grmhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters
    !use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)

    integer                            :: hxCmin1,hxCmin2,hxCmax1,hxCmax2,&
       ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
       ixCmmin1,ixCmmin2,ixCmmax1,ixCmmax2
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    associate(bfaces=>s%prims,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;
    ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;

    fE=zero

    do idim1=1,ndim 
      iwdim1 = Bvec(idim1)
      do idim2=1,ndim
        iwdim2 = Bvec(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmax1=ixCmax1+kr(idim1,1);jxCmax2=ixCmax2+kr(idim1,2);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmax1=ixCmax1+kr(idim2,1);hxCmax2=ixCmax2+kr(idim2,2);
            ! Interpolate to edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwdim1,&
               idim2)+fC(jxCmin1:jxCmax1,jxCmin2:jxCmax2,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwdim2,&
               idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iwdim2,idim1))

            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=qdt*s%dsC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)

            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    !if(associated(usr_set_electric_field)) &
    !  call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmax1=ixCmax1-kr(idim2,1);hxCmax2=ixCmax2-kr(idim2,2);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2);
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=zero
      end where
      ! Time update
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=bfaces(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)
    end do

    end associate
  end subroutine update_faces_average

  !> calculate cell-center values from face-center values
  subroutine mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state)                        :: s

    integer                            :: fxOmin1,fxOmin2,fxOmax1,fxOmax2,&
        gxOmin1,gxOmin2,gxOmax1,gxOmax2, hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
        jxOmin1,jxOmin2,jxOmax1,jxOmax2, kxOmin1,kxOmin2,kxOmax1,kxOmax2, idim

    associate(prim=>s%prim, prims=>s%prims)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the prim arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
      hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      prim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         Bvec(idim))=half/s%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idim)*(prims(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idim)*s%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idim)+prims(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         idim)*s%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,idim))
    end do

    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,Bvec(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,Bvec(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do

    end associate

  end subroutine mhd_face_to_center

end module mod_grmhd_phys_divb_ct
