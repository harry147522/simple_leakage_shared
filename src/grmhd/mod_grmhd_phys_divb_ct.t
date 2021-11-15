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

  subroutine grmhd_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixI^S,1:nprim)
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine grmhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters
    !use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer                            :: hxC^L,ixC^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixI^S,1:ndim)

    associate(bfaces=>s%prims,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    fE=zero

    do idim1=1,ndim 
      iwdim1 = Bvec(idim1)
      do idim2=1,ndim
        iwdim2 = Bvec(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            hxC^L=ixC^L+kr(idim2,^D);
            ! Interpolate to edges
            fE(ixC^S,idir)=quarter*(fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)&
                                   -fC(ixC^S,iwdim2,idim1)-fC(hxC^S,iwdim2,idim1))

            fE(ixC^S,idir)=qdt*s%dsC(ixC^S,idir)*fE(ixC^S,idir)

            if (.not.slab) then
              where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixC^S,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    !if(associated(usr_set_electric_field)) &
    !  call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate
  end subroutine update_faces_average

  !> calculate cell-center values from face-center values
  subroutine mhd_face_to_center(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim

    associate(prim=>s%prim, prims=>s%prims)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the prim arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxO^L=ixO^L-kr(idim,^D);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      prim(ixO^S,Bvec(idim))=half/s%surface(ixO^S,idim)*&
        (prims(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
        +prims(hxO^S,idim)*s%surfaceC(hxO^S,idim))
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
