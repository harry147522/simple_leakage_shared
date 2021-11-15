!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry
  implicit none
  public

  integer :: coordinate=-1
  integer, parameter :: Cartesian          = 0
  integer, parameter :: Cartesian_stretched= 1
  integer, parameter :: cylindrical        = 2
  integer, parameter :: spherical          = 3


contains

  !> Set the coordinate system to be used
  subroutine set_coordinate_system(geom)
    use mod_global_parameters
    character(len=*), intent(in) :: geom !< Name of the coordinate system

    ! Store the geometry name
    geometry_name = geom

    select case (geom)
    case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
      ndir = ndim
      coordinate=Cartesian
    case ("Cartesian_1.5D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
      coordinate=Cartesian
      ndir = 2
    case ("Cartesian_1.75D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
      coordinate=Cartesian
      ndir = 3
    case ("Cartesian_2.5D")
      if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
      coordinate=Cartesian
      ndir = 3
    case ("cylindrical","cylindrical_2D","cylindrical_3D")
      ndir = ndim
      r_   = 1
      z_   = 2
      if(ndir==3) phi_ = 3
      coordinate=cylindrical
    case ("cylindrical_2.5D")
      if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      z_   = 2
      phi_ = 3
      coordinate=cylindrical
    case ("polar","polar_2D","polar_3D")
      ndir = ndim
      r_   = 1
      phi_ = 2
      if(ndir==3) z_ = 3
      coordinate=cylindrical
    case ("polar_1.5D")
       if (ndim /= 1) call mpistop("Geometry polar_1.5D but ndim /= 1")
       ndir = 2
       r_   = 1
       phi_ = 2
       coordinate=cylindrical
    case ("polar_2.5D")
      if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      phi_ = 2
      z_   = 3
      coordinate=cylindrical
    case ("spherical","spherical_2D","spherical_3D")
      ndir = ndim
      r_   = 1
      ! Note that spherical 1.5D is not supported
      if(ndir>=2) theta_ = 2
      if(ndir==3) phi_ = 3
      z_   = -1
      coordinate=spherical
    case ("spherical_2.5D")
      if (ndim /= 2) &
           call mpistop("Geometry spherical_2.5D requires ndim == 2")
      ndir    = 3
      r_      = 1
      theta_  = 2
      phi_    = 3
      z_      = -1
      coordinate=spherical
    case default
      call mpistop("Unknown geometry specified")
    end select
  end subroutine set_coordinate_system

  subroutine set_pole
    use mod_global_parameters

    select case (coordinate)
    case (spherical) {^IFTHREED
      ! For spherical grid, check whether phi-direction is periodic
      if(periodB(ndim)) then
        if(phi_/=3) call mpistop("phi_ should be 3 in 3D spherical coord!")
        if(mod(ng3(1),2)/=0) &
             call mpistop("Number of meshes in phi-direction should be even!")
        if(abs(xprobmin2)<smalldouble) then
          if(mype==0) write(unitterm,*) &
               "Will apply pi-periodic conditions at northpole!"
          poleB(1,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no northpole!"
        end if
        if(abs(xprobmax2-dpi)<smalldouble) then
          if(mype==0) write(unitterm,*) &
               "Will apply pi-periodic conditions at southpole!"
          poleB(2,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no southpole!"
        end if
      end if}
    case (cylindrical)
      {
      if (^D == phi_ .and. periodB(^D)) then
        if(mod(ng^D(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if\}
    end select

  end subroutine set_pole

  !> Deallocate geometry-related variables
  subroutine putgridgeo(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    deallocate(ps(igrid)%surfaceC,ps(igrid)%surface,ps(igrid)%dvolume,ps(igrid)%dx,ps(igrid)%christoffel,&
         psc(igrid)%dx,ps(igrid)%ds,psc(igrid)%dvolume)

  end subroutine putgridgeo

  !> calculate area of surfaces of cells
  subroutine get_surface_area(s,ixG^L)
    use mod_global_parameters

    type(state) :: s
    integer, intent(in) :: ixG^L

    double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)

    select case (coordinate)
    case (Cartesian,Cartesian_stretched)
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)}

      {^IFONED
      s%surfaceC(ixG^S,1)=1.d0
      s%surface(ixG^S,1) =1.d0
      }
      {^IFTWOD
      s%surfaceC(ixG^S,1)=dx2(ixG^S)
      s%surfaceC(ixG^S,2)=drs(ixG^S)
      s%surface(ixG^S,1) =dx2(ixG^S)
      s%surface(ixG^S,2)=drs(ixG^S)
      }
      {^IFTHREED
      s%surfaceC(ixG^S,1)= dx2(ixG^S)*dx3(ixG^S)
      s%surfaceC(ixG^S,2)= drs(ixG^S)*dx3(ixG^S)
      s%surfaceC(ixG^S,3)= drs(ixG^S)*dx2(ixG^S)
      s%surface(ixG^S,1)=s%surfaceC(ixG^S,1)
      s%surface(ixG^S,2)=s%surfaceC(ixG^S,2)
      s%surface(ixG^S,3)=s%surfaceC(ixG^S,3)
      }
      {s%surfaceC(0^D%ixG^S,^D)=s%surfaceC(1^D%ixG^S,^D);\}
    case (spherical)
      x(ixG^S,1)=s%x(ixG^S,1)
      {^NOONED
      x(ixG^S,2)=s%x(ixG^S,2)
      }
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)
      }
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)
      }

      s%surfaceC(ixG^S,1)=(x(ixG^S,1)+0.5d0*drs(ixG^S))**2 {^NOONED &
           *2.0d0*dsin(x(ixG^S,2))*dsin(0.5d0*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}

      {^NOONED
      s%surfaceC(ixG^S,2)=(x(ixG^S,1)**2 +drs(ixG^S)**2/12.0d0)*drs(ixG^S)&
           *dsin(x(ixG^S,2)+0.5d0*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}

      {^IFTHREED
      s%surfaceC(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      }

      {^IFONED
      s%surfaceC(0,1)=dabs(x(1,1)-0.5d0*drs(1))**2
      }
      {^IFTWOD
      s%surfaceC(0,ixGmin2:ixGmax2,1)=(x(1,ixGmin2:ixGmax2,1)-0.5d0*drs(1,ixGmin2:ixGmax2))**2 &
         *2.0d0*dsin(x(1,ixGmin2:ixGmax2,2))*dsin(0.5d0*dx2(1,ixGmin2:ixGmax2))

      s%surfaceC(ixGmin1:ixGmax1,0,2)=(x(ixGmin1:ixGmax1,1,1)**2 +drs(ixGmin1:ixGmax1,1)**2/12.0d0) *drs(ixGmin1:ixGmax1,1)&
                       *dsin(x(ixGmin1:ixGmax1,1,2)-0.5d0*dx2(ixGmin1:ixGmax1,1))
      }
      {^IFTHREED
      s%surfaceC(0,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)=(x(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1)-0.5d0*drs(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3))**2*2.0d0*dsin(x(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2))*dsin(0.5d0*dx2(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))*dx3(1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3)
      s%surfaceC(ixGmin1:ixGmax1,0,ixGmin3:ixGmax3,2)=x(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3,1)*drs(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3)*dsin(x(ixGmin1:ixGmax1,1,ixGmin3:ixGmax3,&
         2)-0.5d0*dx2(ixGmin1:ixGmax1,1,ixGmin3:ixGmax3))*dx3(ixGmin1:ixGmax1,1,&
         ixGmin3:ixGmax3)
      s%surfaceC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,0,3)=&
         s%surfaceC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1,3)
      }

      s%surface(ixG^S,1)=x(ixG^S,1)**2 {^NOONED &
           *2.0d0*dsin(x(ixG^S,2))*dsin(0.5d0*dx2(ixG^S))}{^IFTHREED*dx3(ixG^S)}
      {^NOONED
      s%surface(ixG^S,2)=x(ixG^S,1)*drs(ixG^S)&
           *dsin(x(ixG^S,2))}{^IFTHREED*dx3(ixG^S)}

      {^IFTHREED
      s%surface(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)}

    case (cylindrical)
      x(ixG^S,1)=s%x(ixG^S,1)
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)}

      s%surfaceC(ixG^S,1)=(x(ixG^S,1)+0.5d0*drs(ixG^S)) {^NOONED *dx2(ixG^S) } {^IFTHREED *dx3(ixG^S) }
      {^NOONED
      s%surfaceC(ixG^S,2)=x(ixG^S,1)*drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      }
      {^IFTHREED
      s%surfaceC(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      }

      {^IFONED
      s%surfaceC(0,1)=dabs(x(1,1)-0.5d0*drs(1))
      }
      {^IFTWOD
      s%surfaceC(0,ixGmin2:ixGmax2,1)=dabs(x(1,ixGmin2:ixGmax2,1)-0.5d0*drs(1,&
         ixGmin2:ixGmax2))*dx2(1,ixGmin2:ixGmax2)
      }
      {^IFTHREED
      s%surfaceC(0,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)=dabs(x(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1)-0.5d0*drs(1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))*dx2(1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3)*dx3(1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3)
      }
      {s%surfaceC(0^DE%ixG^S,^DE)=s%surfaceC(1^DE%ixG^S,^DE);\}

      s%surface(ixG^S,1)=(x(ixG^S,1)) {^NOONED *dx2(ixG^S) } {^IFTHREED *dx3(ixG^S) }
      {^NOONED
      s%surface(ixG^S,2)=x(ixG^S,1)*drs(ixG^S){^IFTHREED*dx3(ixG^S)}
      }
      {^IFTHREED
      s%surface(ixG^S,3)=x(ixG^S,1)*drs(ixG^S)*dx2(ixG^S)
      }

    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_surface_area

  !> calculate christoffel symbols of cells after calculating the surface area and volume
  subroutine get_christoffel(s,ixG^L)
    use mod_global_parameters

    type(state) :: s
    integer, intent(in) :: ixG^L

    double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)
    integer          :: h1x^L{^NOONED, h2x^L}

    s%christoffel=0.0d0
    select case (coordinate)
    case (Cartesian,Cartesian_stretched)
       return ! nothing to do here
    case (spherical)
      h1x^L=ixG^L-kr(1,^D); {^NOONED h2x^L=ixG^L-kr(2,^D);}
      x(ixG^S,1)=s%x(ixG^S,1)
      {^NOONED
      x(ixG^S,2)=s%x(ixG^S,2)
      }
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)
      }
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)
      }

      s%christoffel(ixG^S,2,1,2)=0.5d0 * (s%surfaceC(ixG^S,1) - s%surfaceC(h1x^S,1)) / s%dvolume(ixG^S)
      s%christoffel(ixG^S,2,2,1)=s%christoffel(ixG^S,2,1,2)
      s%christoffel(ixG^S,3,3,1)=s%christoffel(ixG^S,2,1,2)
      s%christoffel(ixG^S,3,1,3)=s%christoffel(ixG^S,2,1,2)

      s%christoffel(ixG^S,1,2,2)=-0.25d0 * ((x(ixG^S,1)+0.5d0*drs(ixG^S))**2*s%surfaceC(ixG^S,1) &
                                          - (x(ixG^S,1)-0.5d0*drs(ixG^S))**2*s%surfaceC(h1x^S,1)) / s%dvolume(ixG^S)

      s%christoffel(ixG^S,1,3,3)=-0.25d0 / s%dvolume(ixG^S) &
                           * ((x(ixG^S,1)+0.5d0*drs(ixG^S))**4 - (x(ixG^S,1)-0.5d0*drs(ixG^S))**4){^NOONED &
           *( (dcos(x(ixG^S,2)+0.5d0*dx2(ixG^S))**3/3.0d0 - dcos(x(ixG^S,2)+0.5d0*dx2(ixG^S))) &
             -(dcos(x(ixG^S,2)-0.5d0*dx2(ixG^S))**3/3.0d0 - dcos(x(ixG^S,2)-0.5d0*dx2(ixG^S))) )}{^IFTHREED*dx3(ixG^S)}

      {^NOONED
      s%christoffel(ixG^S,2,3,3)= - ( &
            (dsin(x(ixG^S,2)+0.5d0*dx2(ixG^S))**2)*s%surfaceC(ixG^S,2) &
          - (dsin(x(ixG^S,2)-0.5d0*dx2(ixG^S))**2)*s%surfaceC(h2x^S,2) &
              ) / 3.0d0 / s%dvolume(ixG^S) 

      s%christoffel(ixG^S,3,2,3)=dcos(x(ixG^S,2))/dsin(x(ixG^S,2))
      s%christoffel(ixG^S,3,3,2)=s%christoffel(ixG^S,3,2,3)
      }
    case (cylindrical)
      x(ixG^S,1)=s%x(ixG^S,1)
      drs(ixG^S)=s%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=s%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=s%dx(ixG^S,3)}

      if ( z_ == 2 ) then
         s%christoffel(ixG^S,1,3,3)= - (x(ixG^S,1)**2 +drs(ixG^S)**2/12.0d0) / x(ixG^S,1)
         s%christoffel(ixG^S,3,1,3)= 1.0d0 / x(ixG^S,1)
         s%christoffel(ixG^S,3,3,1) = s%christoffel(ixG^S,3,1,3)
      else
         s%christoffel(ixG^S,1,2,2)= - (x(ixG^S,1)**2 +drs(ixG^S)**2/12.0d0) / x(ixG^S,1)
         s%christoffel(ixG^S,2,1,2)= 1.0d0 / x(ixG^S,1)
         s%christoffel(ixG^S,2,2,1) = s%christoffel(ixG^S,2,1,2)
      end if
    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_christoffel

  !> get the 3-metric in flat space sqrt_gamma_hat 
  subroutine get_sqrt_gamma_hat(x, ixI^L, ixO^L, sqrt_gamma)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out) :: sqrt_gamma(ixI^S)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    select case (coordinate)
    case (Cartesian)
      sqrt_gamma(ixO^S) = 1.0d0
    case (cylindrical)
      sqrt_gamma(ixO^S) = dabs( x(ixO^S, 1) )
    case (spherical)
      sqrt_gamma(ixO^S) = dabs( x(ixO^S, 1)**2 {^NOONED * dsin(x(ixO^S, 2)) } )
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_sqrt_gamma_hat

  !> get the 3-metric in flat space gamma_ij_hat 
  subroutine get_gamma_ij_hat(x, ixI^L, ixO^L, gamma)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(out) :: gamma(ixI^S, 1:3, 1:3)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    gamma = 0.0d0

    select case (coordinate)
    case (Cartesian)
      gamma(ixO^S,1,1) = 1.0d0
      gamma(ixO^S,2,2) = 1.0d0
      gamma(ixO^S,3,3) = 1.0d0
    case (cylindrical)
      gamma(ixO^S,1,1) = 1.0d0
      if ( z_ == 2 ) then
         gamma(ixO^S,2,2) = 1.0d0
         gamma(ixO^S,3,3) = x(ixO^S, 1)**2
      else
         gamma(ixO^S,3,3) = 1.0d0
         gamma(ixO^S,2,2) = x(ixO^S, 1)**2
      end if
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      gamma(ixO^S,1,1) = 1.0d0
      gamma(ixO^S,2,2) = x(ixO^S, 1)**2
      gamma(ixO^S,3,3) =  gamma(ixO^S,2,2) {^NOONED * dsin(x(ixO^S, 2))**2 }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_gamma_ij_hat
! add BHAC
  !> transform vectors between natural and orthonomal basis 
  subroutine get_natural2orthonormal(x, ixI^L, ixO^L, N2R)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: N2R(ixI^S, 1:3)

    N2R(ixO^S,1:3) = 1.0d0

    select case (coordinate)
    case (Cartesian)
      ! nothing to do here
    case (cylindrical)
      if ( z_ == 2 ) then
         N2R(ixO^S,3) = x(ixO^S, 1)
      else
         N2R(ixO^S,2) = x(ixO^S, 1)
      end if
    case (spherical)
      ! note: r_ = 1,theta_=2, phi_ = 3.
      N2R(ixO^S,2) = x(ixO^S, 1)
      N2R(ixO^S,3) = x(ixO^S, 1) {^NOONED * dsin(x(ixO^S, 2)) }
    case default
      call mpistop("Sorry, coordinate unknown")
    end select
  end subroutine get_natural2orthonormal

  !> Calculate dervitives of a scalar q within ixL in direction idir
  subroutine partial_d(q,ixI^L,ixO^L,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(inout) :: gradq(ixI^S)
    double precision                :: x(ixI^S,1:ndim)
    integer                         :: jxO^L, hxO^L

    x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim)

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=0.5d0*(q(jxO^S)-q(hxO^S))/dxlevel(idir)
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,2)-x(hxO^S,2))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,3)-x(hxO^S,3))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,phi_)-x(hxO^S,phi_))!*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine partial_d

  !> Calculate covariant dervitives (based on flat metric) of a scalar q within ixL in direction idir
  !> \nabla_i X^j = \partial_i X^j + \Gamma^j_{ik} X^k
  subroutine covariant_dervitive(qvec,ixI^L,ixO^L,dqvec)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qvec(ixI^S,1:3)
    double precision, intent(inout) :: dqvec(ixI^S,1:3,1:3)
    !double precision                :: x(ixI^S,1:ndim)
    double precision                :: chris(ixI^S,1:3,1:3,1:3)
    integer                         :: jxO^L, hxO^L
    integer                         :: idir, jdir, kdir

    !x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim)
    associate(x=>block%x)
    
    select case(coordinate)

    case(Cartesian)
       dqvec=0.0d0
       do idir = 1, ndim
          do jdir = 1, 3
             hxO^L=ixO^L-kr(idir,^D);
             jxO^L=ixO^L+kr(idir,^D);
             dqvec(ixO^S,idir,jdir)=0.5d0*(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/dxlevel(idir)
          end do
       end do

    case(Cartesian_stretched)
       dqvec=0.0d0
       do idir = 1, ndim
          do jdir = 1, 3
             hxO^L=ixO^L-kr(idir,^D);
             jxO^L=ixO^L+kr(idir,^D);
             dqvec(ixO^S,idir,jdir)=(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/(x(jxO^S,idir)-x(hxO^S,idir))
          end do
       end do

     case(spherical)
       chris = 0.0d0
       chris(ixO^S,1, 2, 2) = - x(ixO^S,1)
       chris(ixO^S,1, 3, 3) = - x(ixO^S,1) {^NOONED * dsin(x(ixO^S,2)) }
       chris(ixO^S,2, 1, 2) = 1.0d0 / x(ixO^S,1) 
       chris(ixO^S,2, 2, 1) = chris(ixO^S,2, 1, 2)
       chris(ixO^S,3, 3, 1) = chris(ixO^S,2, 1, 2)
       chris(ixO^S,3, 1, 3) = chris(ixO^S,2, 1, 2)
       {^NOONED 
       chris(ixO^S,2, 3, 3) = - dsin(x(ixO^S,2)) * dcos(x(ixO^S,2)) 
       chris(ixO^S,3, 2, 3) = dcos(x(ixO^S,2)) / dsin(x(ixO^S,2)) 
       chris(ixO^S,3, 3, 2) = chris(ixO^S,3, 2, 3)
       }
       dqvec=0.0d0
       do idir = 1, 3
          do jdir = 1, 3
             hxO^L=ixO^L-kr(idir,^D);
             jxO^L=ixO^L+kr(idir,^D);
             select case(idir)
             case(1)
               dqvec(ixO^S,jdir,idir)=(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/((x(jxO^S,1)-x(hxO^S,1)))
               {^NOONED
             case(2)
               dqvec(ixO^S,jdir,idir)=(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/(x(jxO^S,2)-x(hxO^S,2))
               }
               {^IFTHREED
             case(3)
               dqvec(ixO^S,jdir,idir)=(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/(x(jxO^S,3)-x(hxO^S,3))
               }
             end select
             ! connection terms
             do kdir = 1, 3
                dqvec(ixO^S,jdir,idir)=dqvec(ixO^S,jdir,idir) &
                              + chris(ixO^S,jdir,idir,kdir)*qvec(ixO^S, kdir)
             end do
          end do
       end do

    case(cylindrical)
       dqvec=0.0d0
       do idir = 1, 3
          do jdir = 1, 3
             hxO^L=ixO^L-kr(idir,^D);
             jxO^L=ixO^L+kr(idir,^D);
             if(idir==phi_) then
               dqvec(ixO^S,idir,jdir)=(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/(x(jxO^S,phi_)-x(hxO^S,phi_))!*x(ixO^S,r_))
             else
               dqvec(ixO^S,idir,jdir)=(qvec(jxO^S,jdir)-qvec(hxO^S,jdir))/(x(jxO^S,idir)-x(hxO^S,idir))
             end if
             do kdir = 1,3
                dqvec(ixO^S,idir,jdir)=dqvec(ixO^S,idir,jdir) &
                              + chris(ixO^S,jdir,idir,kdir)*qvec(ixO^S, kdir)
             end do
             call mpistop('cov der cylindrical not ready yet')
          end do
       end do

    case default
      call mpistop('Unknown geometry')
    end select

    end associate
  end subroutine covariant_dervitive

  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixI^L,ixO^L,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(inout) :: gradq(ixI^S)
    double precision                :: x(ixI^S,1:ndim)
    integer                         :: jxO^L, hxO^L

    x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim)

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=0.5d0*(q(jxO^S)-q(hxO^S))/dxlevel(idir)
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,3)-x(hxO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,phi_)-x(hxO^S,phi_))*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine gradient

  !> Calculate gradient of a scalar q in direction idir at cell interfaces
  subroutine gradientx(q,x,ixI^L,ixO^L,idir,gradq,fourth_order)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S), x(ixI^S,1:ndim)
    double precision, intent(inout) :: gradq(ixI^S)
    logical, intent(in)             :: fourth_order
    integer                         :: jxO^L, hxO^L, kxO^L

    if(fourth_order) then
      ! Fourth order, stencil width is two
      kxO^L=ixO^L^LADD2;
      if(ixImin^D>kxOmin^D.or.ixImax^D<kxOmax^D|.or.) &
           call mpistop("Error in gradientx: Non-conforming input limits")
      hxO^L=ixO^L-kr(idir,^D);
      jxO^L=ixO^L+kr(idir,^D);
      kxO^L=ixO^L+kr(idir,^D)*2;
    else
      hxO^L=ixO^L;
    end if
    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      if(fourth_order) then
        gradq(ixO^S)=(27.d0*(q(jxO^S)-q(ixO^S))-q(kxO^S)+q(hxO^S))/24.d0/dxlevel(idir)
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/dxlevel(idir)
      end if
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,3)-x(hxO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)))
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,phi_)-x(hxO^S,phi_))*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine gradientx

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> first use limiter to go from cell center to edge
  subroutine gradientS(q,ixI^L,ixO^L,idir,gradq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixI^L, ixO^L, idir
    double precision, intent(in)       :: q(ixI^S)
    double precision, intent(inout)    :: gradq(ixI^S)
    double precision ,dimension(ixI^S) :: qC,qL,qR,dqC,ldq,rdq

    double precision :: x(ixI^S,1:ndim)
    double precision :: invdx
    integer          :: hxO^L,ixC^L,jxC^L,gxC^L,hxC^L

    x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim)

    invdx=1.d0/dxlevel(idir)
    hxO^L=ixO^L-kr(idir,^D);
    ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
    jxC^L=ixC^L+kr(idir,^D);
    gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
    hxC^L=gxC^L+kr(idir,^D);

    ! set the gradient limiter here
    qR(gxC^S) = q(hxC^S)
    qL(gxC^S) = q(gxC^S)
    if (typegradlimiter/=limiter_ppm) then
      dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
      call dwlimiter2(dqC,ixI^L,gxC^L,idir,typegradlimiter,ldw=ldq,rdw=rdq)
      qL(ixC^S) = qL(ixC^S) + 0.5d0*ldq(ixC^S)
      qR(ixC^S) = qR(ixC^S) - 0.5d0*rdq(jxC^S)
    else
      call PPMlimitervar(ixI^L,ixO^L,idir,q,q,qL,qR)
    endif

    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))*invdx
    case(Cartesian_stretched)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/block%dx(ixO^S,idir)
    case(spherical)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/block%dx(ixO^S,idir)
      select case(idir)
      case(2)
        gradq(ixO^S)=gradq(ixO^S)/x(ixO^S,1)
        {^IFTHREED
      case(3)
        gradq(ixO^S)=gradq(ixO^S)/(x(ixO^S,1)*dsin(x(ixO^S,2)))
        }
      end select
    case(cylindrical)
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/block%dx(ixO^S,idir)
      if(idir==phi_) gradq(ixO^S)=gradq(ixO^S)/x(ixO^S,1)
    end select

  end subroutine gradientS

  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixI^L,ixO^L,divq, fourthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qvec(ixI^S,1:ndir)
    double precision, intent(inout) :: divq(ixI^S)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical                         :: use_4th_order
    double precision                :: qC(ixI^S), invdx(1:ndim)
    integer                         :: jxO^L, hxO^L, ixC^L, jxC^L
    integer                         :: idims, ix^L, gxO^L, kxO^L

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab_uniform) &
           call mpistop("divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ix^L=ixO^L^LADD2;
    else
      ! Second order, stencil width is one
      ix^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=0.0d0

    if (slab_uniform) then
      do idims=1,ndim
        if (.not. use_4th_order) then
          ! Use second order scheme
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+0.5d0*(qvec(jxO^S,idims) &
               - qvec(hxO^S,idims))*invdx(idims)
        else
          ! Use fourth order scheme
          kxO^L=ixO^L+2*kr(idims,^D);
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          gxO^L=ixO^L-2*kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+&
               (-qvec(kxO^S,idims) + 8.0d0 * qvec(jxO^S,idims) - 8.0d0 * &
               qvec(hxO^S,idims) + qvec(gxO^S,idims))/(12.0d0 * dxlevel(idims))
        end if
      end do
    else
      do idims=1,ndim
        hxO^L=ixO^L-kr(idims,^D);
        ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
        jxC^L=ixC^L+kr(idims,^D);
        if(stretched_dim(idims) .and. stretch_uncentered) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixC^S)=block%surfaceC(ixC^S,idims)*(qvec(ixC^S,idims)+0.5d0*block%dx(ixC^S,idims)*&
               (qvec(jxC^S,idims)-qvec(ixC^S,idims))/(block%x(jxC^S,idims)-block%x(ixC^S,idims)))
        else
          qC(ixC^S)=block%surfaceC(ixC^S,idims)*0.5d0*(qvec(ixC^S,idims)+qvec(jxC^S,idims))
        end if
        divq(ixO^S)=divq(ixO^S)+qC(ixO^S)-qC(hxO^S)
      end do
      divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)
    end if

  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixI^L,ixO^L,divq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: qvec(ixI^S,1:ndir)
    double precision, intent(inout)    :: divq(ixI^S)
    double precision, dimension(ixI^S) :: qL,qR,dqC,ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxO^L,ixC^L,jxC^L,idims,ix^L,gxC^L,hxC^L

    ix^L=ixO^L^LADD2;

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=zero
    do idims=1,ndim
      hxO^L=ixO^L-kr(idims,^D);
      ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
      jxC^L=ixC^L+kr(idims,^D);
      gxCmin^D=ixCmin^D-kr(idims,^D);gxCmax^D=jxCmax^D;
      hxC^L=gxC^L+kr(idims,^D);
      qR(gxC^S) = qvec(hxC^S,idims)
      qL(gxC^S) = qvec(gxC^S,idims)
      if(typegradlimiter/=limiter_ppm) then
        dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
        call dwlimiter2(dqC,ixI^L,gxC^L,idims,typegradlimiter,ldw=ldq,rdw=rdq)
        qL(ixC^S) = qL(ixC^S) + 0.5d0*ldq(ixC^S)
        qR(ixC^S) = qR(ixC^S) - 0.5d0*rdq(jxC^S)
      else
        dqC(ixI^S)=qvec(ixI^S,idims)
        call PPMlimitervar(ixI^L,ixO^L,idims,dqC,dqC,qL,qR)
      endif

      if (slab_uniform) then
        divq(ixO^S)=divq(ixO^S)+0.5d0*(qR(ixO^S)-qL(hxO^S))*invdx(idims)
      else
        qR(ixC^S)=block%surfaceC(ixC^S,idims)*qR(ixC^S)
        qL(ixC^S)=block%surfaceC(ixC^S,idims)*qL(ixC^S)
        divq(ixO^S)=divq(ixO^S)+qR(ixO^S)-qL(hxO^S)
      end if
    end do
    if(.not.slab_uniform) divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)

  end subroutine divvectorS

  !> Calculate curl of a vector qvec within ixL
  !> Options to 
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixI^S,1:ndir0)
    double precision, intent(inout) :: curlvec(ixI^S,idirmin0:3)

    integer          :: ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixI^S),tmp2(ixI^S),xC(ixI^S),surface(ixI^S)

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixO^S,idirmin0:3)=zero

    ! all non-Cartesian cases 
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            tmp(ixI^S)=qvec(ixI^S,kdir)
            hxO^L=ixO^L-kr(jdir,^D);
            jxO^L=ixO^L+kr(jdir,^D);
            ! second order centered differencing 
            tmp(ixO^S)=0.5d0*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
            !> \todo allow for 4th order CD evaluation here as well
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(typecurl) 
              case('central') 
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                ! second order centered differencing 
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,jdir)-block%x(hxO^S,jdir))
              case('Gaussbased') 
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=block%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%dvolume(ixO^S)
              case('Stokesbased') 
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                if(kdir<=ndim)then
                  tmp(ixC^S)=block%ds(ixO^S,kdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                else
                  tmp(ixC^S)=(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%surface(ixO^S,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/(block%ds(ixO^S,jdir)*block%ds(ixO^S,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(typecurl) 
          case('central') ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                select case(jdir)
                case(1)
                tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,1)-block%x(hxO^S,1))*block%x(ixO^S,1))
                {^NOONED    case(2)
                if(idir==1) tmp(ixI^S)=tmp(ixI^S)*dsin(block%x(ixI^S,2))
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,2)-block%x(hxO^S,2))*block%x(ixO^S,1))
                if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/dsin(block%x(ixO^S,2))
                }
                {^IFTHREED  case(3)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,3)-block%x(hxO^S,3))*block%x(ixO^S,1)*dsin(block%x(ixO^S,2)))
                }
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased') 
            do idir=idirmin0,3; 
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=block%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms 
            if(idir==2.and.phi_>0) curlvec(ixO^S,2)=curlvec(ixO^S,2)+qvec(ixO^S,phi_)/block%x(ixO^S,r_)
            {^NOONED
            if(idir==phi_) curlvec(ixO^S,phi_)=curlvec(ixO^S,phi_)-qvec(ixO^S,2)/block%x(ixO^S,r_) &
                 +qvec(ixO^S,r_)*dcos(block%x(ixO^S,2))/(block%x(ixO^S,r_)*dsin(block%x(ixO^S,2)))
            }
            enddo; 
          case('Stokesbased') 
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(dsin(xC(ixO^S))*tmp(ixO^S)-&
                       dsin(xC(hxO^S))*tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)*block%x(ixO^S,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,1)
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(hxO^S)*tmp(hxO^S)-xC(ixO^S)*tmp(ixO^S))*&
                       dsin(block%x(ixO^S,idir))*block%dx(ixO^S,kdir))/block%surface(ixO^S,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  if(ndim==3) then
                    surface(ixO^S)=block%surface(ixO^S,idir)
                  else
                    surface(ixO^S)=block%x(ixO^S,jdir)*block%dx(ixO^S,kdir)*block%dx(ixO^S,jdir)
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /surface(ixO^S)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(typecurl) 
          case('central')  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                if(z_==3.or.z_==-1) then 
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,1)-block%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/block%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR 
                  {^NOONED      case(2)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,2)-block%x(hxO^S,2))*block%x(ixO^S,1))
                  }
                  {^IFTHREED    case(3)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,3)-block%x(hxO^S,3))
                  }
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,1)-block%x(hxO^S,1)) 
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/block%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR
                  {^NOONED      case(2)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,2)-block%x(hxO^S,2))
                  }
                  {^IFTHREED    case(3)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,3)-block%x(hxO^S,3))*block%x(ixO^S,1))
                  }
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased') ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3; 
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=block%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-qvec(ixO^S,z_)/block%x(ixO^S,r_)
              ! polar
              if(phi_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+qvec(ixO^S,z_)/block%x(ixO^S,r_)
            endif
            enddo; 
          case('Stokesbased')
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%x(ixO^S,idir)*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,jdir)
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
    end select

  end subroutine curlvector

  !> Calculate idim's transverse components of curl of a vector qvec within ixL
  !> Options to 
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector_trans(qvec,qvecc,ixI^L,ixO^L,curlvec,idim,idirmin,idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    integer, intent(in)             :: idim, ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixI^S,1:ndir0),qvecc(ixI^S,1:ndir0)
    double precision, intent(inout) :: curlvec(ixI^S,idirmin0:3)

    integer          :: ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixI^S),tmp2(ixI^S),xC(ixI^S),surface(ixI^S)

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixO^S,idirmin0:3)=zero

    ! all non-Cartesian cases 
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3 
          do jdir=1,ndim; do kdir=1,ndir0
            if(lvc(idir,jdir,kdir)/=0)then
              if(jdir/=idim) then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
              else
                ! use two cell-center values for gradient at interface of the two cells
                ! because left and right neighbor interface values is unavailable at block boundary faces
                tmp(ixI^S)=qvecc(ixI^S,kdir)
                hxO^L=ixO^L;
                jxO^L=ixO^L+kr(jdir,^D);
              end if
              ! second order centered differencing 
              tmp(ixO^S)=0.5d0*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
              !> \todo allow for 4th order CD evaluation here as well
              if(lvc(idir,jdir,kdir)==1)then
                curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
              else
                curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
              endif
              if(idir<idirmin)idirmin=idir
            endif
          enddo; enddo;
        end do
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(typecurl) 
              case('central') 
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                ! second order centered differencing 
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,jdir)-block%x(hxO^S,jdir))
              case('Gaussbased') 
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=block%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%dvolume(ixO^S)
              case('Stokesbased') 
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                if(kdir<=ndim)then
                  tmp(ixC^S)=block%ds(ixO^S,kdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                else
                  tmp(ixC^S)=(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%surface(ixO^S,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/(block%ds(ixO^S,jdir)*block%ds(ixO^S,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(typecurl) 
          case('central') ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                select case(jdir)
                case(1)
                tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,1)-block%x(hxO^S,1))*block%x(ixO^S,1))
                {^NOONED    case(2)
                if(idir==1) tmp(ixI^S)=tmp(ixI^S)*dsin(block%x(ixI^S,2))
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,2)-block%x(hxO^S,2))*block%x(ixO^S,1))
                if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/dsin(block%x(ixO^S,2))
                }
                {^IFTHREED  case(3)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,3)-block%x(hxO^S,3))*block%x(ixO^S,1)*dsin(block%x(ixO^S,2)))
                }
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased') 
            do idir=idirmin0,3; 
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=block%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms 
            if(idir==2.and.phi_>0) curlvec(ixO^S,2)=curlvec(ixO^S,2)+qvec(ixO^S,phi_)/block%x(ixO^S,r_)
            {^NOONED
            if(idir==phi_) curlvec(ixO^S,phi_)=curlvec(ixO^S,phi_)-qvec(ixO^S,2)/block%x(ixO^S,r_) &
                 +qvec(ixO^S,r_)*dcos(block%x(ixO^S,2))/(block%x(ixO^S,r_)*dsin(block%x(ixO^S,2)))
            }
            enddo; 
          case('Stokesbased') 
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(dsin(xC(ixO^S))*tmp(ixO^S)-&
                       dsin(xC(hxO^S))*tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)*block%x(ixO^S,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,1)
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(hxO^S)*tmp(hxO^S)-xC(ixO^S)*tmp(ixO^S))*&
                       dsin(block%x(ixO^S,idir))*block%dx(ixO^S,kdir))/block%surface(ixO^S,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  if(ndim==3) then
                    surface(ixO^S)=block%surface(ixO^S,idir)
                  else
                    surface(ixO^S)=block%x(ixO^S,jdir)*block%dx(ixO^S,kdir)*block%dx(ixO^S,jdir)
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /surface(ixO^S)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(typecurl) 
          case('central')  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixI^S)=qvec(ixI^S,kdir)
                hxO^L=ixO^L-kr(jdir,^D);
                jxO^L=ixO^L+kr(jdir,^D);
                if(z_==3.or.z_==-1) then 
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,1)-block%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/block%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR 
                  {^NOONED      case(2)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,2)-block%x(hxO^S,2))*block%x(ixO^S,1))
                  }
                  {^IFTHREED    case(3)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,3)-block%x(hxO^S,3))
                  }
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1) ! R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,1)-block%x(hxO^S,1)) 
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/block%x(ixO^S,1) ! (1/R)*d(R V_phi)/dR
                  {^NOONED      case(2)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,2)-block%x(hxO^S,2))
                  }
                  {^IFTHREED    case(3)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,3)-block%x(hxO^S,3))*block%x(ixO^S,1))
                  }
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased') ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3; 
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                tmp(ixC^S)=block%surfaceC(ixC^S,jdir)*(qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir)))
                tmp2(ixO^S)=(tmp(ixO^S)-tmp(hxO^S))/block%dvolume(ixO^S)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
                else
                  curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-qvec(ixO^S,z_)/block%x(ixO^S,r_)
              ! polar
              if(phi_==2) curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+qvec(ixO^S,z_)/block%x(ixO^S,r_)
            endif
            enddo; 
          case('Stokesbased')
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%x(ixO^S,idir)*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,jdir)
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
    end select

  end subroutine curlvector_trans

end module mod_geometry
