!> Moyt(2:ny-1)dule for eos
module mod_eos_interpolation

  implicit none
  public

  contains

  !---------------------------------------------------------------------
 !
  !     purpose: interpolation of a function of three variables in an
  !              equidistant(!!!) table.
  !
  !     method:  8-point Lagrange linear interpolation formula          
  !
  !     x        input vector of first  variable
  !     y        input vector of second variable
  !     z        input vector of third  variable
  !
  !     f        output vector of interpolated function values
  !
  !
  !     ft       3d array of tabulated function values
  !     nx       x-dimension of table
  !     ny       y-dimension of table
  !     nz       z-dimension of table
  !     xt       vector of x-coordinates of table
  !     yt       vector of y-coordinates of table
  !     zt       vector of z-coordinates of table
  !
  !     d1       centered derivative of ft with respect to x
  !     d2       centered derivative of ft with respect to y
  !     d3       centered derivative of ft with respect to z
  !---------------------------------------------------------------------
  subroutine intep3d(x, y, z, f, ft, nx, ny, nz, xt, yt, zt, ix_out, iy_out,&
      iz_out, d1, d2, d3)
     implicit none

!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     kt       vector length of input and output vectors
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!     d1       centered derivative of ft with respect to x
!     d2       centered derivative of ft with respect to y
!     d3       centered derivative of ft with respect to z
!     Note that d? only make sense when intp3d is called with kt=1
!---------------------------------------------------------------------



     integer, intent(in) ::  nx,ny,nz
     double precision ::  x,y,z
     double precision ::  f
     double precision :: xt(nx),yt(ny),zt(nz)
     double precision :: ft(nx,ny,nz)
     double precision, intent(out), optional ::  d1,d2,d3
!  integer :: ktx = 400
     double precision, allocatable  :: fh(:,:), delx(:), dely(:), delz(:)
     double precision, allocatable  :: a1(:), a2(:), a3(:), a4(:)
     double precision, allocatable  :: a5(:), a6(:), a7(:), a8(:)
     integer, intent(out), optional     ::  ix_out, iy_out, iz_out

     double precision :: dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer :: n,ix,iy,iz
     double precision :: a,b,c

     allocate(fh(1,8))

     allocate(delx(1))
     allocate(dely(1))
     allocate(delz(1))
     allocate(a1(1))
     allocate(a2(1))
     allocate(a3(1))
     allocate(a4(1))
     allocate(a5(1))
     allocate(a6(1))
     allocate(a7(1))
     allocate(a8(1))

! write(*,*)  nx,ny,nz, "x,y,z"
! write(*,*) shape(ft)
     a = float(nx-1)
     b = float(ny-1)
     c = float(nz-1)

! ! !------  determine spacing parameters of (equidistant!!!) table
      dx    = (xt(nx) - xt(1)) / a !FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / b !FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / c !FLOAT(nz-1)


!     write(*,*) dx,dy,dz ,"dxdydz"
      dxi   = 1.0d0 / dx
      dyi   = 1.0d0 / dy
      dzi   = 1.0d0 / dz

      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi

      dxyzi = dxi * dyi * dzi
! ! ------- loop over all points to be interpolated
!      DO  n = 1, kt
! ------- determine location in (equidistant!!!) table
       ix = 2 + INT( (x - xt(1) - 1.d-10) * dxi )
       iy = 2 + INT( (y - yt(1) - 1.d-10) * dyi )
       iz = 2 + INT( (z - zt(1) - 1.d-10) * dzi )

       ix = MAX( 2, MIN( ix, nx ) )
       iy = MAX( 2, MIN( iy, ny ) )
       iz = MAX( 2, MIN( iz, nz ) )

     if(present(ix_out)) ix_out = ix
     if(present(iy_out)) iy_out = iy
     if(present(iz_out)) iz_out = iz

!        write(*,*) iy-1,iy,iy+1
!------- set-up auxiliary arrays for Lagrange interpolation

       delx(1) = xt(ix) - x
       dely(1) = yt(iy) - y
       delz(1) = zt(iz) - z

       fh(1,1) = ft(ix  , iy  , iz  )
       fh(1,2) = ft(ix-1, iy  , iz  )
       fh(1,3) = ft(ix  , iy-1, iz  )
       fh(1,4) = ft(ix  , iy  , iz-1)
       fh(1,5) = ft(ix-1, iy-1, iz  )
       fh(1,6) = ft(ix-1, iy  , iz-1)
       fh(1,7) = ft(ix  , iy-1, iz-1)
       fh(1,8) = ft(ix-1, iy-1, iz-1)

!------ set up coefficients of the interpolation polynomial and
!       evaluate function values

       a1(1) = fh(1,1)
       a2(1) = dxi   * ( fh(1,2) - fh(1,1) )
       a3(1) = dyi   * ( fh(1,3) - fh(1,1) )
       a4(1) = dzi   * ( fh(1,4) - fh(1,1) )
       a5(1) = dxyi  * ( fh(1,5) - fh(1,2) - fh(1,3) + fh(1,1) )
       a6(1) = dxzi  * ( fh(1,6) - fh(1,2) - fh(1,4) + fh(1,1) )
       a7(1) = dyzi  * ( fh(1,7) - fh(1,3) - fh(1,4) + fh(1,1) )
       a8(1) = dxyzi * ( fh(1,8) - fh(1,1) + fh(1,2) + fh(1,3) +fh(1,4) - fh(1,&
          5) - fh(1,6) - fh(1,7) )

      
     if(present(d1)) d1 = -a2(1)
     if(present(d2)) d2 = -a3(1)
     if(present(d3)) d3  = -a4(1)

       f  = a1(1) +  a2(1) * delx(1)  +  a3(1) * dely(1)    +  a4(1) * delz(1) &
             +  a5(1) * delx(1) * dely(1)       +  a6(1) * delx(1) * delz(1)   &
              +  a7(1) * dely(1) * delz(1)       +  a8(1) * delx(1) * dely(1) &
          * delz(1)

!      ENDDO


     deallocate(fh)
     deallocate(delx)
     deallocate(dely)
     deallocate(delz)
     deallocate(a1)
     deallocate(a2)
     deallocate(a3)
     deallocate(a4)
     deallocate(a5)
     deallocate(a6)
     deallocate(a7)
     deallocate(a8)


      RETURN
      END subroutine 

      SUBROUTINE intep3d_many ( x, y, z, f, ft, nx, ny, nz, nvars, xt, yt, zt)
!
      implicit none
!
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     kt       vector length of input and output vectors
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!---------------------------------------------------------------------

      integer nx,ny,nz,iv,nvars
      double precision :: ft(nx,ny,nz,nvars)

      double precision x,y,z,f(1,nvars)
      double precision xt(nx),yt(ny),zt(nz)
      double precision d1,d2,d3
!
!
      double precision  fh(1,8,nvars), delx(1), dely(1), delz(1), a1(1,nvars),&
          a2(1,nvars), a3(1,nvars), a4(1,nvars), a5(1,nvars), a6(1,nvars),&
          a7(1,nvars), a8(1,nvars)

      double precision dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

!
!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)
!
      dxi   = 1. / dx
      dyi   = 1. / dy
      dzi   = 1. / dz
!
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
!
      dxyzi = dxi * dyi * dzi

!------- loop over all points to be interpolated
!
      n = 1
!
!------- determine location in (equidistant!!!) table
!
         ix = 2 + INT( (x - xt(1) - 1.d-10) * dxi )
         iy = 2 + INT( (y - yt(1) - 1.d-10) * dyi )
         iz = 2 + INT( (z - zt(1) - 1.d-10) * dzi )
!
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
!
!         write(*,*) iy-1,iy,iy+1
!
!------- set-up auxiliary arrays for Lagrange interpolation
!
         delx(n) = xt(ix) - x
         dely(n) = yt(iy) - y
         delz(n) = zt(iz) - z
!
         do iv = 1, nvars
            fh(n,1,iv) = ft(ix  , iy  , iz, iv  )
            fh(n,2,iv) = ft(ix-1, iy  , iz, iv  )
            fh(n,3,iv) = ft(ix  , iy-1, iz, iv  )
            fh(n,4,iv) = ft(ix  , iy  , iz-1, iv)
            fh(n,5,iv) = ft(ix-1, iy-1, iz, iv  )
            fh(n,6,iv) = ft(ix-1, iy  , iz-1, iv)
            fh(n,7,iv) = ft(ix  , iy-1, iz-1, iv)
            fh(n,8,iv) = ft(ix-1, iy-1, iz-1, iv)
!
!------ set up coefficients of the interpolation polynomial and
!       evaluate function values
            !
            a1(n,iv) = fh(n,1,iv)
            a2(n,iv) = dxi   * ( fh(n,2,iv) - fh(n,1,iv) )
            a3(n,iv) = dyi   * ( fh(n,3,iv) - fh(n,1,iv) )
            a4(n,iv) = dzi   * ( fh(n,4,iv) - fh(n,1,iv) )
            a5(n,iv) = dxyi  * ( fh(n,5,iv) - fh(n,2,iv) - fh(n,3,iv) + fh(n,1,&
               iv) )
            a6(n,iv) = dxzi  * ( fh(n,6,iv) - fh(n,2,iv) - fh(n,4,iv) + fh(n,1,&
               iv) )
            a7(n,iv) = dyzi  * ( fh(n,7,iv) - fh(n,3,iv) - fh(n,4,iv) + fh(n,1,&
               iv) )
            a8(n,iv) = dxyzi * ( fh(n,8,iv) - fh(n,1,iv) + fh(n,2,iv) + fh(n,3,&
               iv) + fh(n,4,iv) - fh(n,5,iv) - fh(n,6,iv) - fh(n,7,iv) )
!
            f(n,iv)  = a1(n,iv) +  a2(n,iv) * delx(n)                         &
               +  a3(n,iv) * dely(n)                         +  a4(n,&
               iv) * delz(n)                         +  a5(n,&
               iv) * delx(n) * dely(n)               +  a6(n,&
               iv) * delx(n) * delz(n)               +  a7(n,&
               iv) * dely(n) * delz(n)               +  a8(n,&
               iv) * delx(n) * dely(n) * delz(n)
!

           enddo
!
    end SUBROUTINE

end module mod_eos_interpolation

