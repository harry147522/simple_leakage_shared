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
  subroutine intep3d(x, y, z, f, &
                  ft, &
                  nx, ny, nz, &
                  xt, yt, zt, &
                  ix_out, iy_out, iz_out, &
                  d1, d2, d3)
     implicit none
     integer, intent(in)   ::  nx, ny, nz
     double precision, intent(in)                   ::  x, y, z
     double precision, intent(out)                  ::  f
     integer, intent(out), optional                 ::  ix_out, iy_out, iz_out
     double precision, intent(out), optional        ::  d1, d2, d3
     double precision, intent(in)                   ::  xt(nx), yt(ny), zt(nz)
     double precision, intent(in)                   ::  ft(nx, ny, nz)

     integer                                        ::  ix, iy, iz
     integer                                        ::  pm1x, pm1y, pm1z
     double precision                               ::  dx, dy, dz
     double precision                               ::  dxi, dyi, dzi
     double precision                               ::  dxyi, dyzi, dxzi
     double precision                               ::  dxyzi

     double precision                               ::  delx, dely, delz
     double precision                               ::  fh(8), a(8)

!         ix = 2 + INT( (x(1) - xt(1) - 1.e-10) * dxi )
!         iy = 2 + INT( (y(1) - yt(1) - 1.e-10) * dyi )
!         iz = 2 + INT( (z(1) - zt(1) - 1.e-10) * dzi )
!!
!         ix = MAX( 2, MIN( ix, nx ) )
!         iy = MAX( 2, MIN( iy, ny ) )
!         iz = MAX( 2, MIN( iz, nz ) )

!     ix = minloc( dabs( xt(2:nx-1) - x ), dim = 1 )
!     iy = minloc( dabs( yt(2:ny-1) - y ), dim = 1 )
!     iz = minloc( dabs( zt(2:nz-1) - z ), dim = 1 )
!!!!!!!!!!!!! patrick version  
     ix = 1 +  minloc( dabs( xt(2:nx-1) - x ), dim = 1 )
     iy = 1 + minloc( dabs( yt(2:ny-1) - y ), dim = 1 )
     iz = 1 + minloc( dabs( zt(2:nz-1) - z ), dim = 1 )

     if(present(ix_out)) ix_out = ix
     if(present(iy_out)) iy_out = iy
     if(present(iz_out)) iz_out = iz


!      write(*,*) ix,iy,iz

     delx = x - xt(ix)
     dely = y - yt(iy)
     delz = z - zt(iz)
                                            
     if (delx > 0.0d0) then 
        pm1x = 1
     else
        pm1x = -1
     end if

     if (dely > 0.0d0) then 
        pm1y = 1
     else
        pm1y = -1
     end if

     if (delz > 0.0d0) then 
        pm1z = 1
     else
        pm1z = -1
     end if

     ! determine spacing parameters of (equidistant!!!) table

     dx    = ( xt(ix+pm1x)-xt(ix) )
     dy    = ( yt(iy+pm1y)-yt(iy) )
     dz    = ( zt(iz+pm1z)-zt(iz) )

     dxi   = 1. / dx
     dyi   = 1. / dy
     dzi   = 1. / dz

     dxyi  = dxi * dyi
     dxzi  = dxi * dzi
     dyzi  = dyi * dzi

     dxyzi = dxi * dyi * dzi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!! GR1D version
!      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
!      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
!      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)
!
!      dxi   = 1.d0 / dx
!      dyi   = 1.d0 / dy
!      dzi   = 1.d0 / dz
!
!      dxyi  = dxi * dyi
!      dxzi  = dxi * dzi
!      dyzi  = dyi * dzi
!
!      dxyzi = dxi * dyi * dzi
!
!!------- loop over all points to be interpolated
!
!!------- determine location in (equidistant!!!) table
!
!      ix = 2 + INT( (x - xt(1) - 1.d-10) * dxi )
!      iy = 2 + INT( (y - yt(1) - 1.d-10) * dyi )
!      iz = 2 + INT( (z - zt(1) - 1.d-10) * dzi )
!
!      ix = MAX( 2, MIN( ix, nx ) )
!      iy = MAX( 2, MIN( iy, ny ) )
!      iz = MAX( 2, MIN( iz, nz ) )
!
!
!     if(present(ix_out)) ix_out = ix
!     if(present(iy_out)) iy_out = iy
!     if(present(iz_out)) iz_out = iz
!
!         delx = xt(ix) - x
!         dely = yt(iy) - y
!         delz = zt(iz) - z
!
!!     if (delx > 0.0d0) then
!!        pm1x = 1
!!     else
!!        pm1x = -1
!!     end if
!!
!!     if (dely > 0.0d0) then
!!        pm1y = 1
!!     else
!!        pm1y = -1
!!     end if
!!
!!     if (delz > 0.0d0) then
!!        pm1z = 1
!!     else
!!        pm1z = -1
!!     end if

!write(*,*) pm1x,pm1y,pm1z
!! patrick's
     fh(1) = ft(ix  , iy  , iz  )                             
     fh(2) = ft(ix+pm1x, iy  , iz  )                             
     fh(3) = ft(ix  , iy+pm1y, iz  )                             
     fh(4) = ft(ix  , iy  , iz+pm1z)                             
     fh(5) = ft(ix+pm1x, iy+pm1y, iz  )                             
     fh(6) = ft(ix+pm1x, iy  , iz+pm1z)                             
     fh(7) = ft(ix  , iy+pm1y, iz+pm1z)                             
     fh(8) = ft(ix+pm1x, iy+pm1y, iz+pm1z)                             

     ! set up coefficients of the interpolation polynomial and 
     ! evaluate function values 
     a(1) = fh(1)                             
     a(2) = dxi   * ( fh(2) - fh(1) )       
     a(3) = dyi   * ( fh(3) - fh(1) )       
     a(4) = dzi   * ( fh(4) - fh(1) )       
     a(5) = dxyi  * ( fh(5) - fh(2) - fh(3) + fh(1) )
     a(6) = dxzi  * ( fh(6) - fh(2) - fh(4) + fh(1) )
     a(7) = dyzi  * ( fh(7) - fh(3) - fh(4) + fh(1) )
     a(8) = dxyzi * ( fh(8) - fh(1) + fh(2) + fh(3) + &
                      fh(4) - fh(5) - fh(6) - fh(7) )

     ! output values

     if(present(d1)) d1 = -a(2)
     if(present(d2)) d2 = -a(3)
     if(present(d3)) d3 = -a(4)
     f  = a(1) +  a(2) * delx                      &
               +  a(3) * dely                      &
               +  a(4) * delz                      &
               +  a(5) * delx * dely               &
               +  a(6) * delx * delz               &
               +  a(7) * dely * delz               &
               +  a(8) * delx * dely * delz     
!!!!!!
! gr1d
!     fh(1) = ft(ix  , iy  , iz  )                             
!     fh(2) = ft(ix-1, iy  , iz  )                             
!     fh(3) = ft(ix  , iy-1, iz  )                             
!     fh(4) = ft(ix  , iy  , iz-1)                             
!     fh(5) = ft(ix-1, iy-1, iz  )                             
!     fh(6) = ft(ix-1, iy  , iz-1)                             
!     fh(7) = ft(ix  , iy-1, iz-1)                             
!     fh(8) = ft(ix-1, iy-1, iz-1)                             
!
!     ! set up coefficients of the interpolation polynomial and 
!     ! evaluate function values 
!     a(1) = fh(1)                             
!     a(2) = dxi   * ( fh(2) - fh(1) )       
!     a(3) = dyi   * ( fh(3) - fh(1) )       
!     a(4) = dzi   * ( fh(4) - fh(1) )       
!     a(5) = dxyi  * ( fh(5) - fh(2) - fh(3) + fh(1) )
!     a(6) = dxzi  * ( fh(6) - fh(2) - fh(4) + fh(1) )
!     a(7) = dyzi  * ( fh(7) - fh(3) - fh(4) + fh(1) )
!     a(8) = dxyzi * ( fh(8) - fh(1) + fh(2) + fh(3) + &
!                      fh(4) - fh(5) - fh(6) - fh(7) )
!
!     ! output values
!
!     if(present(d1)) d1 = -a(2)
!     if(present(d2)) d2 = -a(3)
!     if(present(d3)) d3 = -a(4)
!     f  = a(1) +  a(2) * delx                      &
!               +  a(3) * dely                      &
!               +  a(4) * delz                      &
!               +  a(5) * delx * dely               &
!               +  a(6) * delx * delz               &
!               +  a(7) * dely * delz               &
!               +  a(8) * delx * dely * delz     

  end subroutine 

  subroutine intep3d_many(x, y, z, f, &
                  ft, &
                  nx, ny, nz, nvars, &
                  xt, yt, zt)
   !               ix_out, iy_out, iz_out)
     implicit none
     integer, intent(in)   ::  nx, ny, nz, nvars
     integer :: iv
     double precision, intent(in)                   ::  x, y, z
     double precision, intent(out)                  ::  f(nvars)
!     integer, intent(out), optional                 ::  ix_out, iy_out, iz_out
!     double precision, intent(out), optional        ::  d1, d2, d3
     double precision, intent(in)                   ::  xt(nx), yt(ny), zt(nz)
     double precision, intent(in)                   ::  ft(nx, ny, nz, nvars)

     integer                                        ::  ix, iy, iz
     integer                                        ::  pm1x, pm1y, pm1z
     double precision                               ::  dx, dy, dz
     double precision                               ::  dxi, dyi, dzi
     double precision                               ::  dxyi, dyzi, dxzi
     double precision                               ::  dxyzi

     double precision                               ::  delx, dely, delz
     double precision                               ::  fh(8,nvars), a(8,nvars)

! patrick version
     ix = 1+minloc( dabs( xt(2:nx-1) - x ), dim = 1 )
     iy = 1+minloc( dabs( yt(2:ny-1) - y ), dim = 1 )
     iz = 1+minloc( dabs( zt(2:nz-1) - z ), dim = 1 )

!     ix = minloc( dabs( xt(2:nx-1) - x ), dim = 1 )
!     iy = minloc( dabs( yt(2:ny-1) - y ), dim = 1 )
!     iz = minloc( dabs( zt(2:nz-1) - z ), dim = 1 )
!     if(present(ix_out)) ix_out = ix
!     if(present(iy_out)) iy_out = iy
!     if(present(iz_out)) iz_out = iz

     delx = x - xt(ix)
     dely = y - yt(iy)
     delz = z - zt(iz)

     if (delx > 0.0d0) then
        pm1x = 1
     else
        pm1x = -1
     end if

     if (dely > 0.0d0) then
        pm1y = 1
     else
        pm1y = -1
     end if

     if (delz > 0.0d0) then
        pm1z = 1
     else
        pm1z = -1
     end if

     ! determine spacing parameters of (equidistant!!!) table

     dx    = ( xt(ix+pm1x)-xt(ix) )
     dy    = ( yt(iy+pm1y)-yt(iy) )
     dz    = ( zt(iz+pm1z)-zt(iz) )

     dxi   = 1. / dx
     dyi   = 1. / dy
     dzi   = 1. / dz

     dxyi  = dxi * dyi
     dxzi  = dxi * dzi
     dyzi  = dyi * dzi

     dxyzi = dxi * dyi * dzi
!!!!!!!!!!!!!!!!
!!!!!GR1D version
!      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
!      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
!      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)
!!
!      dxi   = 1. / dx
!      dyi   = 1. / dy
!      dzi   = 1. / dz
!!
!      dxyi  = dxi * dyi
!      dxzi  = dxi * dzi
!      dyzi  = dyi * dzi
!!
!      dxyzi = dxi * dyi * dzi
!!
!!
!!------- loop over all points to be interpolated
!!
!!
!!------- determine location in (equidistant!!!) table
!!
!         ix = 2 + INT( (x - xt(1) - 1.e-10) * dxi )
!         iy = 2 + INT( (y - yt(1) - 1.e-10) * dyi )
!         iz = 2 + INT( (z - zt(1) - 1.e-10) * dzi )
!!
!         ix = MAX( 2, MIN( ix, nx ) )
!         iy = MAX( 2, MIN( iy, ny ) )
!         iz = MAX( 2, MIN( iz, nz ) )
!!
!!         write(*,*) iy-1,iy,iy+1
!!
!!------- set-up auxiliary arrays for Lagrange interpolation
!!
!         delx = xt(ix) - x
!         dely = yt(iy) - y
!         delz = zt(iz) - z
!
!!     if (delx > 0.0d0) then
!!        pm1x = 1
!!     else
!!        pm1x = -1
!!     end if
!!
!!     if (dely > 0.0d0) then
!!        pm1y = 1
!!     else
!!        pm1y = -1
!!     end if
!!
!!     if (delz > 0.0d0) then
!!        pm1z = 1
!!     else
!!        pm1z = -1
!!     end if


  do iv = 1, nvars
     fh(1,iv) = ft(ix  , iy  , iz  , iv)
     fh(2,iv) = ft(ix+pm1x, iy  , iz , iv )
     fh(3,iv) = ft(ix  , iy+pm1y, iz , iv )
     fh(4,iv) = ft(ix  , iy  , iz+pm1z, iv )
     fh(5,iv) = ft(ix+pm1x, iy+pm1y, iz , iv )
     fh(6,iv) = ft(ix+pm1x, iy  , iz+pm1z, iv)
     fh(7,iv) = ft(ix  , iy+pm1y, iz+pm1z, iv)
     fh(8,iv) = ft(ix+pm1x, iy+pm1y, iz+pm1z, iv)

     ! set up coefficients of the interpolation polynomial and
     ! evaluate function values
     a(1,iv) = fh(1,iv)
     a(2,iv) = dxi   * ( fh(2,iv) - fh(1,iv) )
     a(3,iv) = dyi   * ( fh(3,iv) - fh(1,iv) )
     a(4,iv) = dzi   * ( fh(4,iv) - fh(1,iv) )
     a(5,iv) = dxyi  * ( fh(5,iv) - fh(2,iv) - fh(3,iv) + fh(1,iv) )
     a(6,iv) = dxzi  * ( fh(6,iv) - fh(2,iv) - fh(4,iv) + fh(1,iv) )
     a(7,iv) = dyzi  * ( fh(7,iv) - fh(3,iv) - fh(4,iv) + fh(1,iv) )
     a(8,iv) = dxyzi * ( fh(8,iv) - fh(1,iv) + fh(2,iv) + fh(3,iv) + &
                      fh(4,iv) - fh(5,iv) - fh(6,iv) - fh(7,iv) )

     ! output values

!     if(present(d1)) d1 = -a(2)
!     if(present(d2)) d2 = -a(3)
!     if(present(d3)) d3 = -a(4)
     f(iv)  = a(1,iv) +  a(2,iv) * delx                      &
               +  a(3,iv) * dely                      &
               +  a(4,iv) * delz                      &
               +  a(5,iv) * delx * dely               &
               +  a(6,iv) * delx * delz               &
               +  a(7,iv) * dely * delz               &
               +  a(8,iv) * delx * dely * delz

!     fh(1,iv) = ft(ix  , iy  , iz  , iv)
!     fh(2,iv) = ft(ix-1, iy  , iz , iv )
!     fh(3,iv) = ft(ix  , iy-1, iz , iv )
!     fh(4,iv) = ft(ix  , iy  , iz-1, iv )
!     fh(5,iv) = ft(ix-1, iy-1, iz , iv )
!     fh(6,iv) = ft(ix-1, iy  , iz-1, iv)
!     fh(7,iv) = ft(ix  , iy-1, iz-1, iv)
!     fh(8,iv) = ft(ix-1, iy-1, iz-1, iv)
!
!     ! set up coefficients of the interpolation polynomial and
!     ! evaluate function values
!     a(1,iv) = fh(1,iv)
!     a(2,iv) = dxi   * ( fh(2,iv) - fh(1,iv) )
!     a(3,iv) = dyi   * ( fh(3,iv) - fh(1,iv) )
!     a(4,iv) = dzi   * ( fh(4,iv) - fh(1,iv) )
!     a(5,iv) = dxyi  * ( fh(5,iv) - fh(2,iv) - fh(3,iv) + fh(1,iv) )
!     a(6,iv) = dxzi  * ( fh(6,iv) - fh(2,iv) - fh(4,iv) + fh(1,iv) )
!     a(7,iv) = dyzi  * ( fh(7,iv) - fh(3,iv) - fh(4,iv) + fh(1,iv) )
!     a(8,iv) = dxyzi * ( fh(8,iv) - fh(1,iv) + fh(2,iv) + fh(3,iv) + &
!                      fh(4,iv) - fh(5,iv) - fh(6,iv) - fh(7,iv) )
!
!     ! output values
!
!!     if(present(d1)) d1 = -a(2)
!!     if(present(d2)) d2 = -a(3)
!!     if(present(d3)) d3 = -a(4)
!     f(iv)  = a(1,iv) +  a(2,iv) * delx                      &
!               +  a(3,iv) * dely                      &
!               +  a(4,iv) * delz                      &
!               +  a(5,iv) * delx * dely               &
!               +  a(6,iv) * delx * delz               &
!               +  a(7,iv) * dely * delz               &
!               +  a(8,iv) * delx * dely * delz
  enddo

  end subroutine intep3d_many


end module mod_eos_interpolation
