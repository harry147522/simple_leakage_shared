!-*-f90-*-
module mod_yeofrho

  implicit none
  private
!  logical :: do_ye_of_rho = .true.
!    logical :: do_yeofrhofit = .true. 
!    logical :: do_highcorrection = .false.
!  do_yeofrhofit
  double precision, parameter ::  yeofrho_rho1 = 3.0d7 !lower density cut
  double precision, parameter ::  yeofrho_rho2 = 2.0d13 !upper density cut
  double precision, parameter ::  yeofrho_ye1 = 0.50d0 !upper ye cut
  double precision, parameter ::  yeofrho_ye2 = 0.278d0 !lower ye cut
  double precision, parameter ::  yeofrho_yec = 0.035d0 !correction parameter

!  do_highcorrection use
  double precision, parameter ::  yeofrho_rhoH = 1.0d0
  double precision, parameter ::  yeofrho_yeH = 1.0d0 
!  read ye_profile
  character*128 :: yeprofile_name
 
!  type profile_ye_of_rho
  double precision,allocatable :: yeinlogrho(:)
  double precision,allocatable :: rhoinlogrho(:)
  double precision,allocatable :: logrhoinlogrho(:)
  double precision,allocatable :: slopeinlogrho(:)
!  end type profile_ye_of_rho

  double precision,save :: minrho
  double precision,save :: maxrho
  integer,save :: maxindex
  double precision,save :: yeofrho_logrho1
  double precision,save :: yeofrho_logrho2
!  double precision,save :: yeofrho_ye1
!  double precision,save :: yeofrho_ye2
!  double precision,save :: yeofrho_yec
  logical, save :: do_highcorrection
  double precision,save :: yeofrho_logrhoH
!  double precision,save :: yeofrho_yeH

!public :: adjust_ye
!public :: read_yeprofile
public :: interpolate_ye
public :: fit_ye
public :: init_yeofrho


contains

  subroutine init_yeofrho
     implicit none
  
   write(*,*) 'Using ye_of_rho'
   yeofrho_logrho1 = log10(yeofrho_rho1)  
   yeofrho_logrho2 = log10(yeofrho_rho2)

     if (do_highcorrection) then
         yeofrho_logrhoH = log10(yeofrho_rhoH)
     endif

  end subroutine init_yeofrho


!  subroutine adjust_ye
!    
!!    use gmunu_module
!    use eos    
!    implicit none
!  
!    double precision out_ye, delta_ye, delta_S, out_munu
!    integer eosflag, keytemp, keyerr, i,j,k
!    double precision eosdummy(15)
!    
!    do k=1,n3
!    do j=1,n2
!    do i=1,n1
!       nuchem(i,j,k) = 0.0d0
!       if (rho(i,j,k).gt.1.0d6*rho_gf) then
!          if((.not.bounce).or.(i.gt.ishock(1))) then
!             if (do_yeofrhofit) then
!                call fit_ye(rho(i,j,k)/rho_gf,out_ye)
!             else
!                call interpolate_ye(rho(i,j,k)/rho_gf,out_ye)
!             endif
!             delta_ye = min(0.0d0, out_ye-ye(i,j,k))
!#if (EOS_KEY == 3)
!!             if (eoskey.eq.3) then
!                eosflag = 9 !get mu_nu
!                keytemp = 0
!
!                call eos(out_munu,rho(i,j,k),temp(i,j,k),ye(i,j,k),eps(i,j,k), &             !no need temp ??
!                          keytemp,keyerr,eosflag,eos_rf_prec)
!
!
!
!                nuchem(i,j,k) = out_munu
!                
!                if (out_munu.lt.10.0d0) then
!                   delta_S = 0.0d0
!                elseif (rho(i,j,k).lt.2.0d12*rho_gf) then
!                   delta_S = -delta_ye*(out_munu-10.0d0)/temp(i,j,k)
!                else
!                   delta_S = 0.0d0
!                endif
!
!                !this updates temp using constant entropy, then finds new energy 
!                !associated with that new temp, the one call does both.                
!                keytemp = 2
!                ent(i,j,k) =  ent(i,j,k) + delta_S
!                ye(i,j,k) = ye(i,j,k) + delta_ye
!
!!                write(*,*) "pass good way ye)ofrho.f90 !!!!!!!!!!!!!"
!                keyerr = 0
!                call eos_full(rho(i,j,k),temp(i,j,k),ye(i,j,k),eps(i,j,k),press(i,j,k),eosdummy(14), &
!                              ent(i,j,k), &
!                              cs2(i,j,k), &
!                              eosdummy(2),&
!                              eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
!                              eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
!                              eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(15), &
!                              keytemp,keyerr,eos_rf_prec)
!
!                if (keyerr.ne.0) then
!                   stop "problem in ye_of_rho: eos"
!                endif
!!             else
!!                ye(i) = ye(i) + delta_ye
!!             endif
!#endif
!          endif
!       endif
!    enddo
!    enddo
!    enddo
!
!  end subroutine adjust_ye

  subroutine read_yeprofile

    !in CGS
    implicit none
    
    character(len=100) filename
    integer plen, i, j
    logical proceed
    double precision,allocatable :: profileye(:)
    double precision,allocatable :: profilerho(:)
    double precision tempvar, tempvar2, tempvar3
    double precision slope

    filename = trim(adjustl(yeprofile_name))
    open(unit=666,file=trim(adjustl(filename)))    
    
    read(666,"(I6)") plen
    allocate(profileye(plen))
    allocate(profilerho(plen))

    read(666,"(E21.5,E20.5,E12.5)") profilerho(plen), profileye(plen), tempvar
    maxrho = profilerho(plen)
    minrho = profilerho(plen)
    do i=1,plen-1
       read(666,"(E21.5,E20.5,E12.5)") profilerho(plen-i), profileye(plen-i), tempvar
       if (profilerho(plen-i).gt.maxrho) maxrho=profilerho(plen-i)
       if (profilerho(plen-i).lt.minrho) minrho=profilerho(plen-i)
    enddo

    !limits of lookup table
    if (minrho.le.1.0d0) stop "rho very small, adjust method"
    maxindex = (int(log10(maxrho)*100.0d0)+1)
    !allocate lookup table
    allocate(yeinlogrho(maxindex))
    allocate(rhoinlogrho(maxindex))
    allocate(logrhoinlogrho(maxindex))
    allocate(slopeinlogrho(maxindex))

    !density as a function of index is (10.0d0)**(real(index)/100.0)
    do i=1,maxindex
       rhoinlogrho(i) = (10.0d0)**(real(i)/100.0d0)
       logrhoinlogrho(i) = real(i)/100.0d0
    enddo

    !monotonize the ye
    do i=2,plen
       if (profileye(i).gt.profileye(i-1)) then
          profileye(i) = profileye(i-1)
       endif
    enddo

    !pick first lookup index, interate through profile till we straddle the density we want
    yeinlogrho(1) = profileye(1)
    write(*,*) maxindex,maxrho,rhoinlogrho(1333),rhoinlogrho(1332),profilerho(80)
    do i=2,maxindex
       j=1

       proceed = .true.
       if (rhoinlogrho(i).gt.maxrho) then
          proceed = .false.
          yeinlogrho(i) = yeinlogrho(i-1)
       endif
       do while ((j.lt.(plen+1)).and.proceed)
          if (profilerho(j).gt.rhoinlogrho(i)) then
             proceed = .false.
             !set ye it interpolated value between next profile rho and last one
             if (j.eq.1) then
                !special case, can't interpolate, just assume ye does not change
                yeinlogrho(i) = yeinlogrho(i-1)
             else 
                !here we interpolate
                slope = (profileye(j)-profileye(j-1))/(&
                     (log10(profilerho(j))-log10(profilerho(j-1))))
                yeinlogrho(i) = slope*(logrhoinlogrho(i)-log10(profilerho(j)))+profileye(j)
             endif
          endif
          j=j+1
       enddo
       if (j.eq.(plen+1).and.proceed) then
          write(*,*) "Error",i
          stop
       endif
    enddo
    
    do i=1,maxindex-1
       slopeinlogrho(i) = (yeinlogrho(i+1)-yeinlogrho(i))/(logrhoinlogrho(i+1)-logrhoinlogrho(i))
    enddo

  end subroutine read_yeprofile

  subroutine interpolate_ye(incomingrho, outgoingye)

    ! in CGS
!    use gmunu_module
    implicit none
    
    double precision :: incomingrho,outgoingye
    integer :: index
    double precision :: logincomingrho

    if (incomingrho.lt.1.03d0) then
       stop "Figure something else out"
    end if

    logincomingrho = log10(incomingrho)
    index = int(logincomingrho*100.0d0)


    if (index.gt.maxindex-1) then
       outgoingye = yeinlogrho(maxindex)
    else
       outgoingye = slopeinlogrho(index)*(logincomingrho-logrhoinlogrho(index))+yeinlogrho(index)
    endif

  end subroutine interpolate_ye

  subroutine fit_ye(incomingrho, outgoingye)

    ! in CGS
    implicit none
    
    double precision incomingrho,outgoingye
    double precision logincomingrho
    double precision slope
    
    double precision x,absx
!write(*,*)  "correct inside fit_ye"
    if (incomingrho.lt.1.03d0) then
       stop "Figure something else out"
    end if

    logincomingrho = log10(incomingrho)

    if (logincomingrho.gt.yeofrho_logrho2 .and. do_highcorrection) then

       slope = (yeofrho_yeH - yeofrho_ye2) / (yeofrho_logrhoH - yeofrho_logrho2)

       outgoingye = yeofrho_ye2 + slope*(logincomingrho - yeofrho_logrho2)

    else

       x = max(-1.0d0,min(1.0d0,(2.0d0*logincomingrho - yeofrho_logrho2 - yeofrho_logrho1) &
            /(yeofrho_logrho2-yeofrho_logrho1)))

       absx = abs(x)
    
       outgoingye = 0.5d0*(yeofrho_ye2+yeofrho_ye1) + x/2.0d0*(yeofrho_ye2-yeofrho_ye1) &
            + yeofrho_yec*(1.0d0-absx+4.0d0*absx*(absx-0.5d0)*(absx-1.0d0))

    end if
!  write(*,*) yeofrho_logrho1, "yeofrho_logrho1"
!  write(*,*) yeofrho_logrho2
!  write(*,*) yeofrho_yec
 ! write(*,*) incomingrho, outgoingye, "incoming   outgoing" 


  end subroutine fit_ye

end module mod_yeofrho
    
