!> Module for eos
module mod_eos_tabulated
  use mod_eos_tabulated_parameters

  implicit none
  public

contains

  !> Read this module's parameters from a file
  subroutine eos_tabulated_read_params(files)
    use mod_eos
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /eos_tabulated_list/ eos_table_name, eos_precision

    do n = 1, size(files)
       ! Try to read in the namelists. They can be absent or in a different
       ! order, since we rewind before each read.
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, eos_tabulated_list, end=111)
111    close(unitpar)
    end do

  end subroutine eos_tabulated_read_params

  subroutine eos_tabulated_activate()
    use mod_global_parameters
    use mod_eos
    use mod_eos_readtable

    call eos_tabulated_read_params(par_files)
    eos_type = tabulated

    call readtable(eos_table_name)

    eos_get_pressure_one_grid         => tabulated_get_pressure_one_grid
    eos_get_eps_one_grid              => tabulated_get_eps_one_grid
    eos_get_cs2_one_grid              => tabulated_get_cs2_one_grid

    eos_get_temp_one_grid             => tabulated_get_temp_one_grid
    eos_eps_get_all_one_grid          => tabulated_eps_get_all_one_grid
    eos_temp_get_all_one_grid         => tabulated_temp_get_all_one_grid

    eos_get_eps_range                 => tabulated_get_eps_range

    tabulated_logeps_to_logtemp       => tabulated_logeps_to_logtemp_newton
!    tabulated_logeps_to_logtemp       => tabulated_logeps_to_logtemp_falsi

  end subroutine eos_tabulated_activate

  subroutine tabulated_get_pressure_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    use mod_eos_interpolation
    implicit none
    
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: temp, ye

    double precision             :: log_rho, log_temp, log_eps, leps

    leps = eps
! restrict lower bound eps of the eostable
!    if (leps < eos_epsmin .or. leps < small_eps) then
!      leps = max(small_eps, eos_epsmin)
!    endif


!    if (rho<small_rho_thr) then
    if (rho<eos_rhomin*1.0d0) then
      ! call atmo_get_pressure_one_grid(prs,rho,eps)
       prs = small_press
       return
    end if

    if (rho>eos_rhomax) stop "nuc_eos: rho > rhomax"
    if (rho<eos_rhomin) stop "nuc_eos: rho < rhomin"

    if (present(temp)) then
      if (temp>eos_tempmax) stop "nuc_eos: temp > tempmax"
      if (temp<eos_tempmin) then
        write(*,*) temp, "nuc_eos: temp < tempmin in get_prs"
        stop 
      endif
    endif
    if (present(ye)) then
      if (ye>eos_yemax) stop "nuc_eos: ye > yemax"
      if (ye<eos_yemin) stop "nuc_eos: ye < yemin"
    endif
    if (leps>eos_epsmax) stop "nuc_eos: eps > epsmax"
!    if (leps<eos_epsmin) stop "nuc_eos: eps < epsmin"
    if (leps<eos_epsmin) write(*,*) "nuc_eos: eps < epsmin"

    log_rho = log10(rho)

    if (present(temp))  log_temp = log10(temp)
   
       log_eps = log10(leps)

       call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    call intep3d(log_rho, log_temp, ye, &
                prs, eos_tables(:,:,:,i_logpress), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    prs = 10.0d0**prs

  end subroutine tabulated_get_pressure_one_grid

  subroutine tabulated_get_eps_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    use mod_eos_interpolation

    implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in), optional :: temp, ye
    double precision, intent(inout) :: eps

    double precision             :: log_rho, log_temp
!code test
!    if (rho<small_rho_thr) then
    if (rho<eos_rhomin*1.0d0) then
!       call atmo_get_eps_one_grid(prs,rho,eps)
       eps = small_eps
        write(*,*) "passed rho<eos_rhomin  in get_eps"
       return
    end if
!
    if (rho>eos_rhomax) stop "nuc_eos: rho > rhomax"
    if (rho<eos_rhomin) stop "nuc_eos: rho < rhomin"

    if (present(temp)) then
      if (temp>eos_tempmax) stop "nuc_eos: temp > tempmax"
      if (temp<eos_tempmin) then
        write(*,*) temp, "nuc_eos: temp < tempmin in get_eps"
        stop
      endif


    endif

    if (present(ye)) then
      if (ye>eos_yemax) stop "nuc_eos: ye > yemax"
      if (ye<eos_yemin) stop "nuc_eos: ye < yemin"
    endif
  
    log_rho = log10(rho)
    log_temp = log10(temp)

    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    eps = 10.0d0**eps
    if (eps < eos_epsmin .or. eps < small_eps) then
!       stop "passed low eps in get_eps_eos" 
       eps = max(small_eps, eos_epsmin)
    endif

  end subroutine tabulated_get_eps_one_grid

  subroutine tabulated_get_cs2_one_grid(cs2,rho,eps,temp,ye)

    use mod_eos
    use mod_eos_interpolation
    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: temp, ye

    double precision             :: log_rho, log_eps, log_temp, leps

    leps = eps
! restrict lower bound eps of the eostable
!    if (leps < eos_epsmin .or. leps < small_eps) then
!      leps = max(small_eps, eos_epsmin)
!    endif
    
!code test
!    if (rho<small_rho_thr) then
    if (rho<eos_rhomin*1.2d0) then
!       call atmo_get_cs2_one_grid(cs2,rho,eps)
        cs2 = 0.0d0
       return
    end if

    if (rho>eos_rhomax) stop "nuc_eos: rho > rhomax"
    if (rho<eos_rhomin) stop "nuc_eos: rho < rhomin"
    
    if (present(temp)) then
      if (temp>eos_tempmax) stop "nuc_eos: temp > tempmax"
      if (temp<eos_tempmin) then
        write(*,*) temp, "nuc_eos: temp < tempmin in get_cs2"
        stop
      endif

    endif

    if (present(ye)) then
      if (ye>eos_yemax) stop "nuc_eos: ye > yemax"
      if (ye<eos_yemin) stop "nuc_eos: ye < yemin"
    endif

    if (leps>eos_epsmax) stop "nuc_eos: eps > epsmax"
!    if (leps<eos_epsmin) stop "nuc_eos: eps < epsmin"
    if (leps<eos_epsmin) write(*,*) "nuc_eos: eps < epsmin"


    log_rho = log10(rho)

    if (present(temp)) log_temp = log10(temp)
       log_eps = log10(leps)
       call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    call intep3d(log_rho, log_temp, ye, &
                cs2, eos_tables(:,:,:,i_cs2), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

  end subroutine tabulated_get_cs2_one_grid

  subroutine tabulated_get_temp_one_grid(rho,eps,temp,ye)

    use mod_eos
    use mod_eos_interpolation

    implicit none
    
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps, ye
    double precision, intent(inout) :: temp

    double precision             :: log_rho, log_temp, log_eps, leps
    integer                      :: irho, itemp, iye
    ! for root finding
    integer                      :: it_root
    integer                      :: iter_max = 30
    integer                      :: side=0
    double precision             :: tmp
    double precision             :: fp, fm, f
    double precision             :: zp, zm, z

    leps = eps
! restrict lower bound eps of the eostable
!    if (leps < eos_epsmin .or. leps < small_eps) then
!      leps = max(small_eps, eos_epsmin)
!    endif


    if (rho<eos_rhomin*1.0d0) then
!    if (rho<small_rho_thr) then
       temp = eos_tempmin
       return
    end if

    if (rho>eos_rhomax) stop "nuc_eos: rho > rhomax"
    if (rho<eos_rhomin) stop "nuc_eos: rho < rhomin"
    
      if (temp>eos_tempmax) stop "nuc_eos: temp > tempmax"
      if (temp<eos_tempmin) then
        write(*,*) temp, "nuc_eos: temp < tempmin in get_temp"
        stop
      endif



      if (ye>eos_yemax) stop "nuc_eos: ye > yemax"
      if (ye<eos_yemin) stop "nuc_eos: ye < yemin"

    if (leps>eos_epsmax) stop "nuc_eos: eps > epsmax"
!    if (leps<eos_epsmin) stop "nuc_eos: eps < epsmin"
    if (leps<eos_epsmin) write(*,*) "nuc_eos: eps < epsmin"

    log_rho = log10(rho)
    log_eps = log10(leps)
    log_temp = log10(temp)

    call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    temp = 10.0d0**log_temp

  end subroutine tabulated_get_temp_one_grid


  subroutine tabulated_eps_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
                       dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)

    use mod_eos
    use mod_eos_interpolation
!    use mod_eos_tabulated_parameters
    implicit none

    double precision, intent(in) :: ye
    double precision, intent(inout) :: eps, rho
    double precision, intent(inout), optional :: ent, prs, temp, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

    double precision             :: log_rho, log_temp, log_eps

    double precision             :: ffx(nvars), lprs
    double precision, parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))   !this just like gr1d
    double precision, parameter :: idealgamma = 1.66666666666d0


    if (rho>eos_rhomax) stop "nuc_eos: rho > rhomax"

!    if (rho<eos_rhomin) stop "nuc_eos: rho <rhomin"
!fixme   -->  small_rho = atmo_rho
    if (rho<eos_rhomin*1.0d0) then 
      write(*,*) "nuc_eos: rho < rhomin"
      rho = small_rho
      lprs= idealK1*((small_rho/rho_gf)**(idealgamma))*press_gf
      eps = lprs / rho / ( idealgamma - 1.0d0 )
      if (present(prs)) prs = lprs
      if (present(cs2)) cs2 = idealgamma*lprs/rho
      if (present(ent)) ent = 4.0d0
        write(*,*) lprs, eps, rho, "small prs, eps ,rho  in get_all"
        return
    endif

      if (ye>eos_yemax) stop "nuc_eos: ye > yemax"
      if (ye<eos_yemin) then
        write(*,*) ye
        stop "nuc_eos: ye < yemin at finding all"
      endif


    if (eps>eos_epsmax) stop "nuc_eos: eps > epsmax"
!    if (eps<eos_epsmin) stop "nuc_eos: eps < epsmin"
    if (eps<eos_epsmin) write(*,*) "nuc_eos: eps < epsmin"
    
    if (present(temp)) log_temp = log10(temp)   !this is a must for initialguess
        
    log_rho = log10(rho)
    log_eps = log10(eps)

    call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    call intep3d_many(log_rho, log_temp, ye, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)!, irho, itemp, iye)

   if (present(temp)) then
      temp = 10**log_temp
        
      if (temp>eos_tempmax) stop "nuc_eos: temp > tempmax"

      if (temp<eos_tempmin) then
        write(*,*) temp, "nuc_eos: temp < tempmin in eps_get_all"
        stop
      endif

   endif


   if (present(prs)) prs = 10**ffx(1)
   if (present(ent)) ent  = ffx(3)
   if (present(munu))  munu = ffx(4)
   if (present(cs2)) then
         cs2  = ffx(5)
      if (cs2> 1.0d0) write(*,*) "nuc_eos:cs2>1.0d0" !stop "nuc_eos: cs2 > 1.0d0"
      if (cs2< 0.0d0) write(*,*) "nuc_eos:cs2<1.0d0" !stop "nuc_eos: cs2 < 0.0d0"
   endif
!  derivatives
   if (present(dedt)) dedt = ffx(6)
   if (present(dpdrhoe))  dpdrhoe = ffx(7)
   if (present(dpderho))  dpderho = ffx(8)
!  chemical potentials
   if (present(muhat))  muhat = ffx(9)
   if (present(mu_e))  mu_e = ffx(10)
   if (present(mu_p))  mu_p = ffx(11)
   if (present(mu_n))  mu_n = ffx(12)
!  compositions
   if (present(xa)) xa = ffx(13)
   if (present(xh)) xh = ffx(14)
   if (present(xn)) xn = ffx(15)
   if (present(xp)) xp = ffx(16)
   if (present(abar)) abar = ffx(17)
   if (present(zbar)) zbar = ffx(18)



  end subroutine tabulated_eps_get_all_one_grid


  subroutine tabulated_temp_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
                       dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)

    use mod_eos
    use mod_eos_interpolation
!    use mod_eos_tabulated_parameters
    implicit none

    double precision, intent(in) :: ye
    double precision, intent(inout) :: eps, rho, temp
    double precision, intent(inout), optional :: ent, prs, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

    double precision             :: log_rho, log_temp, log_eps

    double precision             :: ffx(nvars), lprs
    double precision, parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))   !this just like gr1d
    double precision, parameter :: idealgamma = 1.66666666666d0


    if (rho>eos_rhomax) stop "nuc_eos: rho > rhomax"

!    if (rho<eos_rhomin) stop "nuc_eos: rho <rhomin"
!fixme   -->  small_rho = atmo_rho
    if (rho<eos_rhomin*1.0d0) then 
      write(*,*) "nuc_eos: rho < rhomin"
      rho = small_rho
      lprs= idealK1*((small_rho/rho_gf)**(idealgamma))*press_gf
      eps = lprs / rho / ( idealgamma - 1.0d0 )
      if (present(prs)) prs = lprs
      if (present(cs2)) cs2 = idealgamma*lprs/rho
      if (present(ent)) ent = 4.0d0
        write(*,*) lprs, eps, rho, "small prs, eps ,rho  in get_all"
        return
    endif

      if (ye>eos_yemax) stop "nuc_eos: ye > yemax"
      if (ye<eos_yemin) then
        write(*,*) ye
        stop "nuc_eos: ye < yemin at finding all"
      endif

    
    if (temp>eos_tempmax) stop "nuc_eos: temp > tempmax"
      if (temp<eos_tempmin) then
        write(*,*) temp, "nuc_eos: temp < tempmin in temp_get_all"
        stop
      endif



    log_temp = log10(temp)   !this is a must for initialguess
        
    log_rho = log10(rho)
   ! log_eps = log10(eps)

   ! call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    call intep3d_many(log_rho, log_temp, ye, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)!, irho, itemp, iye)

        


   if (present(prs)) prs = 10**ffx(1)
       eps = 10**ffx(2)
   if (present(ent)) ent  = ffx(3)
   if (present(munu))  munu = ffx(4)
   if (present(cs2)) then
         cs2  = ffx(5)
      if (cs2> 1.0d0) write(*,*) "nuc_eos:cs2>1.0d0" !stop "nuc_eos: cs2 > 1.0d0"
      if (cs2< 0.0d0) write(*,*) "nuc_eos:cs2<1.0d0" !stop "nuc_eos: cs2 < 0.0d0"
   endif
!  derivatives
   if (present(dedt)) dedt = ffx(6)
   if (present(dpdrhoe))  dpdrhoe = ffx(7)
   if (present(dpderho))  dpderho = ffx(8)
!  chemical potentials
   if (present(muhat))  muhat = ffx(9)
   if (present(mu_e))  mu_e = ffx(10)
   if (present(mu_p))  mu_p = ffx(11)
   if (present(mu_n))  mu_n = ffx(12)
!  compositions
   if (present(xa)) xa = ffx(13)
   if (present(xh)) xh = ffx(14)
   if (present(xn)) xn = ffx(15)
   if (present(xp)) xp = ffx(16)
   if (present(abar)) abar = ffx(17)
   if (present(zbar)) zbar = ffx(18)



  end subroutine tabulated_temp_get_all_one_grid





  subroutine tabulated_logeps_to_logtemp_falsi(lr,epsin,lt0,y)

    use mod_eos
    use mod_eos_interpolation
!    use mod_rootfinding
!    use mod_eos_tabulated_parameters

    implicit none

    double precision :: lr, epsin, y
    double precision :: lt0

    integer                      :: irho, itemp, iye
    double precision             :: tmp
    double precision             :: logtemp_min, logtemp_max
    ! for root finding
    integer                      :: iter_max = 300
    logical :: root_is_found
    integer                      :: it_root
    integer                      :: side=0
    double precision             :: fp, fm, f
    double precision             :: zp, zm, z


    call intep3d(lr, lt0, y, &
                tmp, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, ix_out=irho, iy_out=itemp, iz_out=iye)

    if (abs(tmp-epsin) < eos_precision*abs(epsin) ) return

    fm = eos_tables(irho,1,iye,i_logenergy) - epsin
    fp = eos_tables(irho,ntemp,iye,i_logenergy) - epsin
    zm = logtemp_table(1)
    zp = logtemp_table(ntemp)
!    logtemp_min = logtemp_table(1)
!    logtemp_max = logtemp_table(ntemp)

    ! get eps from temperature
!    find_root_for_func => func_eps_of_temp
!    call rootfinding_illinois(log_temp, logtemp_min, logtemp_max, eos_precision, iter_max, root_is_found)

!stop
!    call rootfinding_falsi(log_temp, logtemp_min, logtemp_max, eos_precision, iter_max, root_is_found)
 !   if ( .not. root_is_found) then
 !      call mpistop("Fail to find the root in tabulated eos")
 !   endif

    do it_root = 1, iter_max
       z = (fm * zp - fp * zm) / (fm - fp)

       if (abs(zp-zm)  <  (eos_precision * abs(zp+zm)) ) exit

       call intep3d(lr, z, y, &
                f, eos_tables(:,:,:,i_logenergy), &
                 nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
       f = f - epsin
       !if (abs(f) < eos_precision*abs(log_eps) ) exit

       if ( (f * fp) > 0.0d0 ) then
          !f and fp have same sign, copy z to zp
          zp = z
          fp = f
          if (side == 1) fm = fm / 2.0d0
          side = 1
       else if ( (f * fm) > 0.0d0 ) then
          !f and fm have same sign, copy z to zm
          zm = z
          fm = f
          if (side == 2) fp = fp / 2.0d0
          side = 2
       else !it looks like zero
          exit
       end if
    end do

    if (it_root >= iter_max-1) stop "fail to find the root in tabulated eos"

    lt0 = z

!  yt = logtemp_table
!    contains
!       double precision function func_eps_of_temp(log_temp_in)
!         double precision, intent(in) :: log_temp_in
!         call intep3d(log_rho, log_temp_in, ye, &
!                  func_eps_of_temp, eos_tables(:,:,:,i_logenergy), &
!                   nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
!         func_eps_of_temp = func_eps_of_temp - log_eps
!       end function
  end subroutine tabulated_logeps_to_logtemp_falsi




  subroutine tabulated_logeps_to_logtemp_newton(lr,epsin,lt0,y)

    use mod_eos
    use mod_eos_interpolation
!    use mod_rootfinding
!    use mod_eos_tabulated_parameters

    implicit none

  real*8 lr,lt0,y,epsin
  real*8 eps,lt,ldt
  real*8 tol
  real*8 d1,d2,d3
  real*8 eps0,eps1,lt1

  real*8 ltn,ltmax,ltmin
  real*8 tinput,rfeps

  integer :: rl = 0
  integer itmax,i,keyerrt
  integer ii,jj,kk

  !write(*,*) epsin, lr

  keyerrt=0

  tol=1.0d-10 ! need to find energy to less than 1 in 10^-10
  itmax=20 ! use at most 20 iterations, then bomb

  lt=lt0
  lt1=lt

  eps0=epsin
  eps1=eps0


  ltmax=logtemp_table(ntemp)
  ltmin=logtemp_table(1)

  ! Note: We are using Ewald's Lagrangian interpolator here!

  !preconditioning 1: do we already have the right temperature?

    call intep3d(lr, lt, y, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)

  if (abs(eps-eps0).lt.tol*abs(eps0)) then
     return
  endif
  lt1=lt
  eps1=eps

!  write(*,"(i4,1P12E19.10)") 0,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0,d2
!  write(*,"(i4,1P12E19.10)") 0,lr,lt0,y,ltmin,ltmax

  do i=1,itmax
     !d2 is the derivative deps/dlogtemp;
     ldt = -(eps - eps0)/d2
   !write(*,*) d2,  "d2"
!     if(ldt.gt.0.d0) ltmin=lt
!     if(ldt.le.0.d0) ltmax=lt
     ltn = lt+ldt
     ltn = min(ltn,ltmax)
     ltn = max(ltn,ltmin)
     lt1=lt
     lt=ltn
     eps1=eps
!     write(*,"(i4,1P12E19.10)") i,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0,d2,ldt


    call intep3d(lr, lt, y, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)
!     call findthis(lr,lt,y,eps,epst,d1,d2,d3)
     if (abs(eps - eps0).lt.tol*abs(eps0)) then
!        write(*,"(1P12E19.10)") tol,abs(eps-eps0)/eps0
        exit
     endif
     !setup new d2

     ! if we are closer than 10^-2  to the
     ! root (eps-eps0)=0, we are switching to
     ! the secant method, since the table is rather coarse and the
     ! derivatives may be garbage.
     if(abs(eps-eps0).lt.1.0d-3*abs(eps0)) then
        d2 = (eps-eps1)/(lt-lt1)
     endif
!     if(i.ge.10) then
!       write(*,*) "EOS: Did not converge in findtemp!"
!       write(*,*) "rl,logrho,logtemp0,ye,lt,eps,eps0,abs(eps-eps0)/eps0"
!     if(i.gt.5) then
!       write(*,"(i4,1P10E22.14)") i,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0
!     endif
  enddo


 if(i.ge.itmax) then
   !write(*,*)  "i > itmax"
    keyerrt=667
    call bisection(lr,lt0,y,eps0,lt,eos_tables(:,:,:,i_logenergy),keyerrt,1)
    if(keyerrt.eq.667) then
       if(lt.ge.log10(t_max_hack)) then
          ! handling for too high temperatures
          lt = log10(t_max_hack)
          keyerrt=0
          goto 12
       else if(abs(lt-log10(t_max_hack))/log10(t_max_hack).lt.0.025d0) then
          lt0 = min(lt,log10(t_max_hack))
          keyerrt=0
          goto 12
       else
          ! total failure
          write(*,*) "EOS: Did not converge in findtemp!"
          write(*,*) "rl,logrho,logtemp0,ye,lt,eps,eps0,abs(eps-eps0)/eps0"
          write(*,"(i4,i4,1P10E19.10)") i,rl,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0
          write(*,*) "Tried calling bisection... didn't help... :-/"
          write(*,*) "Bisection error: ",keyerrt
       endif
    endif

    lt0=min(lt,log10(t_max_hack))
    return
 endif

12 continue

  lt0=min(lt,log10(t_max_hack))
    


  end subroutine tabulated_logeps_to_logtemp_newton

!  subroutine tabulated_logeps_to_logent(log_rho,log_temp0,ye,sin)
!  
!  ! This routine finds the new temperature based
!  ! on rho, Y_e, entropy, sin = input of entropy
!      use mod_eos
!      use mod_eos_interpolation
!  
!    implicit none
!  
!    double precision, intent(in) :: log_rho, sin, ye
!    double precision, intent(inout) :: log_temp0
!  
!    double precision s,log_temp,ldt
!    double precision tol
!    double precision d1,d2,d3
!    double precision s0,s1,lt1
!  
!    double precision ltn,ltmax,ltmin
!    double precision tinput
!  
!    integer :: irho, itemp, iye
!    integer :: rl = 0
!    integer itmax,i
!    integer ii,jj,kk
!  
!  !  keyerrt=0
!  
!    tol=10.0d-10 !rfeps ! need to find energy to less than 1 in 10^-10
!    itmax=20 ! use at most 20 iterations, then bomb
!  
!    log_temp=log_temp0
!    lt1=log_temp
!  
!    s0=sin
!    s1=s0
!  
!    ltmax=logtemp_table(ntemp)
!    ltmin=logtemp_table(1)
!  
!  
!    !preconditioning 1: do we already have the right temperature?
!      call intep3d(log_rho, log_temp, ye, & 
!                  s, eos_tables(:,:,:,i_entropy), &
!             nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, irho, itemp, iye)
!  
!  
!  !  call findthis(lr,lt,y,s,alltables(:,:,:,3),d1,d2,d3)
!  !
!    if (abs(s-s0).lt.tol*abs(s0)) then
!       return
!    endif
!    lt1=log_temp
!    s1=s
!  
!  
!    do i=1,itmax
!       !d2 is the derivative ds/dlogtemp;
!       ldt = -(s - s0)/d2
!       ltn = log_temp+ldt
!       ltn = min(ltn,ltmax)
!       ltn = max(ltn,ltmin)
!       lt1=log_temp
!       log_temp=ltn
!       s1=s
!  
!  
!      call intep3d(log_rho, log_temp, ye, & 
!                  s, eos_tables(:,:,:,i_entropy), &
!             nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, irho, itemp, iye)
!  
!  !     call findthis(lr,lt,y,s,alltables(:,:,:,3),d1,d2,d3)
!       if (abs(s - s0).lt.tol*abs(s0)) then
!         exit
!       endif
!       !setup new d2
!  
!       ! if we are closer than 10^-2  to the
!       ! root (eps-eps0)=0, we are switching to
!       ! the secant method, since the table is rather coarse and the
!       ! derivatives may be garbage.
!       if(abs(s-s0).lt.1.0d-3*abs(s0)) then
!          d2 = (s-s1)/(log_temp-lt1)
!       endif
!    enddo
!  
!  
!   if(i.ge.itmax) then
!  !    keyerrt=667
!  
!            write(*,*) "EOS: Did not converge in findtemp_entropy!"
!          stop
!  !   call bisection(lr,lt0,y,s0,lt,alltables(:,:,:,3),keyerrt,2)
!  !    if(keyerrt.eq.667) then
!  !          write(*,*) "EOS: Did not converge in findtemp_entropy!"
!  !          write(*,*) "rl,logrho,logtemp0,ye,lt,s,s0,abs(s-s0)/s0"
!  !          write(*,"(i4,i4,1P10E19.10)") i,rl,lr,lt0,y,lt,s,s0,abs(s-s0)/s0
!  !          write(*,*) "Tried calling bisection... didn't help... :-/"
!  !          write(*,*) "Bisection error: ",keyerrt
!  !    endif
!  
!      log_temp0=log_temp
!      return
!   endif
!  
!  
!    log_temp0=log_temp
!  
!  
!  subroutine tabulated_logeps_to_logent

subroutine bisection(lr,lt0,y,eps0,lt,bivar,keyerrt,keybisect)

    use mod_eos
    use mod_eos_interpolation

  implicit none

  double precision lr,lt0,y,eps0,lt
  integer keyerrt

  integer keybisect
  ! 1 -> operate on energy
  ! 2 -> operate on entropy

  !temporary vars
  double precision lt1,lt2,ltmin,ltmax
  double precision f1,f2,fmid,dlt,ltmid
  double precision d1,d2,d3,tol
  double precision f1a,f2a

  integer bcount,i,itmax,maxbcount

  double precision bivar(*)

  bcount = 0
  maxbcount = 80

  tol=1.d-9 ! need to find energy to less than 1 in 10^-9
  itmax=50

  ltmax=logtemp_table(ntemp)
  ltmin=logtemp_table(1)

  lt = lt0
  lt1 = dlog10(min(10.0d0**ltmax,1.10d0*(10.0d0**lt0)))
  lt2 = dlog10(max(10.0d0**ltmin,0.90d0*(10.0d0**lt0)))



    call intep3d(lr, lt1, y, &
                f1a, bivar, &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)
    call intep3d(lr, lt2, y, &
                f2a, bivar, &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)


  f1=f1a-eps0
  f2=f2a-eps0

  keyerrt=0


  do while(f1*f2.ge.0.0d0)
     bcount=bcount+1
     lt1=dlog10(min(10.0d0**ltmax,1.2d0*(10.0d0**lt1)))
     lt2=dlog10(max(10.0d0**ltmin,0.8d0*(10.0d0**lt2)))


    call intep3d(lr, lt1, y, &
                f1a, bivar, &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)

    call intep3d(lr, lt2, y, &
                f2a, bivar, &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)

     f1=f1a-eps0
     f2=f2a-eps0
     if(bcount.ge.maxbcount) then
        keyerrt=667
        return
     endif
  enddo

  if(f1.lt.0.0d0) then
     lt=lt1
     dlt=dlog10((10.0D0**lt2)-(10.0d0**lt1))
  else
     lt=lt2
     dlt=dlog10((10.0D0**lt1)-(10.0d0**lt2))
  endif

  do i=1,itmax
     dlt=dlog10((10.0d0**dlt)*0.5d0)
     ltmid=dlog10(10.0d0**lt+10.0d0**dlt)

    call intep3d(lr, ltmid, y, &
                f2a, bivar, &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table, d1=d1,d2=d2,d3=d3)


     fmid=f2a-eps0
     if(fmid.le.0.0d0) lt=ltmid
     if(abs(1.0d0-f2a/eps0).lt.tol) then
        lt=ltmid
        return
     endif
  enddo



end subroutine bisection

!!idealgasatmophere   get all
!subroutine tabulated_get_all_atmo(rho,eps,prs,cs2)
!
!  implicit none
!    double precision, intent(in)     :: rho
!    double precision, intent(in)     :: eps
!    double precision, intent(inout), optional  :: prs,cs2
!    double precision :: dpde, dpdrho
!
!!     eps = prs / rho / ( atmo_gamma - 1.0d0 )
!
!      if (present(prs)) prs = ( atmo_gamma - 1.0d0 ) * rho * eps
!
!      dpde = (atmo_gamma - 1.0d0 ) * rho
!      dpdrho = (atmo_gamma - 1.0d0 ) * eps
!
!      if (present(cs2)) then
!        cs2 = dpdrho+dpde*prs/rho**2
!        cs2 = cs2 / ( (1.0d0 + prs/rho) + eps )
!      endif
!
!
!end subroutine tabulated_get_all_atmo




  subroutine tabulated_get_eps_range(rho, eps_min, eps_max, ye)
    use mod_eos
    use mod_eos_interpolation
!    use mod_eos_tabulated_parameters

    implicit none

    double precision, intent(in) :: rho
    double precision, intent(in), optional :: ye
    double precision, intent(out) :: eps_max, eps_min
    double precision             :: log_rho, log_temp, eps, lye

!! in normal case
!    eps_max = eos_epsmax
!    eps_min = eos_epsmin

  ! find eps_min    
!    log_rho = log10(eos_rhomin)
!    log_temp = log10(eos_tempmin)
!    lye = eos_yemin
!
!    call intep3d(log_rho, log_temp, lye, &
!                eps, eos_tables(:,:,:,i_logenergy), &
!           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
   
!   eps_min = max(eos_epsmin, eps)
!code test
!    eps_min = eos_epsmin
    eps_min = small_eps



!  write(*,*) eos_epsmin, eps, "eos_epsmin, eps"
!    eps_min = max(eos_epsmin, eps) 

    eps = 0.0d0
  ! find eps_max
    log_rho = log10(eos_rhomax)
    log_temp = log10(eos_tempmax)
    lye = eos_yemax

    call intep3d(log_rho, log_temp, lye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

!  write(*,*) eos_epsmax, eps, "eos_epsmax, eps"
    eps_max = min(eos_epsmax, eps)
!write(*,*) eps_max, eps_min, "max min "    
!stop "at first time of eps_range"


  end subroutine tabulated_get_eps_range  

  subroutine set_eos_epsmin_max

    use mod_eos
    use mod_eos_interpolation

    implicit none

    double precision             :: log_rho, log_temp, eps, ye


    eps = 0
    log_rho = log10(eos_rhomin)
    log_temp = log10(eos_tempmin)
    ye = eos_yemin

    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

!code test    this method find epsmin  too big!
!    eos_epsmin = 10.0d0**eps  

    eps = 0
    log_rho = log10(eos_rhomax)
    log_temp = log10(eos_tempmax)
    ye = eos_yemax

    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    eos_epsmax = 10.0d0**eps
!if (mype == 0 ) then
  write(*,*) "################################################"
  write(6,*) " In set_eos_epsmin_max"
  write(6,*) "Done reading eos tables ,  all in code unit"
  write(*,*)  eos_yemin, eos_yemax, "eos_yemin  and max"
  write(*,*)  eos_tempmin, eos_tempmax, "eos_tempmin  and max"
  write(*,*)  eos_rhomin, eos_rhomax, "eos_rhomin  and max"
  write(*,*)  small_rho, small_rho_thr, "small_rho,  small_rho_thr"
eos_epsmin = small_eps
  write(*,*)  eos_epsmin, eos_epsmax, "eos_epsmin  and max"
  write(*,*)  small_eps, "small_eps" 
  write(*,*) "################################################"
  


!endif
!   stop "set_eos_epsmin max"
  end subroutine set_eos_epsmin_max

end module mod_eos_tabulated
