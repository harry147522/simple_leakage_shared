!   assumed spherical symmetric:  ns_lcoation
!                                 heat_rad (absorption
!  single core to do leakage -->   n1 = nM1 + ghostzones1 * 2

      !optical depth included in cons
!> (only in CFC now)
module mod_grhd_ccsn_leakage_simple_leakage

!  single core to do leakage -->   n1 = nM1 + ghostzones1 * 2

!                         v, v1,v2,v3,time,&          !time
!                         nt,t_bounce
  use mod_eos
  use mod_global_parameters
  use mod_physics
!  use mod_grhd_ccsn_leakage_phys
!  use mod_usr, only: bounce, t_bounce, shock_radius, pns_radius 
  implicit none
  private

      !constants
      integer, parameter :: ghostzones^D = 4

      double precision, parameter :: me_mev = 0.510998910d0 !mass of electron in MeV
      double precision, parameter :: sigma_0 = 1.76d-44 ! in units of cm^2
      double precision, parameter :: alpha = 1.25d0 ! dimension(ixI^S)less
      double precision, parameter :: Qnp = 1.293333d0 !neutron proton mass difference in MeV
      double precision, parameter :: hc_mevcm = 1.23984172d-10 !hc in units of MeV*cm
      double precision, parameter :: rho_min_leak = 1.0d4 !in g/cm^3
      double precision, parameter :: Cv = 0.5d0 + 2.0d0*0.23d0 ! dimensionless
      double precision, parameter :: Ca = 0.5d0 !dimensionless
      double precision, parameter :: gamma_0 = 5.565d-2 ! dimensionless
      double precision, parameter :: fsc = 1.0d0/137.036d0 ! fine structure constant, dimensionless
      double precision, parameter :: pi = 3.14159265358979d0

!!!!!!!!!!!!terminology!!!!!!!!!!!!
!
! neutrino species ---- electron type (e)     - 1
!                       antielectron type (a) - 2
!                       all others (x)        - 3
!
! matter species   ---- neutrons (n)          - 1
!                       protones (p)          - 2
!                       heavy nuclei (h)      - 3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision,allocatable  :: lprim(:^D&, :)   ! cgs primvarialbes
      double precision,allocatable  :: lleak_tau(:^D&, :)   ! read from cons(:,leak_tau(:)) 

!for mapping
      double precision,allocatable  :: lprim_local(:^D&, :)   ! cgs primvarialbes
      double precision,allocatable  :: lleak_tau_local(:^D&, :)   ! read from cons(:,leak_tau(:)) 

      double precision,allocatable  :: lx(:,:)    ! cgs coordinate local     ! only  depend 1:ndim to see r,theta,phi
      double precision,allocatable  :: lxi(:,:)   ! cgs cell face coordinate local

      double precision,allocatable  :: lx_map(:,:)    ! cgs coordinate local  for mapping, as length_gf will chane a little bit

      double precision,allocatable  :: lfac(:^D&)    ! Lorentz factor
      double precision,allocatable  :: gamma_ij(:^D&,:,:)
      double precision,allocatable  :: lcoolingsource(:^D&,:)   ! ndir are mom, ndir + 1 : tau, ndir +2 : Dye

      !number emission rate used in determining change in ye
      double precision,allocatable  :: R_eff(:^D&,:) !effective number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
      double precision,allocatable  :: R_loc(:^D&,:) !local number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3 
      double precision,allocatable  :: R_diff(:^D&,:) !diffusive number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3 
      double precision,allocatable  :: R_tot(:^D&)

      double precision,allocatable  :: Q_eff(:^D&,:) !effective energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
      !energy emission rate used in determining change in energy
      double precision,allocatable  :: Q_loc(:^D&,:) !local energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3 
      double precision,allocatable  :: Q_diff(:^D&,:) !diffusive energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3 
      double precision,allocatable  :: Q_tot(:^D&)

      double precision,allocatable  :: chi_single(:^D&,:) !optical depth with E^2 factored out,indices <radial zones:neutrino species> - 1 / MeV^2
      double precision,allocatable  :: zeta_single(:^D&,:) !Eq. A21, mean free path with E^2 factored out, indices <radial zones:neutrino species> - 1 / MeV^2 / cm

      double precision,allocatable  :: kappa_tilde_nu_scat(:^D&,:,:) !opacity without energy dependence for neutrino-matter scattering
                                                        !indices <radial zones:neutrino species:matter species> - 1 / MeV^2 / cm
      double precision,allocatable  :: kappa_tilde_nu_abs(:^D&,:,:) !opacity without energy dependence for neutrino absorption
                                                        !indices <radial zones:neutrino species:matter species> - 1 / MeV^2 / cm
   
      double precision,allocatable  :: lvol(:^D&) !cell volume in cm^3, indices <radial position>
      double precision,allocatable  :: lIarea(:^D&) !inverse cell area in cm^-2, indices <radial position>
      double precision,allocatable  :: lmass(:^D&) !cell mass in g, indices <radial position>
    
      !degeneracy parameters
         ! 1: electrons  2: p  3: n
      double precision,allocatable  :: eta_nucleons(:^D&,:) !e degeneracy, indices <radial position>, including rest mass - dimension(:^D&)less

      double precision,allocatable  :: eta_hat(:^D&) !n degeneracy - p dengeneracy - Qnp/temp, removes mass difference, indices <radial position> - dimension(:^D&)less
         ! 1: nue  2: nua  3: nux
      double precision,allocatable  :: eta_nu(:^D&,:) !nu_e degeneracy, from eta_e - eta_n + eta_p, indices <radial position> - dimension(:^D&)less

      double precision,allocatable  :: lepton_blocking(:^D&,:) !blocking terms for electrons and positrons, indices <radial position,lepton type> - dimensionless
      
      !mass fractions
      ! 1. xxn   2. xxp   3.  xxa   4. xxh   5. xabar   6. xzbar
      double precision,allocatable  :: mass_fraction(:^D&,:) 
    
      double precision,allocatable  :: eta_pn(:^D&) !proton number density corrected for neutron final state blocking, indices <radial position> - dimension(:^D&)less
      double precision,allocatable  :: eta_np(:^D&) !neutron number density corrected for proton final state blocking, indices <radial position> - dimension(:^D&)less


      !neutrino sphere
      integer          :: ns_location^D(3) !neutrino sphere location, indices <neutrino species> - radial index
      double precision :: ns_energy(3) !neutrino sphere energies, indices <neutrino species> - in MeV
      logical          :: have_ns_energies !flag that if true uses old energies as starting point, otherwise 18.0MeV is used

      double precision :: heat_erms(3) !root mean squared energy at neutrinosphere (F5/F3), indices <neutrino species> in MeV
      double precision :: heat_em(3) !mean energy at neutrinosphere (F4/F3), indices <neutrino species> in MeV
      double precision :: fermi_eta_3(2)     ! = get_fermi_integral(3, eta_nu)) 1: nue, 2: nua
      double precision :: fermi_eta_4(2)     ! = get_fermi_integral(4, eta_nu))
      double precision :: fermi_eta_5(2)     ! = get_fermi_integral(5, eta_nu))

      double precision, public :: t_dump = 0.0d0
       
      logical :: do_heating = .false.
      double precision :: heat_fac
      logical :: do_NNBrem = .false. 
      logical, public :: first_iteration
      logical :: output = .true.
 !Public methods
      logical, public  :: bounce = .true.   !false
      double precision, public :: shock_radius, pns_radius   ! in cm 
      double precision, public :: t_bounce
  public :: grhd_simple_leakage_init
  public :: grhd_simple_leakage 
!  public :: grhd_simple_leakage_add_source 
  public :: grhd_simple_leakage_activate
contains

  subroutine grhd_simple_leakage_activate()
     call grhd_simple_leakage_init()
  end subroutine grhd_simple_leakage_activate


  subroutine grhd_simple_leakage_read_params(files)
    use mod_global_parameters
    implicit none
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grhd_simple_leakage_list/ do_heating, heat_fac, do_NNBrem

    write(*,*) "Simple leakage is activated !"

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grhd_simple_leakage_list, end=111)
111    close(unitpar)
    end do

  !first_iteration = .true.
   !  bounce = .true.

  end subroutine grhd_simple_leakage_read_params


  !> Read this module's parameters from a file
  subroutine grhd_simple_leakage_init()
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    implicit none

!   write(*,*) "u passed a wrong route fku"
!   stop
   call grhd_simple_leakage_read_params(par_files)

   if (coordinate /= spherical) then 
        stop  "only spherical coordinate can use simple leakage"
   endif

   usr_before_addsource => grhd_simple_leakage_nompi
!   usr_source => grhd_simple_leakage_add_source 

  end subroutine grhd_simple_leakage_init
!######################################################################

  subroutine grhd_simple_leakage_nompi(iit,qt)
    use mod_global_parameters
    use mod_geometry
    implicit none

    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                         :: idir
    integer :: iigrid, igrid, ix^D

!    dx^D=rnode(rpdx^D_,igrid);
!    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
!    saveigrid=igrid

    if (.not. bounce) then
        return
    endif


    do iigrid=1,igridstail; igrid=igrids(iigrid);
    {do ix^D = ixM^LL \}

       do idir = 1, ndir+2
       ps(igrid)%cons(ix^D, coolingsource(idir)) = 0.0d0
       enddo
    {enddo^D&\}

    enddo
 
    if (bounce) then
!        write(*,*) "leakage"
!        write(*,*) "wrong path b4 leakage"
!    single core to do it
!code test
       call grhd_simple_leakage   !inside cons(leak_tau(:)) changed
    endif
  end subroutine grhd_simple_leakage_nompi
!######################

  subroutine grhd_simple_leakage_add_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,primCT,qt,cons,x)
    use mod_global_parameters
    use mod_geometry
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qtC, qdt, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: primCT(ixI^S, 1:nprim), qt
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
!    logical, intent(in)             :: qsourcesplit
!    logical, intent(inout)          :: active !< Needs to be set to true when active
    integer                         :: idir,jdir,kdir

   write(*,*) "passed add_source_leakage"
stop

   write(*,*) cons(5,coolingsource(2)), x(5,r_), qdt
!stop
!write(*,*) "passed add_source"
!for code test
    if (bounce) then
       {^NOONED
       do idir =2, ndir
       ! 1D core treatment
          where ( (oneDcore).and.(x(ixO^S,r_) <= r_core) )
            cons(ixO^S, coolingsource(idir)) = 0.0d0
          end where
       end do
       }
!       write(*,*) "wrong path" 
       cons(ixO^S, tau_) = cons(ixO^S, tau_) + qdt * cons(ixO^S, coolingsource(ndir+1))
       cons(ixO^S, ye_con_) = cons(ixO^S, ye_con_) + qdt * cons(ixO^S, coolingsource(ndir+2))
       do idir =1, ndir
         cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) + qdt * cons(ixO^S, coolingsource(idir))
       end do
    endif

  end subroutine grhd_simple_leakage_add_source




!  subroutine grhd_simple_leakage(ixO^L, ixI^L, prim, cons, x, dx^D, coolingsource)
  subroutine grhd_simple_leakage
    use mod_global_parameters
    use mod_geometry
    !use mod_eos
      implicit none

      integer :: idir, ix^D, ix_next^D
      integer :: n^D,  ixI^L, ixO^L

    integer :: iigrid, igrid
    
    !  logical,save :: have_old_tau

      character*1024 :: filename
      double precision :: xmu_e, xmu_n, xmu_p
!######################################################################


!   write(*,*)  ps(1)%cons(5, mom(1)), "cons b4 simple leakage "
!   write(*,*)  ps(1)%prim(5, temp_), "temp b4 simple leakage "
    !finding total array size 
   !i can use ixG^LL%prim(ixG^T) 
!     write(*,*) ps(5)%cons(:,tau_), "b4 inverse mamp"

    call finding_array_size(n^D)

    call initialize_arrays(n^D)

    call mapping_to_single_core(n^D)

   ! have_old_tau = .false

   
   ixImin^D = 1
   ixImax^D = n^D
   ixOmin^D = 1+ghostzones^D
   ixOmax^D = n^D-ghostzones^D
!write(*,*) ixImin1, ixImax1, ixOmin1, ixOmax1
   call finding_intermediate_variables(ixI^L, ixO^L)


   ! only need change radius unit
    lxi(:, r_) = lxi(:, r_)/length_gf
    lx(:, r_) = lx(:, r_)/length_gf

!   ! gamma_ij has lx(:,1) needed to change unit 
!    {do ix^D = ixO^LIM^D \}
!       gamma_ij(ix^D,2,2) = gamma_ij(ix^D,2,2)/length_gf/length_gf 
!       gamma_ij(ix^D,3,3) = gamma_ij(ix^D,3,3)/length_gf/length_gf 
!    {enddo^D&\} ! end loop for every points

    ! copy over density, energy, temperature, ye
    ! convert units

!  NOT calculate ghostzone
    {do ix^D = ixO^LIM^D \}
       lvol(ix^D) = 4.0d0*pi/3.0d0* (lxi(ix1+1, r_)**3 - lxi(ix1, r_)**3 )  ! using inner surface x1    !approxed for 1D 2D 3D 
    {enddo^D&\} ! end loop for every points

    {do ix^D = ixO^LIM^D \}
       lIarea(ix^D) = 1.0d0/(4.0d0*pi*lx(ix1, r_)**2)                ! using cell center   1d 2d approximately the same
    {enddo^D&\} ! end loop for every points

!    ldt = dt/time_gf


!   do k=ghosts3+1,n3-ghosts3
!     do j=ghosts2+1,n2-ghosts2
!        do i=ghosts1+1,n1-ghosts1
!!#if (DIMENSIONS == 1)
!       lmass(i,j,k) = lvol(i,j,k)*lprim(ix^D, rho_)      ! change in 2d
!!#else
!!
!!#endif
!    enddo
!   enddo
!  enddo

!#############             EOS
  xmu_e = 0.0d0
  xmu_p = 0.0d0
  xmu_n = 0.0d0
! input is in code unit EOS


    {do ix^D = ixI^LIM^D \}
 !  {do ix^D = ixO^LIM^D \}

       !call eos to get quantities we need
        call eos_temp_get_all_one_grid(rho = lprim(ix^D, rho_), eps = lprim(ix^D, eps_),  temp = lprim(ix^D, temp_),&
                                ye = lprim(ix^D, ye_), &                                
                                xa = mass_fraction(ix^D,3), xh = mass_fraction(ix^D,4),&
                                xn = mass_fraction(ix^D,1),xp = mass_fraction(ix^D,2),&
                                abar = mass_fraction(ix^D,5), zbar = mass_fraction(ix^D,6), &
                                mu_e = xmu_e, mu_n = xmu_n, mu_p = xmu_p)

       ! in our EOS the rest mass difference is in the chemical potentials of
       ! the neucleons

       eta_nucleons(ix^D,1) = xmu_e/lprim(ix^D, temp_)
       eta_nucleons(ix^D,2) = xmu_p/lprim(ix^D, temp_)
       eta_nucleons(ix^D,3) = xmu_n/lprim(ix^D, temp_)

       eta_hat(ix^D) = eta_nucleons(ix^D,3)-eta_nucleons(ix^D,2) - Qnp/lprim(ix^D, temp_)

       eta_nu(ix^D,1) = eta_nucleons(ix^D,1) - eta_nucleons(ix^D,3) + eta_nucleons(ix^D,2) !fully includes effects of rest masses
       eta_nu(ix^D,2) = -eta_nu(ix^D,1)
       eta_nu(ix^D,3) = 0.0d0

    {enddo^D&\} ! end loop for every points

    {do ix^D = ixI^LIM^D \}
!need map ghostzones as well
    lprim(ix^D,rho_)  = lprim(ix^D, rho_)/rho_gf
    lprim(ix^D,eps_)  = lprim(ix^D, eps_)/eps_gf
    lprim(ix^D,cs2_)  = lprim(ix^D, cs2_) * clight**2
    lprim(ix^D,press_)  = lprim(ix^D, press_)/press_gf
    {enddo^D&\} ! end loop for every points


 !  write(*,*) lprim(5,temp_),lprim(50,temp_),lprim(150,temp_), "temp"

!without this, so god damn wrong       
! 2d need change
    ix1 = 1
        do while(lprim(ix1,rho_).gt.rho_min_leak.and.ix1.lt. ixOmax1)  
                ix1=ix1+1
        enddo

! later we should all use ixO not ixI
!   dont add ghostzone1
    ixOmax1 = ix1  !+ghostzones1


!passed  they are different
!write(*,*) shape(lprim(1:ixOmax1,rho_)), shape(lprim(:,rho_))

!#############s
    !find tau and interpolate degeneracy factors
    call find_ruf_tau(ixI^L, ixO^L) 

!write(*,*) lleak_tau(:,1), shape(lleak_tau(:,1)) 
    !find Q_diff and R_diff
    call compute_diffusion(ixI^L, ixO^L)

    !find Q_loc and R_loc
    call compute_emission(ixI^L, ixO^L)

    !compute the rate of energy and ye loss/gain in each cell for RK
    call compute_update(ixI^L, ixO^L)

!stop "b4 inverse_mapping"
    call inverse_mapping(n^D)
       
    call deallocate_arrays

   
!    do iigrid=1,igridstail; igrid=igrids(iigrid);
!     write(*,*) ps(igrid)%prim(:,rho_)
!     write(*,*) "####################################################"
!     write(*,*) ps(5)%cons(:,tau_), "after inverse mamp"
!    enddo
!stop
!   write(*,*)  ps(1)%cons(5, mom(1)), "cons after simple leakage "
!   write(*,*)  ps(1)%prim(5, temp_), "temp after simple leakage "
  end subroutine grhd_simple_leakage
!######################################################################
  subroutine compute_update(ixI^L, ixO^L)
    use mod_global_parameters
    implicit none

    integer, intent(in) :: ixI^L, ixO^L
    double precision :: heat_rad(ixI^S,3),heat_net(ixI^S)
    double precision :: heat_net_total
    double precision :: heat_const, F(2), dtout
    double precision get_fermi_integral
    integer :: arraymax, idir

!    double precision :: lum_dr(n1,3)  ! for 2D only   a spherical surface luminsoity
    double precision :: dcosj,dr,dtheta

    double precision :: area_2D

!    character(len=100) filename

    integer gain_radius_nue(1)
    integer gain_radius_nua(1)

    integer :: ix^D, ix_next^D, ix_prev^D, neu_num, theta, phi
!output
    double precision, dimension(ixI^S)           :: depsdt !energy change in cell due to neutrinos, indices <radial position> in erg/g/s
    double precision, dimension(ixI^S)           :: dyedt !change in electron number fraction, indices <radial position> in number fraction / s
    double precision, dimension(ixI^S,3)         :: ave_enr_nu !average neutrino energy, indices <radial position:neutrino species> in MeV
    double precision, dimension(ixI^S,3,2)       :: gamma_eff 
    double precision, dimension(ixI^S,3)         :: lum
    double precision, dimension(3)         :: lum_500
      character*1024 filename
      character*1024 filename2


    lum = 0.0d0
    lum_500 = 0.0d0

    lepton_blocking = 0.0d0
    gamma_eff = 0.0d0

    Q_tot = 0.0d0
    R_tot = 0.0d0

    depsdt = 0.0d0
    dyedt = 0.0d0
    R_eff = 0.0d0
    Q_eff = 0.0d0
    heat_rad = 0.0d0
    ave_enr_nu = 0.0d0
    heat_net = 0.0d0
    heat_net_total = 0.0d0


    heat_const = (1.0d0+3.0d0*alpha**2) * sigma_0 / (me_mev**2 * massn_cgs) *0.25d0
!write(*,*) ixOmax1, lx(ixOmax1,r_)
  !  write(*,*) Q_loc(ixOmax1,:), R_loc(ixOmax1,:)
!stop "start of update"

!for test only
!  shock_radius = lx(ixOmax1,r_)
!   write(*,*) shock_radius, "shockradius"
!   write(*,*)  ixOmax1,  "ixOmax1" 

   {do ix^D = ixO^LIM^D \}
       !determine effective rates
       R_eff(ix^D,:) = R_loc(ix^D,:)/(1.0d0+R_loc(ix^D,:)/R_diff(ix^D,:))
       Q_eff(ix^D,:) = Q_loc(ix^D,:)/(1.0d0+Q_loc(ix^D,:)/Q_diff(ix^D,:))

       gamma_eff(ix^D,:,1) = 1.d0/(1.0d0+R_loc(ix^D,:)/R_diff(ix^D,:))
       gamma_eff(ix^D,:,2) = 1.d0/(1.0d0+Q_loc(ix^D,:)/Q_diff(ix^D,:))
 
       do neu_num=1,3 
          if (Q_diff(ix^D,neu_num).eq.0.0d0) then
             Q_eff(ix^D,neu_num) = 0.0d0
             gamma_eff(ix^D,neu_num,2) = 0.0d0
             ave_enr_nu(ix^D,neu_num) = 0.0d0
          else
             ave_enr_nu(ix^D,neu_num) = Q_eff(ix^D,neu_num)/R_eff(ix^D,neu_num)
          endif

          if (R_diff(ix^D,neu_num).eq.0.0d0) then
             R_eff(ix^D,neu_num) = 0.0d0
             gamma_eff(ix^D,neu_num,1) = 0.0d0
             ave_enr_nu(ix^D,neu_num) = 0.0d0
          else
             ave_enr_nu(ix^D,neu_num) = Q_eff(ix^D,neu_num)/R_eff(ix^D,neu_num)
          endif
       enddo

        Q_tot(ix^D) =  (sum(-Q_eff(ix^D,:)))
        R_tot(ix^D) =  -(R_eff(ix^D,1)-R_eff(ix^D,2))
    {enddo^D&\} ! end loop for every points
!
!   write(*,*) Q_eff(:,1) , "Q_eff"
!stop
!    {do ix^D = ixO^LIM^D \}
!         WRITE(*,*) Q_tot(ix^D),"Q_TOT"
!    {enddo^D&\} ! end loop for every points

!        filename = trim(base_filename)//".Q_eff"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*)  global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), Q_eff(ix^D,1:3)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667)
!    endif
!
!
!        filename = trim(base_filename)//".Q_loc"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*)  global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), Q_loc(ix^D,1:3)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667)
!    endif

!        filename = trim(base_filename)//".time_diff"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*)  global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), time_diff(ix^D,1:3,2)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667)
!    endif





   {do ix^D = ixO^LIM^D \}

      ix_prev^D = ix^D - kr(1,^D)
!   shock radius is in cm
!  write(*,*)  shock_radius, "sock" 
!  write(*,*)  "###############"
!  write(*,*) lx(:,r_)
!stop

      if (do_heating .and. lx(ix^D, r_) .lt. shock_radius) then   !shock_radius  is in cm
          !here we sacrifice computational efficiency for clarity...
!        write(*,*) "pass do_heating"

          lepton_blocking(ix^D,1) = 1.0d0/(1.0d0 + exp(eta_nucleons(ix^D,1) - &
               get_fermi_integral(5,eta_nu(ns_location^D(1),1)  )/ &
               get_fermi_integral(4,eta_nu(ns_location^D(1),1))))
          lepton_blocking(ix^D,2) = 1.0d0/(1.0d0 + exp(-eta_nucleons(ix^D,1) - &
               get_fermi_integral(5,eta_nu(ns_location^D(2),2) )/ &
               get_fermi_integral(4,eta_nu(ns_location^D(2),2))))

          F(1:2) = (4.275d0*lleak_tau(ix^D,1:2)+1.15d0)*exp(-2.0d0*lleak_tau(ix^D,1:2))*lepton_blocking(ix^D,1:2)


             heat_rad(ix^D,1) = heat_fac*heat_const * lprim(ix^D, rho_) * mass_fraction(ix^D,1)* lum(ix_prev^D,1) * &
                  lIarea(ix^D) * heat_erms(1)**2 * F(1) * lvol(ix^D) /(lprim(ix^D, alp_)**2)
             heat_rad(ix^D,2) = heat_fac*heat_const * lprim(ix^D, rho_) * mass_fraction(ix^D,2)* lum(ix_prev^D,2) * &
                  lIarea(ix^D) * heat_erms(2)**2 * F(2) * lvol(ix^D) /(lprim(ix^D, alp_)**2)

              heat_rad(ix^D,:) =min(lum(ix_prev^D,:)/lprim(ix^D, alp_)**2,heat_rad(ix^D,:))

          !in MeV/cm^3/s (same units as Q)
          heat_rad(ix^D,1) = heat_rad(ix^D,1) / lvol(ix^D) * erg_to_mev
          heat_rad(ix^D,2) = heat_rad(ix^D,2) / lvol(ix^D) * erg_to_mev
            !red shifted value at infinity
!             lum(ix^D,:) = lum(i-1,j,k,:) +
!             (Q_eff(ix^D,:)-heat_rad(ix^D,:))* &
!                  mev_to_erg*lvol(ix^D)*alp(ix^D)**2*W(ix^D)*(1.0d0+v1(ix^D)*lx1(i))*lpsi(ix^D)**6
            lum(ix^D,:) = lum(ix_prev^D,:) + (Q_eff(ix^D,:)-heat_rad(ix^D,:))*&
                 mev_to_erg*lvol(ix^D)*lprim(ix^D, alp_)**2*lprim(ix^D, psi_)**6

          depsdt(ix^D) = (sum(-Q_eff(ix^D,:)+heat_rad(ix^D,:)))/lprim(ix^D, rho_)*mev_to_erg
          dyedt(ix^D) = -(R_eff(ix^D,1)-R_eff(ix^D,2)-heat_rad(ix^D,1)/heat_em(1)+heat_rad(ix^D,2)/heat_em(2))* &
               massn_cgs/lprim(ix^D, rho_)
          
          heat_net(ix^D) =  (-Q_eff(ix^D,1)-Q_eff(ix^D,2)+heat_rad(ix^D,1)+heat_rad(ix^D,2))*mev_to_erg*lvol(ix^D)

          heat_net_total = heat_net_total + max(heat_net(ix^D),0.0d0)
        
! no need for now
!          !in ergs/g/s (for output)
!          heat_rad(i,1) = heat_rad(i,1) / lrho(i) * mev_to_erg
!          heat_rad(i,2) = heat_rad(i,2) / lrho(i) * mev_to_erg

       else 
          if (lx(ix^D, r_) .lt. shock_radius) then
!          !just cooling
              lum(ix^D,:) = lum(ix_prev^D,:) + Q_eff(ix^D,:)* &
                   mev_to_erg*lvol(ix^D)*lprim(ix^D, alp_)**2*lprim(ix^D, psi_)**6

              depsdt(ix^D) = -sum(Q_eff(ix^D,1:3))*mev_to_erg/lprim(ix^D, rho_)
              dyedt(ix^D) = -(R_eff(ix^D,1)-R_eff(ix^D,2))*massn_cgs/lprim(ix^D, rho_)
          else
!          !outside shocked region
              lum(ix^D,:) = lum(ix_prev^D,:)
          endif
       endif


!    {enddo^D&\} ! end loop for every points
!
!PASS 
!   {do ix^D = ixO^LIM^D \}
!write(*,*)  lum(ix1,1), lx(ix1,r_)/1.d5
! 
!    {enddo^D&\} ! end loop for every points
!
!write(*,*)"###############################"

!   {do ix^D = ixO^LIM^D \}
!
!       if ((do_heating) .and. ((lx(ix^D, r_) .lt. shock_radius) )) then   !shock_radius  is in cm
!
!        need change in 2D###################
!
!          depsdt(ix^D) = (sum(-Q_eff(ix^D,:)+heat_rad(ix^D,:)))/lprim(ix^D, rho_)*mev_to_erg
!          dyedt(ix^D) = -(R_eff(ix^D,1)-R_eff(ix^D,2)-heat_rad(ix^D,1)/heat_em(1)+heat_rad(ix^D,2)/heat_em(2))* &
!               massn_cgs/lprim(ix^D, rho_)
!          
!          heat_net(ix^D) =  (-Q_eff(ix^D,1)-Q_eff(ix^D,2)+heat_rad(ix^D,1)+heat_rad(ix^D,2))*mev_to_erg*lvol(ix^D)
!
!          heat_net_total = heat_net_total + max(heat_net(ix^D),0.0d0)
!       else
!          if (lx(ix^D, r_) .lt. shock_radius) then
!             depsdt(ix^D) = -sum(Q_eff(ix^D,1:3))*mev_to_erg/lprim(ix^D, rho_)
!             dyedt(ix^D) = -(R_eff(ix^D,1)-R_eff(ix^D,2))*massn_cgs/lprim(ix^D, rho_)
!           else
!             depsdt(ix^D) = 0.0d0
!             dyedt(ix^D)  = 0.0d0
!           endif
!       endif


! put into imrpoeved leakage

       if (lprim(ix^D, ye_).le.eos_yemin*1.01d0) then
          !need to surpress any cooling/heating if dyedt.lt.0 near boundary
          if (dyedt(ix^D).lt.0.0d0) then
             dyedt(ix^D) = 0.0d0
             depsdt(ix^D) = 0.0d0
          endif
       endif

      if (lprim(ix^D, ye_).ge.eos_yemax*0.99d0) then
          !need to surpress any cooling/heating if dyedt.lt.0 near boundary
          if (dyedt(ix^D).gt.0.0d0) then
             dyedt(ix^D) = 0.0d0
             depsdt(ix^D) = 0.0d0
          endif
       endif
!
!##############################from gr1d###########################################!
       !set variables for RK, get units right
!          coolingsource(ix^D,2) = W(ix^D)*alp(ix^D)*v(ix^D)*rho(ix^D)*depsdt(ix^D)*eps_gf/time_gf
!          coolingsource(ix^D,3) = W(ix^D)*alp(ix^D)*rho(ix^D)*depsdt(ix^D)*eps_gf/time_gf
!          coolingsource(ix^D,4) = alp(ix^D)*X(ix^D)*rho(ix^D)*dyedt(ix^D)/time_gf

! lcoolingsource is code unit!  gamma_ij is code unit
         do idir = 1, ndir        
!       coolingsource for S_r, S_theta, S_phi
          lcoolingsource(ix^D, idir) = lprim(ix^D, alp_)*lprim(ix^D, psi_)**6 *&
                                   lprim(ix^D, rho_) * rho_gf *depsdt(ix^D)*(eps_gf/time_gf)*&   ! this is the lQ_tot
                                   (lprim(ix^D, beta(idir))*lfac(ix^D)/lprim(ix^D, alp_)  + &
                                    gamma_ij(ix^D, idir, idir)*lfac(ix^D)*& 
                                   ((lprim(ix^D, W_vel(idir))/lfac(ix^D)) - lprim(ix^D, beta(idir))/lprim(ix^D, alp_)) )
         enddo
!                      for tau
          lcoolingsource(ix^D, ndir+1) = lprim(ix^D, psi_)**6*lfac(ix^D)*lprim(ix^D, alp_)*lprim(ix^D, rho_) * rho_gf *depsdt(ix^D)*eps_gf/time_gf
!                      for ye_con
          lcoolingsource(ix^D, ndir+2) = lprim(ix^D, psi_)**6*lprim(ix^D, alp_)*lprim(ix^D, rho_) * rho_gf *dyedt(ix^D)/time_gf 

   {enddo^D&\} ! end loop for every points


!   {do ix^D = ixO^LIM^D \}
!         WRITE(*,*) Q_eff(ix^D,:),"Q_eff"
!         WRITE(*,*) heat_rad(ix^D,:),"heat_rad"
!   {enddo^D&\} ! end loop for every points
     filename = trim(base_filename)//".lumnu"

!code tets
!if (bounce) then
!! in seconds!
!dtout = 5d-5 
!if ( global_time/time_gf .ge. t_dump) then
!     t_dump=t_dump+dtout
! 
!!if (mype==0) then
!!!      write(*,*) lum(ixOmax1,:) ,t_bounce/time_gf ,"lum  and bounce time"
!!
!!              open(unit=666,file=filename,status="unknown",form='formatted',position="append")
!!               write(666,*) global_time/time_gf-t_bounce , dt/time_gf, &
!!                           lum(ixOmax1,1), lleak_tau(5,1), shock_radius/1.0d5,&
!!                                 lprim(5,rho_), lprim(5,temp_), lprim(5,ye_)
!!
!!              close(666)
!!
!!
!!!    if (global_time/time_gf .ge. t_bounce/time_gf+1.0d-4) then
!!!        stop
!!!    endif
!!endif
!endif 
!
!endif
!all pass
!write(*,*) gamma_ij(5,1,1), "gamma11"
!write(*,*)  depsdt(5) ,"depsdt"
!write(*,*)  lfac(5) ,"W"
!
!write(*,*)  lprim(5,beta(1)), lprim(5,alp_) ,"beta1 alp"
!write(*,*) lcoolingsource(5,:) ,"lcoolingsourceall" 
!stop "after coolingsource"
 
!  as using 1d  method to cal, 2d lum    is wrong!   only absorption approxed to
!  be 1d
!#if (DIMENSIONS == 2)
!if (do_heating) then
!   lum(:,:,:,:) = 0.0d0
!    do k=ghosts3+1,leak_global_kmax
!    do j=ghosts2+1,leak_global_jmax
!    do i=ghosts1+1,leak_global_imax
!       if (lx(ix^D, r_) .lt. shock_radius) then
!          lum_dr(:,:) = 0.0d0
!          do theta = ghosts2+1, leak_global_jmax              !   the max theta is pi/2  in gmunu
!               dtheta = lx2i(theta+1) - lx2i(theta)
!!              dcosj = dcos(lx2i(theta+1)) - dcos(lx2i(theta))
!
!              lum_dr(i,:) = lum_dr(i,:) + 2*pi*lx1(i)**2*&
!                                       lalp(i,theta,k)**2 * lpsi(i,theta,k)**6 *&
!!                                      (Q_eff(i,theta,k,:)-heat_rad(i,theta,k,:)) * (-dcosj)  !-ve sign is sin() d(theta) = - dcos
!                                      (Q_eff(i,theta,k,:)-heat_rad(i,theta,k,:)) * dsin(lx2(theta)) *dtheta  !-ve sign is sin() d(theta) = - dcos
!
!          enddo
!
!             lum_dr(i,:) = lum_dr(i,:) *2.0d0            ! include the part of pi/2 to pi
!             dr = lx1i(i+1) - lx1i(i)
!         ! copy for different j
!             lum(i,j,k,:) = lum(i-1,j,k,:) + lum_dr(i,:)*dr * mev_to_erg
!
!!             lum(i,j,k,:) = lum(i,j,k,:) * mev_to_erg
!!             lum(i,j,k,:) = lum(i-1,j,k,:) + (Q_eff(i,j,k,:)-heat_rad(i,j,k,:))* &
!!                  mev_to_erg*lvol(i,j,k)*lprim(ix^D, alp_)**2*lprim(ix^D, psi_)**6
!       else
!             lum(i,j,k,:) = lum(i-1,j,k,:)
!       endif
!    enddo   
!    enddo
!    enddo
!endif
!#endif



!    arraymax = min(max(n1/2,ishock(1)+10),n1)
!!    arraymax = leak_global_imax
!    output = .not.output !alternates output so RK doesn't lead to double output with RK=2
!    
!    if (output) then
!       if (mod(nt,500).eq.0) then
!
!!!          if (.not.small_output) then
!!             filename = trim(adjustl(outdir))//"/blocking.xg"
!!             open(667,file=filename,status='unknown',position='append')
!!             write(667,*) '"Time = ',time-t_bounce
!!             do i=ghosts1+1,arraymax
!!                write(667,"(i8,1P10E15.6)") i,lx1(i),lepton_blocking(i,ghosts2+1,1,1), &
!!                     lepton_blocking(i,ghosts2+1,1,2)
!!             enddo
!!             write(667,*) " "
!!             write(667,*) " "
!!             close(667)
!!!          endif
!
!
!          filename = trim(adjustl(outdir))//"/eta_neutrino.xg"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,*) '"Time = ',time-t_bounce , "ishock = ", ishock(1)
!          do i=ghosts1+1,arraymax
!             write(667,"(i8,1P10E15.6)") i,lx1(i),eta_nue(i,ghosts2+1,1),eta_nua(i,ghosts2+1,1),eta_nux(i,ghosts2+1,1)
!          enddo
!          write(667,*) " "
!          write(667,*) " "
!          close(667)
!
!
!          filename = trim(adjustl(outdir))//"/tau.xg"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,*) '"Time = ',time-t_bounce, "ishock = ", ishock(1)
!          do i=ghosts1+1,arraymax
!             write(667,"(i8,1P10E15.6)") i,lx1(i),leak_tau(i,ghosts2+1,1,1), &
!                  leak_tau(i,ghosts2+1,1,2),leak_tau(i,ghosts2+1,1,3), xxp(i,ghosts2+1,1), xxn(i,ghosts2+1,1),&
!                  xxa(i,ghosts2+1,1), xabar(i,ghosts2+1,1), xzbar(i,ghosts2+1,1)
!          enddo
!          write(667,*) " "
!          write(667,*) " "
!          close(667)
!          
!
!         filename = trim(adjustl(outdir))//"/coolingsource.xg"
!         open(667,file=filename,status='unknown',position='append')
!         write(667,*) "   (2)  - (4)         (5)               (6)", '"Time = ',time-t_bounce, "ishock = ", ishock(1)
!         do i=ghosts1+1,arraymax
!            write(667,"(i8,1P10E15.6)") i,coolingsource(i,ghosts2+1,1,2:4), coolingsource(i,ghosts2+1,1,5),&
!                                             coolingsource(i,ghosts2+1,1,6)
!         enddo
!         write(667,*) " "
!         write(667,*) " "
!         close(667)
!
!
!
!
!!             filename = trim(adjustl(outdir))//"/leak_Q_total.xg"
!!             open(667,file=filename,status='unknown',position='append')
!!             write(667,*) '"Time = ',time-t_bounce, "ishock = ", ishock(1)
!!             do i=ghosts1+1,arraymax
!!                write(667,"(i8,1P10E15.6)") i,lx1(i),Q_tot(i,ghosts2+1,1),R_tot(i,ghosts2+1,1) !Q_diff(i,1)
!!             enddo
!!             write(667,*) " "
!!             write(667,*) " "
!!             close(667)
!!             
!!             filename = trim(adjustl(outdir))//"/leak_Q_nue_pole.xg"  !1D POLE =
!!                                                                      ! equator
!!             open(667,file=filename,status='unknown',position='append')
!!             write(667,*) '"Time = ',time-t_bounce, "ishock = ", ishock(1)
!!             do i=ghosts1+1,arraymax
!!                write(667,"(i8,1P10E15.6)") i,lx1(i),Q_loc(i,ghosts2+1,1,1),Q_eff(i,ghosts2+1,1,1) !Q_diff(i,1)
!!             enddo
!!             write(667,*) " "
!!             write(667,*) " "
!!             close(667)
!!             
!!
!!             filename = trim(adjustl(outdir))//"/gamma_eff.xg"
!!             open(667,file=filename,status='unknown',position='append')
!!             write(667,*) '"Time = ',time-t_bounce,"ishock = ", ishock(1)
!!             do i=ghosts1+1,arraymax
!!                write(667,"(i8,1P10E15.6)") i,lx1(i),gamma_eff(i,ghosts2+1,1,1,2),gamma_eff(i,ghosts2+1,1,2,2),&
!!                                                     gamma_eff(i,ghosts2+1,1,3,2) !Q_diff(i,1)
!!             enddo
!!             write(667,*) " "
!!             write(667,*) " "
!!             close(667)
!!
!!
!          
!       endif
!
!       if (mod(nt,100).eq.0) then
!
!          filename = trim(adjustl(outdir))//"/tau_pole.dat"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,"(1P10E15.6)") time-t_bounce,leak_tau(5,ghosts2+1,1,:)
!          close(667) 
!
!
   {do ix^D = ixO^LIM^D \}
        if (lx(ix^D,r_) .le. 5.0d7) then
                lum_500(1) = max(lum(ix^D,1), lum_500(1))
                lum_500(2) = max(lum(ix^D,2), lum_500(2))
                lum_500(3) = max(lum(ix^D,3), lum_500(3))
        endif
    {enddo^D&\} ! end loop for every points

!    write(*,*) lum_500(:), "lum_500"


        filename = trim(base_filename)//".lum_nu"
!          filename = trim(/users/ho-yin.ng/patrick_gmunu/gmunu/tests/collapse/1D_collapse_leakage/output)//"/lum_nu.dat"
    if (mype == 0) then
          open(667,file=filename,status='unknown',position='append')
          write(667,*) global_time,lum(ixOmax1,1),lum(ixOmax1,2), &
               lum(ixOmax1,3), shock_radius/1.0d5, pns_radius/1.0d5
          close(667) 
    endif


        filename = trim(base_filename)//".leak_tau"
    if (mype == 0) then
          open(667,file=filename,status='unknown',position='append')
          write(667,*) global_time, lleak_tau(5,:), ns_location1(1:3)
          close(667)
    endif

        filename = trim(base_filename)//".eta_nu"
!@    if (mype == 0) then
!@          open(667,file=filename,status='unknown',position='append')
!@          write(667,*) global_time, eta_nu(5,1:3)
!@          close(667)
!@    endif
!        filename = trim(base_filename)//".lum_everygrid"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*)  global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), lum(ix^D,1:3)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667)
!    endif


!
!          filename = trim(adjustl(outdir))//"/nu_spheres.dat"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,"(1P10E15.6)") time-t_bounce,lx(ns_location1(1),1),lx(ns_location1(2),1), &
!               lx(ns_location1(3),1)
!          close(667)
!
!
!       endif    !correspond to  nt100 
!    endif
  end subroutine compute_update

!######################################################################
  subroutine compute_emission(ixI^L, ixO^L)

    use mod_global_parameters

    implicit none
    
    integer, intent(in)             :: ixI^L, ixO^L
    double precision :: beta
    double precision :: pair_const, R_pair, Q_pair
    double precision :: gamma, gamma_const, R_gamma
    double precision :: block_factor_e, block_factor_a, block_factor_x
    double precision :: enr_m, enr_p, enr_tilde_m, enr_tilde_p, const_brems, ksi_brems
    double precision get_fermi_integral

    integer :: ix^D

    !electron & positron capture
    ksi_brems = 0.5d0 
    beta = pi*clight*(1.0d0+3.0d0*alpha**2)*sigma_0/(hc_mevcm**3*me_mev**2)
    const_brems      = 2.08d2*ksi_brems*erg_to_mev        !B30     wtf erg_to_mev, paper doesnt mention

    R_loc = 0.0d0
    Q_loc = 0.0d0
   
         
    {do ix^D = ixO^LIM^D \}

    !electron & positron capture
       R_loc(ix^D,1) = beta*eta_pn(ix^D)*lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D,1))
       Q_loc(ix^D,1) = beta*eta_pn(ix^D)*lprim(ix^D, temp_)**6*get_fermi_integral(5,eta_nucleons(ix^D,1))
       R_loc(ix^D,2) = beta*eta_np(ix^D)*lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D,1))
       Q_loc(ix^D,2) = beta*eta_np(ix^D)*lprim(ix^D, temp_)**6*get_fermi_integral(5,-eta_nucleons(ix^D,1))

       !e-e+ pair processes from Ruffert et al.
       block_factor_e = 1.0d0+exp(eta_nu(ix^D,1)-0.5d0*( &
            get_fermi_integral(4,eta_nucleons(ix^D,1))/get_fermi_integral(3,eta_nucleons(ix^D,1)) + &
            get_fermi_integral(4,-eta_nucleons(ix^D,1))/get_fermi_integral(3,-eta_nucleons(ix^D,1)) &
            ))
       block_factor_a = 1.0d0+exp(eta_nu(ix^D,2)-0.5d0*( &
            get_fermi_integral(4,eta_nucleons(ix^D,1))/get_fermi_integral(3,eta_nucleons(ix^D,1)) + &
            get_fermi_integral(4,-eta_nucleons(ix^D,1))/get_fermi_integral(3,-eta_nucleons(ix^D,1)) &
            ))
       block_factor_x = 1.0d0+exp(eta_nu(ix^D,3)-0.5d0*( &
            get_fermi_integral(4,eta_nucleons(ix^D,1))/get_fermi_integral(3,eta_nucleons(ix^D,1)) + &
            get_fermi_integral(4,-eta_nucleons(ix^D,1))/get_fermi_integral(3,-eta_nucleons(ix^D,1)) &
            ))
       
       enr_m = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D,1))
       enr_p = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D,1))

       enr_tilde_m = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D,1))
       enr_tilde_p = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D,1))

       pair_const = sigma_0*clight/me_mev**2*enr_m*enr_p
       
       R_pair =  pair_const/(36.0d0*block_factor_e*block_factor_a)* &
            ((Cv-Ca)**2+(Cv+Ca)**2)

       R_loc(ix^D,1) = R_loc(ix^D,1) + R_pair
       Q_loc(ix^D,1) = Q_loc(ix^D,1) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
       R_loc(ix^D,2) = R_loc(ix^D,2) + R_pair
       Q_loc(ix^D,2) = Q_loc(ix^D,2) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
       
       R_pair =  pair_const/(9.0d0*block_factor_x**2)*((Cv-Ca)**2+(Cv+Ca-2.0d0)**2)

       R_loc(ix^D,3) = R_loc(ix^D,3) + R_pair
       Q_loc(ix^D,3) = Q_loc(ix^D,3) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

       !plasmon decay from Ruffert et al.
       gamma = gamma_0*sqrt((pi**2+3.0d0*eta_nucleons(ix^D,1)**2)/3.0d0)
       block_factor_e = 1.0d0 + exp(eta_nu(ix^D,1)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
       block_factor_a = 1.0d0 + exp(eta_nu(ix^D,2)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
       block_factor_x = 1.0d0 + exp(eta_nu(ix^D,3)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))

       gamma_const = pi**3*sigma_0*clight*lprim(ix^D, temp_)**8/(me_mev**2*3.0d0*fsc*hc_mevcm**6)* &
            gamma**6*exp(-gamma)*(1.0d0+gamma)
       

       R_gamma = Cv**2*gamma_const/(block_factor_e*block_factor_a)
       R_loc(ix^D,1) = R_loc(ix^D,1) + R_gamma
       Q_loc(ix^D,1) = Q_loc(ix^D,1) + R_gamma*0.5d0*lprim(ix^D, temp_)*(2.0d0+gamma**2/(1.0d0+gamma))
       R_loc(ix^D,2) = R_loc(ix^D,2) + R_gamma
       Q_loc(ix^D,2) = Q_loc(ix^D,2) + R_gamma*0.5d0*lprim(ix^D, temp_)*(2.0d0+gamma**2/(1.0d0+gamma))

       R_gamma = (Cv-1.0d0)**2*4.0d0*gamma_const/block_factor_x**2
       R_loc(ix^D,3) = R_loc(ix^D,3) + R_gamma 
       Q_loc(ix^D,3) = Q_loc(ix^D,3) + R_gamma*0.5d0*lprim(ix^D, temp_)*(2.0d0+gamma**2/(1.0d0+gamma))
       
       !NN Bremsstrahlung (non degenerate limit, BRT06 (with fix to constant out front 1.04 -> 2.0778, Burrows)


!       if (do_NNBrem) then
!          R_pair = 0.3333d0*const_brems * ( mass_fraction(ix^D,1)**2 + mass_fraction(ix^D,2)**2 +&
!                         (28.d0/3.d0)*mass_fraction(ix^D,1)*mass_fraction(ix^D,2))*lprim(ix^D, rho_)**2*lprim(ix^D, temp_)**4.5d0
!          Q_pair = R_pair * lprim(ix^D, temp_)*3.0d0
!  
!          Q_pair = Q_pair *4.0d0
!          R_pair = R_pair *4.0d0
!          
!          Q_loc(ix^D, 3) =  Q_loc(ix^D,3) + Q_pair 
!          R_loc(ix^D, 3) =  R_loc(ix^D,3) + R_pair 
!        endif        


       if (do_NNBrem) then
          R_pair = 0.231d0*(2.0778d2*erg_to_mev)*0.5d0*&
                  (mass_fraction(ix^D,1)**2+mass_fraction(ix^D,2)**2+28.0d0/3.0d0*mass_fraction(ix^D,1)*mass_fraction(ix^D,2))* &
               lprim(ix^D, rho_)**2*lprim(ix^D, temp_)**(4.5d0)
          Q_pair = R_pair*lprim(ix^D, temp_)/0.231d0*0.504d0
          
          R_loc(ix^D,1) = R_loc(ix^D,1) + R_pair
          Q_loc(ix^D,1) = Q_loc(ix^D,1) + Q_pair
          
          R_loc(ix^D,2) = R_loc(ix^D,2) + R_pair
          Q_loc(ix^D,2) = Q_loc(ix^D,2) + Q_pair
          R_loc(ix^D,3) = R_loc(ix^D,3) + 4.0d0*R_pair
          Q_loc(ix^D,3) = Q_loc(ix^D,3) + 4.0d0*Q_pair
       endif

    {enddo^D&\} ! end loop for every points


!   write(*,*)  Q_loc(5,:), R_loc(5,:), "Q_lco  R_lxco"
!stop
  end subroutine compute_emission
!######################################################################
  subroutine compute_diffusion(ixI^L, ixO^L)

    use mod_global_parameters
    implicit none

    double precision :: scattering_kappa,abs_kappa,rate_const,dr
    double precision :: block_factor
    double precision :: get_fermi_integral
    integer          :: ix^D, ix_next^D
    integer, intent(in)             :: ixI^L, ixO^L

    kappa_tilde_nu_scat = 0.0d0
    kappa_tilde_nu_abs = 0.0d0
    eta_pn = 0.0d0
    eta_np = 0.0d0 
    zeta_single = 0.0d0
    chi_single  = 0.0d0
    R_diff = 0.0d0   
    Q_diff = 0.0d0   
!stop "start of diffusion" 
    {do ix^D = ixO^LIM^D \}

       !scattering
       scattering_kappa = lprim(ix^D, rho_)*avo*0.25d0*sigma_0/me_mev**2
       kappa_tilde_nu_scat(ix^D,1,1) = mass_fraction(ix^D,1)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,1,2) = mass_fraction(ix^D,2)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,2,1) = mass_fraction(ix^D,1)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,2,2) = mass_fraction(ix^D,2)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,3,1) = mass_fraction(ix^D,1)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,3,2) = mass_fraction(ix^D,2)*scattering_kappa
       
       scattering_kappa = lprim(ix^D, rho_)*avo*0.0625d0*sigma_0/me_mev**2* &
            mass_fraction(ix^D,5)*(1.0d0-mass_fraction(ix^D,6)/mass_fraction(ix^D,5))**2 ! only have 1 factor of A because kappa multiples the number fraction, not mass fractions
       kappa_tilde_nu_scat(ix^D,1,3) = mass_fraction(ix^D,4)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,2,3) = mass_fraction(ix^D,4)*scattering_kappa
       kappa_tilde_nu_scat(ix^D,3,3) = mass_fraction(ix^D,4)*scattering_kappa

       eta_pn(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,1)-mass_fraction(ix^D,2))/(exp(eta_hat(ix^D))-1.0d0)
       eta_pn(ix^D) = max(eta_pn(ix^D),0.0d0)
       eta_np(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,2)-mass_fraction(ix^D,1))/(exp(-eta_hat(ix^D))-1.0d0)
       eta_np(ix^D) = max(eta_np(ix^D),0.0d0)
      !  write(*,*) eta_np(ix^D), lx(ix^D,r_),"eta_np"
!       if (eta_pn(ix^D) .le. 0.0d0  .or.  eta_np(ix^D) .le. 0.0d0) then 
!           stop "eta_pn, np   -ve"
!       endif


       if (lprim(ix^D, rho_).lt.1.0d11) then
          !non degenerate here, use mass fractions as chemical potentials fail at low densities
          eta_pn(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,2)
          eta_np(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,1)
       endif

       !absorption
!       abs_kappa = lrho(i)*avo*
       abs_kappa = (1.0d0+3.0d0*alpha**2)*0.25d0*sigma_0/me_mev**2
       block_factor = 1.0d0 + exp(eta_nucleons(ix^D,1)-get_fermi_integral(5,eta_nu(ix^D,1))/ &
            get_fermi_integral(4,eta_nu(ix^D,1)))
       kappa_tilde_nu_abs(ix^D,1,1) = eta_np(ix^D)*abs_kappa/block_factor
       kappa_tilde_nu_abs(ix^D,2,1) = 0.0d0 !no absorption of a-type on neutrons
       kappa_tilde_nu_abs(ix^D,3,1) = 0.0d0 !no absorption of x-type neutrinos
       kappa_tilde_nu_abs(ix^D,1,2) = 0.0d0 !no absorption of e-type on protons
       block_factor = 1.0d0 + exp(-eta_nucleons(ix^D,1)-get_fermi_integral(5,eta_nu(ix^D,2))/ &
            get_fermi_integral(4,eta_nu(ix^D,2)))
       kappa_tilde_nu_abs(ix^D,2,2) = eta_pn(ix^D)*abs_kappa/block_factor
       kappa_tilde_nu_abs(ix^D,3,2) = 0.0d0 !no absorption of x-type neutrinos
       kappa_tilde_nu_abs(ix^D,1,3) = 0.0d0 !no absorption on nuclei
       kappa_tilde_nu_abs(ix^D,2,3) = 0.0d0 !no absorption on nuclei
       kappa_tilde_nu_abs(ix^D,3,3) = 0.0d0 !no absorption on nuclei

       !sum up opacities to get zeta_single (again, factoring out energy dependence)
       zeta_single(ix^D,1) = kappa_tilde_nu_scat(ix^D,1,1) + kappa_tilde_nu_scat(ix^D,1,2) + &
            kappa_tilde_nu_scat(ix^D,1,3) + kappa_tilde_nu_abs(ix^D,1,1) + &
            kappa_tilde_nu_abs(ix^D,1,2) + kappa_tilde_nu_abs(ix^D,1,3)

       zeta_single(ix^D,2) = kappa_tilde_nu_scat(ix^D,2,1) + kappa_tilde_nu_scat(ix^D,2,2) + &
            kappa_tilde_nu_scat(ix^D,2,3) + kappa_tilde_nu_abs(ix^D,2,1) + &
            kappa_tilde_nu_abs(ix^D,2,2) + kappa_tilde_nu_abs(ix^D,2,3)

       zeta_single(ix^D,3) = kappa_tilde_nu_scat(ix^D,3,1) + kappa_tilde_nu_scat(ix^D,3,2) + &
            kappa_tilde_nu_scat(ix^D,3,3) + kappa_tilde_nu_abs(ix^D,3,1) + &
            kappa_tilde_nu_abs(ix^D,3,2) + kappa_tilde_nu_abs(ix^D,3,3)
    {enddo^D&\} ! end loop for every points

!write(*,*)  eta_np(5),eta_pn(5), "eta_np"
!stop
!  write(*,*) kappa_tilde_nu_abs(:,1,1), "kappa_tilde_nu_abs_diffusiobn"
!  write(*,*) kappa_tilde_nu_scat(:,1,1), "kappa_tilde_nu_scat_diffusiobn"
!stop

    {do ix^D = ixOmax^D, ixOmin^D, -1 \}
!    do k=ghosts3+1,n3-ghosts3
!    do j=ghosts2+1,n2-ghosts2
!    do i=n1-ghosts1,ghosts1+1,-1
        ix_next^D = ix^D + kr(1,^D) 
        dr = lxi(ix1+1, r_) - lxi(ix1, r_)  
       !integrate zeta to get chi_single, tau with energy dependence factored out
       chi_single(ix^D,1) = chi_single(ix_next^D,1)  + zeta_single(ix^D,1)* dr
       chi_single(ix^D,2) = chi_single(ix_next^D,2)  + zeta_single(ix^D,2)* dr
       chi_single(ix^D,3) = chi_single(ix_next^D,3)  + zeta_single(ix^D,3)* dr
    {enddo^D&\} ! end loop for every points
    
    {do ix^D = ixO^LIM^D \}
       !now we can determine diffusion rates
       rate_const = 4.0d0*pi*clight*zeta_single(ix^D,1)/(hc_mevcm**3*6.0d0*chi_single(ix^D,1)**2)
       R_diff(ix^D,1) = rate_const*lprim(ix^D, temp_)*get_fermi_integral(0,eta_nu(ix^D,1))
       Q_diff(ix^D,1) = rate_const*lprim(ix^D, temp_)**2*get_fermi_integral(1,eta_nu(ix^D,1))
       
       rate_const = 4.0d0*pi*clight*zeta_single(ix^D,2)/(hc_mevcm**3*6.0d0*chi_single(ix^D,2)**2)
       R_diff(ix^D,2) = rate_const*lprim(ix^D, temp_)*get_fermi_integral(0,eta_nu(ix^D,2))
       Q_diff(ix^D,2) = rate_const*lprim(ix^D, temp_)**2*get_fermi_integral(1,eta_nu(ix^D,2))
       
       rate_const = 16.0d0*pi*clight*zeta_single(ix^D,3)/(hc_mevcm**3*6.0d0*chi_single(ix^D,3)**2)
       R_diff(ix^D,3) = rate_const*lprim(ix^D, temp_)*get_fermi_integral(0,eta_nu(ix^D,3))
       Q_diff(ix^D,3) = rate_const*lprim(ix^D, temp_)**2*get_fermi_integral(1,eta_nu(ix^D,3))

    {enddo^D&\} ! end loop for every points

  end subroutine compute_diffusion



  subroutine finding_array_size(n^D)
   use mod_global_parameters
   implicit none
!   integer, intent(in)             :: ixI^L, ixO^L
!      ixG^LL --> ixI^L   ;   ixM^LL --> ixO^L
   integer :: nM_local^D, iigrid, igrid, nM^D
   integer :: nG_local^D, nnG^D
   integer, intent(inout) :: n^D
{^IFONED
   nM_local1 = 0
   nG_local1 = 0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       !without ghostzones
!        nM_local1 = SIZE(ps(igrid)%x(ixMlo1:ixMhi1,1))

        nM_local1 = SIZE(ps(igrid)%x(ixMlo1:ixMhi1,1))+nM_local1
!        nG_local1 = SIZE(ps(igrid)%x(ixGlo1:ixGhi1,1))+nG_local1
  
 !       write(*,*) nM_local1 , "nM_local1"
    enddo
        call MPI_ALLREDUCE(nM_local1, nM1, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)


!        call MPI_ALLREDUCE(nG_local1, nnG1, 1, mpi_integer, &
!                           MPI_SUM, icomm, ierrmpi)
!  write(*,*) nM1, nnG1,  "nM1, nG1"
!        n1 = nM1 + ghostzones1*2
        n1 = nM1 + ghostzones1
!write(*,*) n1,'n1'
!stop  "n1 value"
!      write(*,*) n1 , "n1"
!stop "n1 value "
}

{^IFTWOD
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       !without ghostzones
        nM_local1 = SIZE(ps(igrid)%x(ixMlo1:ixMhi1,1, r_))
        nM_local2 = SIZE(ps(igrid)%x(1,ixMmin2:ixMmax2, theta_))
    enddo
        call MPI_ALLREDUCE(nM_local1, nM1, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)
        call MPI_ALLREDUCE(nM_local2, nM2, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)
!        n1 = nM1 + ghostzones1*2
        n1 = nM1 + ghostzones1
        n2 = nM2 + ghostzones2*2
} 
    
{^IFTHREED
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       !without ghostzones
        nM_local1 = SIZE(ps(igrid)%x(ixMlo1:ixMhi1,1,1,r_))
        nM_local2 = SIZE(ps(igrid)%x(1,ixMmin2:ixMmax2,1,theta_))
        nM_local3 = SIZE(ps(igrid)%x(1,1,ixMmin3:ixMmax3,phi_))
    enddo
        call MPI_ALLREDUCE(nM_local1, nM1, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)
        call MPI_ALLREDUCE(nM_local2, nM2, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)
        call MPI_ALLREDUCE(nM_local3, nM3, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)
!        n1 = nM1 + ghostzones1*2
        n1 = nM1 + ghostzones1
        n2 = nM2 + ghostzones2*2
        n3 = nM3 + ghostzones3*2
} 

  end subroutine finding_array_size

  subroutine initialize_arrays(n^D)
    use mod_global_parameters
    implicit none
    integer, intent(in) :: n^D
    integer :: iigrid, igrid


    allocate(lprim(n^D,nprim)) 
    allocate(lleak_tau(n^D,3)) 

    allocate(lx(n1,ndim)) 
    allocate(lxi(n1,ndim)) 
    lx = 0.0d0
    lxi = 0.0d0
    allocate(lx_map(n1,ndim)) 
    lx_map = 0.0d0
    allocate(lfac(n^D)) 
    allocate(gamma_ij(n^D,3,3)) 
    allocate(lcoolingsource(n^D,ndir+2)) 


    allocate(R_tot(n^D)) 
    allocate(Q_tot(n^D)) 

    allocate(R_eff(n^D,3)) 
    allocate(R_loc(n^D,3)) 
    allocate(R_diff(n^D,3)) 
    allocate(Q_eff(n^D,3)) 
    allocate(Q_loc(n^D,3)) 
    allocate(Q_diff(n^D,3)) 
    allocate(kappa_tilde_nu_scat(n^D,3,3)) 
    allocate(kappa_tilde_nu_abs(n^D,3,3)) 
    allocate(lvol(n^D)) 
    allocate(lIarea(n^D)) 
    allocate(lmass(n^D)) 

    allocate(eta_nucleons(n^D,3)) 
    allocate(eta_hat(n^D)) 
    allocate(eta_nu(n^D,3)) 
    allocate(lepton_blocking(n^D,2)) 

    allocate(mass_fraction(n^D,6)) 
    allocate(eta_pn(n^D)) 
    allocate(eta_np(n^D)) 

    allocate(zeta_single(n^D,3)) 
    allocate(chi_single(n^D,3)) 

    allocate(lprim_local(n^D, nprim))
    allocate(lleak_tau_local(n^D, 3))
       
    chi_single = 0.0d0
    zeta_single = 0.0d0
 
    lprim = 0.0d0
    lleak_tau = 0.0d0   


    lprim_local = 0.0d0
    lleak_tau_local = 0.0d0   

    lfac = 0.0d0 
    gamma_ij = 0.0d0 
    lcoolingsource = 0.0d0 

    R_tot = 0.0d0   
    Q_tot = 0.0d0   
 
    R_eff = 0.0d0
    R_loc = 0.0d0
    R_diff = 0.0d0
    Q_eff = 0.0d0
    Q_loc = 0.0d0
    Q_diff = 0.0d0

    kappa_tilde_nu_scat = 0.0d0
    kappa_tilde_nu_abs = 0.0d0

    lvol = 0.0d0
    lIarea = 0.0d0
    lprim = 0.0d0
    lmass = 0.0d0

    eta_nucleons = 0.0d0
    eta_hat = 0.0d0
    eta_nu = 0.0d0
    lepton_blocking = 0.0d0

    mass_fraction = 0.0d0
    eta_pn = 0.0d0
    eta_np = 0.0d0
  end subroutine initialize_arrays

  subroutine mapping_to_single_core(n^D)
    use mod_global_parameters
  
   implicit none
   integer, intent(in)  ::  n^D
!   double precision, intent(in)    :: prim(ixI^S, 1:nprim)    
!   double precision, intent(in)    :: cons(ixI^S, 1:ncons)
   integer ::  iigrid, igrid, ix^D, i,j, index_prim, flavor, idir
   double precision :: coor(ndim), coor_local(ndim)
!   double precision :: lprim_local(n1, nprim), lleak_tau_local(n1, 3)

!!  find the r_min 
!{^IFONED
!   do iigrid=1,igridstail; igrid=igrids(iigrid);
!        coor_local(1) = minval(ps(igrid)%x(ixMlo1:ixMhi1, 1)
!   enddo  
!       call MPI_ALLREDUCE(coor_local(1), coor(1), 1, mpi_double_precision, &
!                MPI_MIN, icomm, ierrmpi)
!}
!
!{^IFTWOD
!   do iigrid=1,igridstail; igrid=igrids(iigrid);
!        coor_local(1) = minval(ps(igrid)%x(ixMlo1:ixMhi1,1, 1)
!        coor_local(2) = minval(ps(igrid)%x(1, ixMlo1:ixMhi1, 2)
!   enddo  
!       call MPI_ALLREDUCE(coor_local(1:2), coor(1:2), 2, mpi_double_precision, &
!                MPI_MIN, icomm, ierrmpi)
!}
!{^IFTHREED
!   do iigrid=1,igridstail; igrid=igrids(iigrid);
!        coor_local(1) = minval(ps(igrid)%x(ixMlo1:ixMhi1,1, 1)
!        coor_local(2) = minval(ps(igrid)%x(1, ixMlo1:ixMhi1, 2)
!        coor_local(3) = minval(ps(igrid)%x(1,1, ixMlo1:ixMhi1, 3)
!   enddo  
!       call MPI_ALLREDUCE(coor_local(1:3), coor(1:3), 3, mpi_double_precision, &
!                MPI_MIN, icomm, ierrmpi)
!}
!
!
!    {do ix^D = ixO^LIM^D \}

!!  mapping values of coordinates 
!   {do i = ghostzones1 +1 , nG1-ghostzones1  \}
!        do ndir = 1, ndim
!         lx(i,ndir) = coor(ndir)
!        !lxi(i,ndir)....
!        enddo
!
!      do iigrid=1,igridstail; igrid=igrids(iigrid);
!         {do ix^D = ixMmin^D, ixMmax^D \}
!               if (ps(igrid)%x(ix^D,1) > coor(1))) then 
!                   coor_local(1) = ps(igrid)%x(ix1, 1)
!                exit
!               endif 
!         {enddo^D&\}
!         call MPI_ALLREDUCE(coor_local(1:ndim), coor(1:ndim), ndim, mpi_double_precision, &
!                 MPI_MIN, icomm, ierrmpi)
!      enddo 
!   {enddo^D&\}
!
!
!!  mapping all required hydro/metric/leak_tau 
!   do i = ghostzones1+1, nG1-ghostzones1
!      do iigrid=1,igridstail; igrid=igrids(iigrid);
!         do ix1 = ixMlo1, ixMhi1
!
!           if (ps(igrid)%x(ix1,1) .eq. lx1(i)) 
!
!             do index_prim = 1, nprim
!                     lprim(i,index_prim) = prim(ix1, index_prim)     
!             enddo           
!
!                     leak_tau(i,1:3) = cons(ix1, leak_tau(1:3))     
!           endif
!
!         enddo
!
!      enddo
!   enddo
!! check there are any error unmapped positions
!   do i = ghostzones1+1, nG1-ghostzones1
!      do index_prim = 1, nprim
!        if (lprim(i, index_prim) .eq. 0.0d0)
!          write(*,*) "error! this is not mapped", i, index_prim
!          stop
!        endif
!      enddo
!      
!      do flavor = 1, 3
!        if (leak_tau(i, flavor) .eq. 0.0d0)
!          write(*,*) "error! this is not mapped", i, flavor
!          stop
!        endif
!      enddo
!   enddo 
!!##############1D template
   coor_local(1) = 1.d99 
!  find the r_min 
   do iigrid=1,igridstail; igrid=igrids(iigrid);
        coor_local(1) = min(minval(ps(igrid)%x(ixMlo1:ixMhi1, r_)), coor_local(1))
   enddo  
       call MPI_ALLREDUCE(coor_local(1), coor(1), 1, mpi_double_precision, &
                MPI_MIN, icomm, ierrmpi)
!  mapping values of coordinates 
   do i = ghostzones1 +1 , n1 
         lx(i,r_) = coor(1)
         lx_map(i,r_) = coor(1)
        ! x1i
       coor_local(1) = 1.d99
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         do ix1 = ixMhi1, ixMlo1, -1
               if (ps(igrid)%x(ix1,r_) > coor(1)) then
                   coor_local(1) = min(ps(igrid)%x(ix1, r_),coor_local(1))
!EXIT THIS KIND CANT AR , NEED RETHINK THE ALGO 
!               exit
               endif
         enddo
      enddo 
         call MPI_ALLREDUCE(coor_local(1), coor(1), 1, mpi_double_precision, &
                 MPI_MIN, icomm, ierrmpi)

!        write(*,*) coor(1), "coor(1)"
   enddo

!pass 1d
!   if (mype==0) then
!   
!      do i = ghostzones1 +1 , n1-ghostzones1 
!           write(*,*) lx(i,r_), lx(i+1,r_)-lx(i,r_), i  
!      enddo
!   stop "all local x(r_) inside leakage "
!   endif

!checking lx(i, r_)
   do i = ghostzones1 +1 , n1 
     if (lx(i,r_) .le. lx(i-1,r_)) then
        stop "lx1(i,r_) cannot be less than lx1(i-1,r_), there is sth wrong"
     endif
   enddo

!  once u get lx --> find lxi
!  correct lx(r_)  (from old gmunu Grid.F90)
    lxi(ghostzones1+1,r_) = 0.0d0   !the first inner cell must be zero for integral and differential
!   lxi(ghostzones1+2,r_) = lx(ghostzones1+1,r_)*2.0d0
     do i = ghostzones1+2, n1
         lxi(i,r_) = (lx(i-1,r_) + lx(i,r_) )/2.0d0
     enddo


!pass
!!   if (mype==0) then
!    !in parallel write 
!!   do iigrid=1,igridstail; igrid=igrids(iigrid);
!      do i = ghostzones1 +1 , n1-ghostzones1 
!           write(*,*) lx(i,r_), lxi(i,r_),  i  
!      enddo
!!   enddo
!     stop "all local xi(r_) inside leakage "
!!   endif

!  for code test       
!   do iigrid=1,igridstail; igrid=igrids(iigrid);
!        
!      do ix1 = ixMlo1, ixMhi1
!        if (mype==0) then
!              write(*,*) ps(igrid)%cons( ix1, leak_tau(1))
!        endif
!      enddo
!   enddo

!  mapping all required hydro/metric/leak_tau 
   do iigrid=1,igridstail; igrid=igrids(iigrid);

      do ix1 = ixMlo1, ixMhi1
        do i = ghostzones1+1, n1-ghostzones1
!           write(*,*) ps(igrid)%prim(ix1,rho_),ps(igrid)%x(ix1,r_), lx(i,r_)   
           if (ps(igrid)%x(ix1,r_) .eq. lx(i,r_)) then 
!             stop
             do index_prim = 1, nprim
                     lprim_local(i,index_prim) = ps(igrid)%prim(ix1, index_prim)     
             enddo           
!                     lleak_tau_local(i,1:3) = ps(igrid)%cons(ix1, leak_tau(1:3))     

!                write(*,*) "enter right loop["
           endif

         enddo

      enddo
   enddo

!mapping all the lprim_local into lprim from all cores (as some cores have zero values)
! as some primitive quantities are -ve ,  need take MIP_MIN

!     do index_prim = 1, nprim
 !      call MPI_ALLREDUCE(lprim_local(:, rho_), lprim(:, rho_), &
 !                         n1, mpi_double_precision, &
 !                         MPI_MAX, icomm, ierrmpi)

!write(*,*) lprim_local(:,rho_)

       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, rho_), lprim(ghostzones1+1:n1-ghostzones1, rho_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)

!stop
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, ye_), lprim(ghostzones1+1:n1-ghostzones1, ye_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, temp_), lprim(ghostzones1+1:n1-ghostzones1, temp_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, cs2_), lprim(ghostzones1+1:n1-ghostzones1, cs2_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, press_), lprim(ghostzones1+1:n1-ghostzones1, press_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, eps_), lprim(ghostzones1+1:n1-ghostzones1, eps_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, ent_), lprim(ghostzones1+1:n1-ghostzones1, ent_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, alp_), lprim(ghostzones1+1:n1-ghostzones1, alp_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, psi_), lprim(ghostzones1+1:n1-ghostzones1, psi_), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_MAX, icomm, ierrmpi)

!!  take 
      do idir = 1, ndir
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, W_vel(idir)), lprim(ghostzones1+1:n1-ghostzones1, W_vel(idir)), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_SUM, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, beta(idir)), lprim(ghostzones1+1:n1-ghostzones1, beta(idir)), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_SUM, icomm, ierrmpi)
       call MPI_ALLREDUCE(lprim_local(ghostzones1+1:n1-ghostzones1, vecX(idir)), lprim(ghostzones1+1:n1-ghostzones1, vecX(idir)), &
                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
                          MPI_SUM, icomm, ierrmpi)
      enddo
!!     enddo 


!     do i = 1, 3
!       call MPI_ALLREDUCE(lleak_tau_local(ghostzones1+1:n1-ghostzones1,i), lleak_tau(ghostzones1+1:n1-ghostzones1,i), &
!                          n1-ghostzones1-ghostzones1, mpi_double_precision, &
!                          MPI_MAX, icomm, ierrmpi)
!     enddo

!passed for MPI!
  ! leak_tau for next iteration not test yet
!      do i = ghostzones1 +1 , n1-ghostzones1
!           write(*,*) lx(i,r_),  lprim(i,rho_),i
!      enddo
!     stop "all local lprim inside leakage "


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################mapping ghostzones
!mapping inner ghostzones for coordinates
   j = ghostzones1 + 1
  do i = ghostzones1, 1, -1
     lx(i,r_) = -lx(j,r_)
     lx_map(i,r_) = -lx_map(j,r_)
     lxi(i,r_) = -lxi(j+1,r_)
     j = j + 1
  enddo
!mapping outer ghostzones for coordinates
! No need,  because ghostzones1 are (n1-ghostzones1:n1), already mapped
!pass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################mapping ghostzones




! mapping ghostzones values |1 |2 |3 |4 | copy real zones |4 |3 |2 |1 | first 4
! cells
  do i = 1, ghostzones1
      do index_prim = 1, nprim
                     lprim(i,index_prim) = lprim(ghostzones1*2+1-i, index_prim)
      enddo
!                     lleak_tau(i,1:3) = lleak_tau(ghostzones1*2+1-i, 1:3)
  enddo

!pass
!if (mype==0) then
!write(*,*) lprim(:,rho_)
!stop "mapping lprims"
!endif

!        lprim(:,W_vel(1)) = lprim_local(:,W_vel(1))
      
!      do index_prim = 1, nprim
!        lprim(:,index_prim) = lprim_local(:,index_prim)
!      enddo
!        lleak_tau(:,1:3) = lleak_tau_local(:,1:3)

!passed for MPI!
  ! leak_tau for next iteration not test yet
!      do i = ghostzones1 +1 , n1-ghostzones1 
!           write(*,*) lx(i,r_),  lprim(i,alp_), lprim(i, vecX(1)),i  
!      enddo
!     stop "all local lprim inside leakage "



! mapping ghostzones values |93 |94 |95 |96 | copy real zones |96 |95 |94 |93 | last 4
! cells
   j = 1
  do i = n1, n1-ghostzones1+1, -1
      do index_prim = 1, nprim
                     lprim(i,index_prim) = lprim(n1-ghostzones1*2+j, index_prim)     
      enddo
!                     lleak_tau(i,1:3) = lleak_tau(n1-ghostzones1*2+j, 1:3)     
     j = j + 1
  enddo
!if (mype==0) then
!  write(*,*) "ghostzone mapping checker"
!!  write(*,*) lx(1:4,r_), lx(5:8,r_)
!  write(*,*) lprim(1:4,rho_), lprim(5:8,rho_)
!!  write(*,*) lprim(:,rho_)  pass

!  write(*,*) lprim(10,ent_), lprim(10,temp_), lprim(10,ye_) pass!
!stop "mapping ghostzone"
!endif


! check there are any error unmapped positions
  do i = ghostzones1+1, n1
      do index_prim = 1, nprim
       ! if (lprim(i, index_prim) .eq. 0.0d0) then
        if (lprim(i, rho_) .eq. 0.0d0) then
          write(*,*) "error! this is not mapped for prim", i
!          write(*,*) "error! this is not mapped for prim", i, index_prim
          stop
        endif
      enddo


!    if (.not. first_iteration) then      
!      do flavor = 1, 3
!        if (lleak_tau(i, flavor) .eq. 0.0d0) then
!          write(*,*) "error! this is not mapped for leak", i, flavor
!          stop
!        endif
!      enddo 
!     endif
  enddo 
!for code test
!       if (.not. first_iteration) then
!        write(*,*) lleak_tau(5,:), "leak_tau after mapping"
!        stop
!        endif
!stop "end of mapping"
  end subroutine mapping_to_single_core


  subroutine inverse_mapping(n^D)
    use mod_global_parameters
   implicit none
   integer, intent(in) :: n^D
   integer :: igrid, iigrid, i, ix^D, idir
    double precision ::　r_min, r_min_local
!   double precision, dimension(ixI^S) 

!  change back to code unit to match x(:,r_)
 !  lx(:,r_) = lx(:,r_)*length_gf
 
!       do i = ghostzones1+1, n1-ghostzones1
!   write(*,*) lcoolingsource(i,1), lx(i,r_)
!     enddo
!stop



!  note that lx has changed to cgs and then back to code unit => not accuate, so we use lx_map
!  mapping back leak_tau and coolingsource to cons 
  do i = ghostzones1+1, n1-ghostzones1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      do ix1 = ixMlo1, ixMhi1
 
           if (ps(igrid)%x(ix1,r_) .eq. lx_map(i,r_)) then 
             do idir = 1, ndir+2
             ps(igrid)%cons(ix1, coolingsource(idir)) = lcoolingsource(i,idir)
             enddo 
!             ps(igrid)%cons(ix1, leak_tau(1:3)) = lleak_tau(i,1:3)
           endif

      enddo
    enddo
  enddo

!write(*,*) "passed inverse mapping"
!!code tets 
!       call MPI_ALLREDUCE(r_min_local, r_min, &
!                          1, mpi_double_precision, &
!                          MPI_MIN, icomm, ierrmpi)

!pass
!  do i = ghostzones1+1, n1-ghostzones1
!        write(*,*) lcoolingsource(i,1), lx_map(i,r_), "af inverse map"
!  enddo
!stop 

!pass for MPI all cores are valued
! proved twice for coolingsource, passed for leak_tau 
! pass for MPI all cores --> 1D
!   if (mype==0) then
!      do i = ghostzones1+1, n1-ghostzones1
!         write(*,*) lcoolingsource(i,1), lx_map(i,r_), "lcoolingsource(1)" 
!      enddo
!   endif
!!
!write(*,*) "local coolingsource and leak"
!write(*,*) "local coolingsource and leak"
!write(*,*) "local coolingsource and leak"
!write(*,*) "local coolingsource and leak"
!   do iigrid=1,igridstail; igrid=igrids(iigrid);
!     do ix1 = ixMlo1, ixMhi1
!!   if (mype==0) then
!           write(*,*) ps(igrid)%cons(ix1,leak_tau(1)), ps(igrid)%x(ix1,r_) ,"cons(leak_tau)"
!!   endif
!    enddo
!   enddo
!!
!stop "inverse mapping"

  end subroutine inverse_mapping


  subroutine deallocate_arrays
    !deallocate leakage variables
    deallocate(lx)
    deallocate(lxi)
    deallocate(lx_map)
        
    deallocate(lfac)
    deallocate(gamma_ij)
    deallocate(lcoolingsource)


    deallocate(R_tot)
    deallocate(R_eff)
    deallocate(R_loc)
    deallocate(R_diff)
    deallocate(Q_tot)
    deallocate(Q_eff)
    deallocate(Q_loc)
    deallocate(Q_diff)

    deallocate(chi_single)
    deallocate(zeta_single)
    deallocate(kappa_tilde_nu_scat)
    deallocate(kappa_tilde_nu_abs)

    deallocate(lvol)
    deallocate(lIarea)
    deallocate(lmass)

    deallocate(lprim)
    deallocate(lleak_tau)
    deallocate(lprim_local)
    deallocate(lleak_tau_local)

    deallocate(eta_nucleons)
    deallocate(eta_hat)
    deallocate(eta_nu)
    deallocate(lepton_blocking)

    deallocate(mass_fraction)

    deallocate(eta_pn)
    deallocate(eta_np)

  end subroutine deallocate_arrays

!######################################################################
  subroutine find_ruf_tau(ixI^L, ixO^L)

    use mod_global_parameters

    implicit none
    integer, intent(in)             :: ixI^L, ixO^L
    integer jm1,ix^D, ix_next^D
    integer icount
    integer,parameter :: icount_max = 2000
    ! **** EOS ****
    double precision, dimension(ixI^S) :: kappa_const_scat_n
    double precision, dimension(ixI^S) :: kappa_const_scat_p
!    double precision, dimension(ixI^S) :: kappa_const_scat_h
    double precision, dimension(ixI^S) :: kappa_const_abs
    double precision :: csn_0,csp_0,t1,t2
    double precision :: xerr
    double precision, parameter :: xerr_out = 1.0d-12
    double precision, dimension(ixI^S,3) :: kappa_tot   
    double precision, dimension(ixI^S,3) :: kappa_tot_p 
    double precision, dimension(ixI^S,3) :: kappa_scat_n ! 1/cm
    double precision, dimension(ixI^S,3) :: kappa_scat_p ! 1/cm
!    double precision, dimension(ixI^S) :: kappa_scat_h ! 1/cm
    double precision, dimension(ixI^S) :: kappa_abs_n  ! 1/cm
    double precision, dimension(ixI^S) :: kappa_abs_p  ! 1/cm
    double precision, dimension(ixI^S) :: local_eta_nue
    double precision, dimension(ixI^S) :: local_eta_nua
    double precision, dimension(ixI^S) :: local_eta_nux
    double precision, dimension(ixI^S) :: eta_nue_eq


    character*1024:: filename

    double precision :: dr
    double precision :: xlye,xyn,xynp,xyp,xypn
    double precision :: get_fermi_integral

    ! initialize some suff:
    kappa_tot   = 1.0d0
    kappa_tot_p = 1.0d0
    kappa_scat_n = 1.0d-5 ! 1/cm
    kappa_scat_p = 1.0d-5 ! 1/cm
    kappa_abs_n  = 1.0d-5 ! 1/cm
    kappa_abs_p  = 1.0d-5 ! 1/cm

    local_eta_nux = 0.0d0
    local_eta_nue = 0.0d0
    local_eta_nua = 0.0d0

    ! Neutrino-nucleon scattering transport cross-section
    
    ! C_s,N in equation (A1) of Ruffert et al.
    !
    csn_0 = (1.0d0 + 5.0D0*alpha**2) / 24.0d0
    csp_0 = (4.0d0*(Cv-1.0d0)**2 + 5.0d0*alpha**2) / 24.0d0
    !
!write(*,*) ixO^LIM1
!stop 

    {do ix^D = ixO^LIM^D \}
       ! constant parts of kappa (A6)
       t1 = sigma_0 * avo * lprim(ix^D, rho_) * (lprim(ix^D, temp_)/me_mev)**2
       kappa_const_scat_n(ix^D) = csn_0 * t1
       kappa_const_scat_p(ix^D) = csp_0 * t1
       ! (A11) constant part
       kappa_const_abs(ix^D) = (1.0d0+3.0d0*alpha**2)/4.0d0 * t1
    {enddo^D&\} ! end loop for every points

    ! Loop to get converged result for tau.
    ! This is discussed in the text between equations (A5) and
    ! (A6). Note that for the initial iteration the kappas are set to 1.0d-5
    ! unless we have tau from previous time, then use it as starting point
    icount = 1
    xerr = 1.0d0

    do while(xerr.gt.xerr_out .and. icount.lt.icount_max)
       ! copy over current into previous kappa
       kappa_tot_p = kappa_tot
        write(*,*)  icount, 'icount'
       ! set up new kappa based on individual
       ! contributions
       ! nu_e; (A17)
       kappa_tot(ixO^S,1) = &
              kappa_scat_p(ixO^S,1) &
            + kappa_scat_n(ixO^S,1) &
            + kappa_abs_n(ixO^S)
       ! antis; (A18)
       kappa_tot(ixO^S,2) = &
              kappa_scat_p(ixO^S,2) &
            + kappa_scat_n(ixO^S,2) &
            + kappa_abs_p(ixO^S)

       ! nu_xl (A19)
       kappa_tot(ixO^S,3) = &
            + kappa_scat_p(ixO^S,3) &
            + kappa_scat_n(ixO^S,3)

      

!   write(*,*)  kappa_tot(5,1), " kapp_tot"
!stop
       ! Integrate optical depths: Equation (A20)
       ! Note that this is done for energy transport
       if(icount.gt.2) then
!       if(icount.gt.2.or. first_iteration) then
!         call cal_leak_tau(ixI^L, ixO^L, kappa_tot, lx, lxi, leak_tau, lprim)

!        write(*,*)  "pass right way"
         lleak_tau = 0.0d0
         {do ix^D = ixOmax^D, ixOmin^D, -1 \}
               dr = lxi(ix1+1,r_) - lxi(ix1,r_)
               ix_next^D = ix^D + kr(1,^D)
               lleak_tau(ix^D, 1:3) = lleak_tau(ix_next^D,1:3) + &
                     kappa_tot(ix^D,1:3) * dr * lprim(ix^D, psi_)**2
         {enddo^D&\} ! end loop for every points
       endif
!        write(*,*) lleak_tau(5,:), "lleak_tau"


       jm1=1 !energy optical depth, switch to 0 for number
       {do ix^D = ixO^LIM^D \}
          local_eta_nux(ix^D) = 0.0d0   ! (A2)
          ! (A5) equilibrium eta, we have rest masses in our chemical potentials
          ! no need to include mass difference in eta
          eta_nue_eq(ix^D) = eta_nu(ix^D,1)
          ! (A3); note that the ^0 etas are set to 0.0d0
          local_eta_nue(ix^D) =  eta_nue_eq(ix^D) * (1.0d0-exp(-lleak_tau(ix^D,1))) 
          ! (A4)
          local_eta_nua(ix^D) = -eta_nue_eq(ix^D) * (1.0d0-exp(-lleak_tau(ix^D,2))) 
         {enddo^D&\} ! end loop for every points

!        write(*,*)  local_eta_nue(5), local_eta_nua(5), "lcoal_eta_Nue, and nua"

       {do ix^D = ixO^LIM^D \}
          !assuming completely dissociated, valid in side shock
          xlye = lprim(ix^D, ye_)
          ! (A8)
          xyn = (1.0d0-xlye) / (1.0d0 + 2.0d0/3.0d0 * max(eta_nucleons(ix^D,3),0.0d0))
          xyp = xlye / (1.0d0 + 2.0d0/3.0d0*max(eta_nucleons(ix^D,2),0.0d0))
          t1 = exp(-eta_hat(ix^D))
          ! (A13)
          xynp = max((2.0d0*xlye-1.0d0)/ (t1-1.0d0),0.0d0)
          ! (A14)
          xypn = max(xynp * t1,0.0d0)
 
          ! electron neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nue(ix^D)) / & 
               get_fermi_integral(2+jm1,local_eta_nue(ix^D))
          ! (A15)
          t2 = 1.0d0 + exp(eta_nucleons(ix^D,1)-get_fermi_integral(5,local_eta_nue(ix^D)) / &
               get_fermi_integral(4,local_eta_nue(ix^D)))
          ! (A6)

          kappa_scat_n(ix^D,1) = kappa_const_scat_n(ix^D) * xyn  * t1
          kappa_scat_p(ix^D,1) = kappa_const_scat_p(ix^D) * xyp  * t1

          ! (A11)
          kappa_abs_n(ix^D) = kappa_const_abs(ix^D) * xynp * t1 / t2 

          ! anti-electron neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nua(ix^D)) / & 
               get_fermi_integral(2+jm1,local_eta_nua(ix^D))
          ! (A16)
          t2 = 1.0d0 + exp(-eta_nucleons(ix^D,1)-get_fermi_integral(5,local_eta_nua(ix^D)) / &
               get_fermi_integral(4,local_eta_nua(ix^D)))
          ! (A6)
          kappa_scat_n(ix^D,2) = kappa_const_scat_n(ix^D) * xyn  * t1
          kappa_scat_p(ix^D,2) = kappa_const_scat_p(ix^D) * xyp  * t1

          ! (A12)
          kappa_abs_p(ix^D) = kappa_const_abs(ix^D) * xypn * t1 / t2 
          ! nux neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nux(ix^D)) / & 
               get_fermi_integral(2+jm1,local_eta_nux(ix^D))
          ! (A6)
          kappa_scat_n(ix^D,3) = kappa_const_scat_n(ix^D) * xyn * t1
          kappa_scat_p(ix^D,3) = kappa_const_scat_p(ix^D) * xyp * t1

       
       {enddo^D&\} ! end loop for every points

!   write(*,*) kappa_scat_n(5,:), "kappa_scat_n" 
!   write(*,*) kappa_scat_p(5,:), "kappa_scatp" 
       
!pass
!       write(*,*)  kappa_scat_n(5,1), kappa_scat_p(5,1),kappa_abs_n(5), kappa_abs_p(5), "kappa_scat kappa_abs"
!       write(*,*)  (1.0d0 + exp(-eta_nucleons(5,1)-get_fermi_integral(5,local_eta_nua(5)) / &
!               get_fermi_integral(4,local_eta_nua(5)))), local_eta_nua(5), "lep_blocking factor  , local_eta_nua"
!
!         t1 = exp(-eta_hat(5))
!          ! (A13)
!          xynp = max((2.0d0*lprim(5,ye_)-1.0d0)/ (t1-1.0d0),0.0d0)
!
!
!       write(*,*)  xynp*avo*lprim(5,rho_),  "eta_np"

       ! compute relative change xerr
       xerr = 0.0d0
       {do ix^D = ixO^LIM^D \}
          xerr = max(xerr,abs(kappa_tot(ix^D,1)/kappa_tot_p(ix^D,1)-1.0d0))
          xerr = max(xerr,abs(kappa_tot(ix^D,2)/kappa_tot_p(ix^D,2)-1.0d0))
          xerr = max(xerr,abs(kappa_tot(ix^D,3)/kappa_tot_p(ix^D,3)-1.0d0))
       {enddo^D&\} ! end loop for every points

       icount = icount + 1
!write(*,*) xerr, "icount+1"
    enddo


    
    if(icount.ge.icount_max) then
       write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
       stop "icount > icount_max in leakage; simple_leakage.F90"
    endif


    ! Recompute tau based on the most recent kappa_tot
         lleak_tau = 0.0d0
         {do ix^D = ixOmax^D, ixOmin^D, -1 \}
               dr = lxi(ix1+1,r_) - lxi(ix1,r_)
               ix_next^D = ix^D + kr(1,^D)
               lleak_tau(ix^D, 1:3) = lleak_tau(ix_next^D,1:3) + &
                     kappa_tot(ix^D,1:3) * dr * lprim(ix^D, psi_)**2
         {enddo^D&\} ! end loop for every points
!         call cal_leak_tau(ixI^L, ixO^L, kappa_tot, lx, lxi, leak_tau, lprim)


 !   first_iteration = .false.    
!###################  i took equatorial  ns_lcoation  for 2D


  ns_location^D(1:3) = ixOmin^D

  {do ix^D = ixO^LIM^D \}
       !nu-spheres for heating          ns  assumes  shperical symmetric for all
       !                                dimensions  (ns_location_i, 1, 1)
       if (lleak_tau(ix^D,1).gt.0.66666d0) then
          ns_location1(1) = ix1
       endif

       if (lleak_tau(ix^D,2).gt.0.66666d0) then
          ns_location1(2) = ix1
       endif

       if (lleak_tau(ix^D,3).gt.0.66666d0) then
          ns_location1(3) = ix1
       endif
  {enddo^D&\} 

    heat_erms(1) = lprim(ns_location^D(1), temp_)*&
                   sqrt(get_fermi_integral(5,local_eta_nue(ns_location^D(1)))/&
                   get_fermi_integral(3,local_eta_nue(ns_location^D(1))))

    heat_erms(2) = lprim(ns_location^D(2), temp_)*&
                   sqrt(get_fermi_integral(5,local_eta_nua(ns_location^D(2)))/&
                   get_fermi_integral(3,local_eta_nua(ns_location^D(2))))

    heat_erms(3) = lprim(ns_location^D(3), temp_)*&
                   sqrt(get_fermi_integral(5,0.0d0)/ &
                   get_fermi_integral(3,0.0d0))

    heat_em(1) = lprim(ns_location^D(1), temp_)*get_fermi_integral(5,local_eta_nue(ns_location^D(1)))/&
         get_fermi_integral(4,local_eta_nue(ns_location^D(1)))
    heat_em(2) = lprim(ns_location^D(2), temp_)*get_fermi_integral(5,local_eta_nua(ns_location^D(2)))/&
         get_fermi_integral(4,local_eta_nua(ns_location^D(2)))
    heat_em(3) = lprim(ns_location^D(3), temp_)*get_fermi_integral(5,0.0d0)/ &
         get_fermi_integral(4,0.0d0) !not used

    !set degeneracy factors to interpolated values
    eta_nu(ixO^S,1) = local_eta_nue(ixO^S)
    eta_nu(ixO^S,2) = local_eta_nua(ixO^S)
    eta_nu(ixO^S,3) = local_eta_nux(ixO^S)
    
!  write(*,*) lprim(5,temp_), lprim(
!  write(*,*)  eta_nu(5,:), lleak_tau(5,:)
!stop "after leaktau"
 
!        write(*,*) lleak_tau(5,:),"leak_tau after ruf_tau 5 - 10"
!        stop
  end subroutine find_ruf_tau
!######################################################################


  subroutine finding_intermediate_variables(ixI^L, ixO^L)
   ! use mod_global_parameters
    implicit none

    integer, intent(in) :: ixI^L, ixO^L
    double precision    :: W2v2(ixI^S)
    double precision    :: lfac2(ixI^S)    
    integer :: idir, ix^D

     ! calculate lfac  (spherical coor)
         {do ix^D = ixO^LIM^D \}
           lfac(ix^D) = 1.0d0
           gamma_ij(ix^D, 1,1) = 1.0d0*lprim(ix^D, psi_)**4
           gamma_ij(ix^D, 2,2) = lx(ix1,r_)**2*lprim(ix^D, psi_)**4
           gamma_ij(ix^D, 3,3) = gamma_ij(ix^D,2,2) {^NOONED * dsin(lx(ix2, theta_))**2}
       
               W2v2 = 0.0d0
            do idir = 1, ndir
               W2v2(ix^D) = W2v2(ix^D) + gamma_ij(ix^D,idir,idir)*lprim(ix^D, W_vel(idir))**2
            enddo
     
            lfac2(ix^D) = W2v2(ix^D) + 1
     
            lfac(ix^D) = dsqrt( lfac2(ix^D) )
         {enddo^D&\} ! end loop for every points
    
!write(*,*) lfac(:)
!stop "lfac" 
       ! gamma_ij  is code unit!!!!!!!!

 end subroutine finding_intermediate_variables



end module mod_grhd_ccsn_leakage_simple_leakage


!######################################################################
function get_fermi_integral(ifermi,eta)
  implicit none
  integer ifermi
  real*8 get_fermi_integral
  real*8 eta
  real*8 fermi_integral_analytical
  
  fermi_integral_analytical = 0.0d0
  
  ! Expressions for Fermi integrals given in Takahashi et al. 1978 
  if (eta.gt.1.D-3) then  
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             (eta**2/2.0D0 + 1.6449d0)/(1.0D0+EXP(-1.6855d0*eta))
     case (2)
        fermi_integral_analytical = &
             (eta**3/3.0D0 + 3.2899d0*eta)/(1.0D0-EXP(-1.8246d0*eta))
     case (3)
        fermi_integral_analytical = & 
             (eta**4/4.0D0 + 4.9348d0*eta**2+11.3644d0) / &
             (1.0D0+EXP(-1.9039d0*eta))        
     case (4)
        fermi_integral_analytical = &
             (eta**5/5.0D0 + 6.5797d0*eta**3+45.4576d0*eta) / &
             (1.0D0-EXP(-1.9484d0*eta))        
     case (5)
        fermi_integral_analytical = &
             (eta**6/6.0D0 + 8.2247d0*eta**4 + 113.6439d0*eta**2 + &
             236.5323d0)/(1.0D0+EXP(-1.9727d0*eta))
     end select
     
  else
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             EXP(eta)/(1.0D0+0.2159d0*EXP(0.8857d0*eta))
     case (2)
        fermi_integral_analytical = & 
             2.0D0*EXP(eta)/(1.0D0+0.1092d0*EXP(0.8908d0*eta))
     case (3)
        fermi_integral_analytical = & 
             6.0D0*EXP(eta)/(1.0D0+0.0559d0*EXP(0.9069d0*eta))
     case (4)
        fermi_integral_analytical = & 
             24.0D0*EXP(eta)/(1.0D0+0.0287d0*EXP(0.9257d0*eta))
     case (5)
        fermi_integral_analytical = &
             120.0D0*EXP(eta) / (1.0D0 + 0.0147d0*EXP(0.9431d0*eta))
     end select
     
  endif
  get_fermi_integral = fermi_integral_analytical
  
  return
end function get_fermi_integral


