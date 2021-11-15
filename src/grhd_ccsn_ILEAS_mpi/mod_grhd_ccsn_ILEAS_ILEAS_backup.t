!   assumed spherical symmetric:  ns_lcoation
!                                 heat_rad (absorption
!  single core to do ILEAS -->   n1 = nM1 + ghostzones1 * 2

      !optical depth included in cons
!> (only in CFC now)
module mod_grhd_ccsn_ILEAS_ILEAS

!  single core to do ILEAS -->   n1 = nM1 + ghostzones1 * 2
!                         v, v1,v2,v3,time,&          !time
!                         nt,t_bounce
  use mod_eos
  use mod_global_parameters
  use mod_physics
!  use mod_grhd_ccsn_ILEAS_phys
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
      double precision, parameter:: g_A = 1.25d0
      double precision, parameter:: Cp = (4.0d0*(Cv -1.0d0)**2 + 5.0d0*g_A**2)/24.0d0
      double precision, parameter:: Cn = (1.0d0 +5.0d0*g_A**2)/24.0d0
      double precision, parameter:: m_b_mev = 931.494           !Mev/c^2    atmoic mass unit
      double precision, parameter:: m_b_g = 1.66053904d-24      !g    atmoic mass unit
      double precision, parameter:: mp_mev = 938.28             ! rest mass    in Mev/c^2
      double precision, parameter:: mn_mev = 939.57             !
      double precision, parameter   :: planck_mevs =  4.13566769d-21      !planck h    in Mevs
      double precision, parameter   :: ksi_brems = 0.5d0
 
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

      double precision,allocatable  :: chi(:^D&,:) !optical depth with E^2 factored out,indices <radial zones:neutrino species> - 1 / MeV^2
      double precision,allocatable  :: zeta(:^D&,:) !Eq. A21, mean free path with E^2 factored out, indices <radial zones:neutrino species> - 1 / MeV^2 / cm

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

      double precision,allocatable  :: eta_nu_eq(:^D&,:) 

      double precision,allocatable  :: lepton_blocking(:^D&,:) !blocking terms for electrons and positrons, indices <radial position,lepton type> - dimensionless
      
      !mass fractions
      ! 1. xxn   2. xxp   3.  xxa   4. xxh   5. xabar   6. xzbar
      double precision,allocatable  :: mass_fraction(:^D&,:) 
    


      !neutrino sphere
      integer          :: ns_location^D(3) !neutrino sphere location, indices <neutrino species> - radial index
      double precision :: ns_energy(3) !neutrino sphere energies, indices <neutrino species> - in MeV
      logical          :: have_ns_energies !flag that if true uses old energies as starting point, otherwise 18.0MeV is used

      double precision :: heat_erms(3) !root mean squared energy at neutrinosphere (F5/F3), indices <neutrino species> in MeV
      double precision :: heat_em(3) !mean energy at neutrinosphere (F4/F3), indices <neutrino species> in MeV
      double precision :: fermi_eta_3(2)     ! = get_fermi_integral(3, eta_nu)) 1: nue, 2: nua
      double precision :: fermi_eta_4(2)     ! = get_fermi_integral(4, eta_nu))
      double precision :: fermi_eta_5(2)     ! = get_fermi_integral(5, eta_nu))


!ILEAS new varialbes:   find_kappa_opa 
      double precision, allocatable:: kappa_opa_total_j(:^D&,:,:) !total average opacity without energy dependence for neutrino absorption or opacity depth
                                                    !indices <radial zones: ang: phi : neutrino species: : j> - 1 / MeV^2 / cm (B19)
                           !!!!!!!ALL kappa_opa in       -  1 / MeV^2 / cm
      double precision, allocatable :: kappa_opa_abs_j(:^D&,:,:) !leakage module(for opacity): (energy independent) averaged absorption opacities for nu_e,anu_e, but not nu_x
!indices < radial zones: angular zones: phi zones : neutrino species: j> - 1 / MeV^2 / cm    (B13,B14)

      
      ! compute_zeta_n_N      !number den   < rad: ang: phi : p=1,n=2,alpha=3,heavy=4>
      double precision, allocatable :: n_N(:^D&, :)
                ! different from zeta
      double precision, allocatable :: zeta_N(:^D&, :)

      double precision, allocatable :: ksi_NN(:^D&, :)
      double precision, allocatable :: n_b(:^D&)

      double precision,allocatable  :: eta_pn(:^D&) !proton number density corrected for neutron final state blocking, indices <radial position> - dimension(:^D&)less
      double precision,allocatable  :: eta_np(:^D&) !neutron number density corrected for proton final state blocking, indices <radial position> - dimension(:^D&)less

      ! compute timescale   
      double precision, allocatable :: bar_E_j(:^D&,:,:)  ! integrated neu spectrum neu energy eq(22)  <rad: ang: phi : neu species: j>

      ! production timescale   
      double precision, allocatable :: time_prod(:^D&,:,:) !prod time scale    < rad: ang: phi neu species: j>
      ! diffusion timescale   
      double precision, allocatable :: time_diff(:^D&,:,:) !prod time scale    < rad: ang: phi neu species: j>
   
      ! compute Diffusion
      double precision, allocatable :: kappa_diff_total(:^D&,:,:) !Leakage module (for diffusion time scale):(energy dependent) total opacity for each neu, for all j < rad : ang: phi  : neu species: energy space >

!      integer,parameter   :: n_eps = 3                   ! number of centered energy bins
!      integer,parameter   :: n_eps = 14                   ! number of centered energy bins
!code test
      integer,parameter   :: n_eps = 19                   ! number of centered energy bins
!      integer,parameter   :: n_eps = 27                   ! number of centered energy bins
      integer,parameter   :: epsbins_init = 1             !first energy bin
!      integer,parameter   :: epsbins_max = 16            !final energy bin
!code test
!      integer,parameter   :: epsbins_max = 14            !final energy b4
      integer,parameter   :: epsbins_max = 19            !final energy b4
!      integer,parameter   :: epsbins_max = 27            !final energy b4
!      integer,parameter   :: epsbins_max = 3            !final energy b4


      double precision :: eps_binsi(n_eps+1)
      double precision :: eps_bins(n_eps)
      double precision, public :: t_dump = 0.0d0
       
      logical :: do_heating = .false.
      double precision :: heat_fac
      logical :: do_NNBrem = .false. 
      logical, public :: first_iteration
      logical :: output = .true.
 !Public methods
      logical, public  :: bounce = .false.
      double precision, public :: shock_radius, pns_radius   ! in cm 
      double precision, public :: t_bounce
  public :: grhd_ILEAS_init
  public :: grhd_ILEAS 
!  public :: grhd_ILEAS_add_source 
  public :: grhd_ILEAS_activate
contains

  subroutine grhd_ILEAS_activate()
     call grhd_ILEAS_init()
  end subroutine grhd_ILEAS_activate


  subroutine grhd_ILEAS_read_params(files)
    use mod_global_parameters
    implicit none
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grhd_ILEAS_list/ do_heating, heat_fac, do_NNBrem

    write(*,*) " ILEAS is activated !"

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grhd_ILEAS_list, end=111)
111    close(unitpar)
    end do


  end subroutine grhd_ILEAS_read_params


  !> Read this module's parameters from a file
  subroutine grhd_ILEAS_init()
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    implicit none

!   write(*,*) "u passed a wrong route fku"
!   stop
   call grhd_ILEAS_read_params(par_files)

   if (coordinate /= spherical) then 
        stop  "only spherical coordinate can use ILEAS"
   endif

   usr_nompi => grhd_ILEAS_nompi
!   usr_source => grhd_ILEAS_add_source 

  end subroutine grhd_ILEAS_init
!######################################################################

  subroutine grhd_ILEAS_nompi(iit,qt)
    use mod_global_parameters
    use mod_geometry
    implicit none

    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                         :: idir
    integer :: iigrid, igrid, ix^D, m

!    dx^D=rnode(rpdx^D_,igrid);
!    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
!    saveigrid=igrid



    do iigrid=1,igridstail; igrid=igrids(iigrid);
    {do ix^D = ixM^LL \}

       do idir = 1, ndir+2
       ps(igrid)%cons(ix^D, coolingsource(idir)) = 0.0d0
       enddo
    {enddo^D&\}

    enddo
 
    if (bounce) then
!        write(*,*) "ILEAS"
!        write(*,*) "wrong path b4 ILEAS"
!    single core to do it
!code test
       call grhd_ILEAS   !inside cons(leak_tau(:)) changed
    else 
         return
    endif
  end subroutine grhd_ILEAS_nompi
!######################

!  subroutine grhd_ILEAS_add_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,primCT,qt,cons,x)
!    use mod_global_parameters
!    use mod_geometry
!    implicit none
!
!    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
!    double precision, intent(in)    :: qtC, qdt, x(ixI^S, 1:ndim)
!    double precision, intent(in)    :: primCT(ixI^S, 1:nprim), qt
!    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
!!    logical, intent(in)             :: qsourcesplit
!!    logical, intent(inout)          :: active !< Needs to be set to true when active
!    integer                         :: idir,jdir,kdir
!
!   write(*,*) "passed add_source_ILEAS"
!stop
!
!   write(*,*) cons(5,coolingsource(2)), x(5,r_), qdt
!!stop
!!write(*,*) "passed add_source"
!!for code test
!    if (bounce) then
!       {^NOONED
!       do idir =2, ndir
!       ! 1D core treatment
!          where ( (oneDcore).and.(x(ixO^S,r_) <= r_core) )
!            cons(ixO^S, coolingsource(idir)) = 0.0d0
!          end where
!       end do
!       }
!!       write(*,*) "wrong path" 
!       cons(ixO^S, tau_) = cons(ixO^S, tau_) + qdt * cons(ixO^S, coolingsource(ndir+1))
!       cons(ixO^S, ye_con_) = cons(ixO^S, ye_con_) + qdt * cons(ixO^S, coolingsource(ndir+2))
!       do idir =1, ndir
!         cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) + qdt * cons(ixO^S, coolingsource(idir))
!       end do
!    endif
!
!  end subroutine grhd_ILEAS_add_source




!  subroutine grhd_ILEAS(ixO^L, ixI^L, prim, cons, x, dx^D, coolingsource)
  subroutine grhd_ILEAS
    use mod_global_parameters
    use mod_geometry
    !use mod_eos
      implicit none

      integer :: idir, ix^D, ix_next^D
      integer :: n^D,  ixI^L, ixO^L

    integer :: iigrid, igrid, m
    
    !  logical,save :: have_old_tau

      character*1024 :: filename
      double precision :: xmu_e, xmu_n, xmu_p, xtemp,xrho,xye,xeps
      double precision :: xeta_nucleons(1:3), xeta_nu_eq(1)
!######################################################################


!   write(*,*)  ps(1)%cons(5, mom(1)), "cons b4  ILEAS "
!   write(*,*)  ps(1)%prim(5, temp_), "temp b4  ILEAS "
    !finding total array size 
   !i can use ixG^LL%prim(ixG^T) 
!     write(*,*) ps(5)%cons(:,tau_), "b4 inverse mamp"

    n^D = 0
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
        ! 1d only
    {do ix^D = ixImin^D, ixImax^D-1 \}
!    {do ix^D = ixO^LIM^D \}
       lvol(ix^D) = 4.0d0*pi/3.0d0* (lxi(ix1+1, r_)**3 - lxi(ix1, r_)**3 )  ! using inner surface x1    !approxed for 1D 2D 3D 
    {enddo^D&\} ! end loop for every points

        ! 1d only
    {do ix^D = ixI^LIM^D \}
!    {do ix^D = ixO^LIM^D \}
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

!  check how ye affects eta_nue
!        xye =  0.28d0 
!        xtemp = 10.0d0
!        xrho = 5.9d-4
!        call eos_temp_get_all_one_grid(rho = xrho,eps= xeps,  temp = xtemp,&
!                                ye = xye, &                                
!                                mu_e = xmu_e, mu_n = xmu_n, mu_p = xmu_p)
!
!
!       xeta_nucleons(1) = xmu_e/xtemp
!       xeta_nucleons(2) = xmu_p/xtemp
!       xeta_nucleons(3) = xmu_n/xtemp
!       xeta_nu_eq(1) = xeta_nucleons(1) - xeta_nucleons(3) + xeta_nucleons(2)
!write(*,*)   xeta_nu_eq(1), "xeta_nu_eq"
!stop




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

       eta_hat(ix^D) = eta_nucleons(ix^D,3)-eta_nucleons(ix^D,2) - 1.293333d0/lprim(ix^D, temp_)

       eta_nu(ix^D,1) = eta_nucleons(ix^D,1) - eta_nucleons(ix^D,3) + eta_nucleons(ix^D,2) !fully includes effects of rest masses
       eta_nu(ix^D,2) = -eta_nu(ix^D,1)
       eta_nu(ix^D,3) = 0.0d0

!      ILEAS uses,  calculate [1 - f(eta_eq)]
       eta_nu_eq(ix^D,1) = eta_nucleons(ix^D,1) - eta_nucleons(ix^D,3) + eta_nucleons(ix^D,2)
       eta_nu_eq(ix^D,2) = -eta_nu_eq(ix^D,1)
       eta_nu_eq(ix^D,3) = 0.0d0


!     ILEAS needs

        n_b(ix^D) = lprim(ix^D, rho_)/m_b_g

        n_N(ix^D,3) = mass_fraction(ix^D,3)*n_b(ix^D)/4.0d0  ! number den of alpha
        n_N(ix^D,4) = mass_fraction(ix^D,4)*n_b(ix^D)/mass_fraction(ix^D,5) ! number den of heavy       !i checked again and againa must correct for these two
!  code test
!@        n_N(ix^D,3) = 0.0d0
!@        n_N(ix^D,4) = 0.0d0
        


    {enddo^D&\} ! end loop for every points

!  write(*,*)  sum(lprim(:,rho_)), 'total density'
 !     write(*,*) lprim(5, temp_), lprim(5, rho_), lprim(5,ye_), eta_nu_eq(5,1)
   !write(*,*) lprim(5,temp_),lprim(50,temp_),lprim(150,temp_), "temp"

!write(*,*) lprim(5, temp_), "temp"
!eta_hat passed

    {do ix^D = ixI^LIM^D \}
!need map ghostzones as well
    lprim(ix^D,rho_)  = lprim(ix^D, rho_)/rho_gf
    lprim(ix^D,eps_)  = lprim(ix^D, eps_)/eps_gf
    lprim(ix^D,cs2_)  = lprim(ix^D, cs2_) * clight**2
    lprim(ix^D,press_)  = lprim(ix^D, press_)/press_gf
    {enddo^D&\} ! end loop for every points

!without this, so god damn wrong       
! 2d need change
    ix1 = 1
        do while(lprim(ix1,rho_).gt.rho_min_leak.and.ix1.lt. ixOmax1)  
                ix1=ix1+1
        enddo
! write(*,*) ix1, 'ix1'
! write(*,*) sum(lprim(:,rho_)), 'sum of rho'
!   dont add ghostzone1

  ! warning
  !  as deviatives needs ghostzone that we used ixI^LIM^D to do "doloop", but calculate "compute_update", we used ixO^LIM^D
    ixOmax1 = ix1  !+ghostzones1

! set epsbins and epsbinsi
! lower eps,  higher contributions in low temp regio
!    eps_binsi(1) = 5.d0
!    eps_binsi(2) = 6.4d0
!    eps_binsi(3) = 8.4d0
!    eps_binsi(4) = 11.2d0
!    eps_binsi(5) = 15.2d0
!    eps_binsi(6) = 20.7d0
!    eps_binsi(7) = 28.4d0
!    eps_binsi(8) = 39.2d0
!    eps_binsi(9) = 54.3d0
!    eps_binsi(10) = 75.5d0
!    eps_binsi(11) = 105.2d0
!    eps_binsi(12) = 146.7d0
!    eps_binsi(13) = 204.8d0
!    eps_binsi(14) = 286.1d0
!    eps_binsi(15) = 400.0d0  

! code test including 5more bins from 0 MeV to 5 MeV
    eps_binsi(1) = 0.0d0
    eps_binsi(2) = 1.0d0
    eps_binsi(3) = 2.0d0
    eps_binsi(4) = 3.0d0
    eps_binsi(5) = 4.0d0
    eps_binsi(6) = 5.d0
    eps_binsi(7) = 6.4d0
    eps_binsi(8) = 8.4d0
    eps_binsi(9) = 11.2d0
    eps_binsi(10) = 15.2d0
    eps_binsi(11) = 20.7d0
    eps_binsi(12) = 28.4d0
    eps_binsi(13) = 39.2d0
    eps_binsi(14) = 54.3d0
    eps_binsi(15) = 75.5d0
    eps_binsi(16) = 105.2d0
    eps_binsi(17) = 146.7d0
    eps_binsi(18) = 204.8d0
    eps_binsi(19) = 286.1d0
    eps_binsi(20) = 400.0d0  


!    eps_binsi(1) = 0.0d0
!    eps_binsi(2) = 1.0d0
!    eps_binsi(3) = 2.0d0
!    eps_binsi(4) = 3.0d0
!    eps_binsi(5) = 4.0d0
!    eps_binsi(6) = 5.d0
!    eps_binsi(7) = 6.4d0
!    eps_binsi(8) = 8.4d0
!    eps_binsi(9) = 11.2d0
!    eps_binsi(10) = 13.2d0
!    eps_binsi(11) = 15.2d0
!    eps_binsi(12) = 17.5d0
!    eps_binsi(13) = 20.7d0
!    eps_binsi(14) = 24.5d0
!    eps_binsi(15) = 28.4d0
!    eps_binsi(16) = 34.8d0
!    eps_binsi(17) = 39.2d0
!    eps_binsi(18) = 47.0d0
!    eps_binsi(19) = 54.3d0
!    eps_binsi(20) = 65.5d0
!    eps_binsi(21) = 75.5d0
!    eps_binsi(22) = 90.0d0
!    eps_binsi(23) = 105.2d0
!    eps_binsi(24) = 126.3d0
!    eps_binsi(25) = 146.7d0
!    eps_binsi(26) = 175.5d0
!    eps_binsi(27) = 204.8d0
!    eps_binsi(28) = 240.0d0
!    eps_binsi(29) = 286.1d0
!    eps_binsi(30) = 343.d0
!    eps_binsi(31) = 400.0d0  



    eps_bins = 0.0d0
    do m = epsbins_init, epsbins_max
      eps_bins(m) = (eps_binsi(m+1) + eps_binsi(m))/2.0d0
    enddo


!passed  they are different
!write(*,*) shape(lprim(1:ixOmax1,rho_)), shape(lprim(:,rho_))

!#############s
    !find tau and interpolate degeneracy factors
    call find_ruf_tau(ixI^L, ixO^L) 

    !find Q_loc and R_loc             !dependent: lleak_tau
    call compute_emission(ixI^L, ixO^L)

    ! find kappa_opa_diffusion for diffusion timescale   !dependent: lleak_tau, zeta_N, n_N, eps_bins
    ! simple leakage will not use this
!    call find_kappa_opa_diffusion(ixI^L, ixO^L)

    call time_prod_and_bar_E(ixI^L, ixO^L)  !R/Q_loc

    !find Q_diff and R_diff   !bar_E_j, eps_bins 
    call compute_diffusion(ixI^L, ixO^L)

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
!   write(*,*)  ps(1)%cons(5, mom(1)), "cons after  ILEAS "
!   write(*,*)  ps(1)%prim(5, temp_), "temp after  ILEAS "
  end subroutine grhd_ILEAS


!######################################################################
! update will use ixOmin to ixOmax
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

!    character(len=100) filename

    integer gain_radius_nue(1)
    integer gain_radius_nua(1)

    integer :: ix^D, ix_next^D, ix_prev^D, neu, theta, phi, neu_num
!output
    double precision, dimension(ixI^S)           :: depsdt !energy change in cell due to neutrinos, indices <radial position> in erg/g/s
    double precision, dimension(ixI^S)           :: dyedt !change in electron number fraction, indices <radial position> in number fraction / s
    double precision, dimension(ixI^S,3)         :: ave_enr_nu !average neutrino energy, indices <radial position:neutrino species> in MeV
    double precision, dimension(ixI^S,3,2)       :: gamma_eff 
    double precision, dimension(ixI^S,3)         :: lum
    double precision, dimension(3)         :: lum_500
      character*1024 filename

    lum_500 = 0.0d0
    lum = 0.0d0

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

  ! ILEAS  
       do neu = 1, 3   
    !             {do ix^D = ixI^LIM^D \}
                {do ix^D = ixO^LIM^D \}
                    gamma_eff(ix^D,neu,1) = 1.0d0/((1.0d0 + time_diff(ix^D,neu,1)/time_prod(ix^D,neu,1)))            
                    gamma_eff(ix^D,neu,2) = 1.0d0/((1.0d0 + time_diff(ix^D,neu,2)/time_prod(ix^D,neu,2)))
                        
                    R_eff(ix^D,neu) = gamma_eff(ix^D,neu,1) * R_loc(ix^D, neu)
                    Q_eff(ix^D,neu) = gamma_eff(ix^D,neu,2) * Q_loc(ix^D, neu)

                    !    if (gamma_eff(ix^D,neu,1) .eq. 0.0d0) then
                    !          R_eff(ix^D,neu) = 0.0d0
                    !          ave_enr_nu(ix^D,neu) = 0.0d0
                    !    else
                    !          ave_enr_nu(ix^D,neu) = Q_eff(ix^D,neu)/R_eff(ix^D,neu)                        
                    !    endif

                    !    if (gamma_eff(ix^D,neu,2) .eq. 0.0d0) then
                    !          Q_eff(ix^D,neu) = 0.0d0
                    !          ave_enr_nu(ix^D,neu) = 0.0d0
                    !    else
                    !          ave_enr_nu(ix^D,neu) = Q_eff(ix^D,neu)/R_eff(ix^D,neu)                        
                    !    endif
                
                {enddo^D&\} ! end loop for every points
        
      enddo



              !{do ix^D = ixI^LIM^D \}
                {do ix^D = ixO^LIM^D \}
                   Q_tot(ix^D) = sum(-Q_eff(ix^D,:))
                   R_tot(ix^D) =  -(R_eff(ix^D,1)-R_eff(ix^D,2))
             !     write(*,*) Q_eff(ix^D,1) ,"Q_eff_max"
                {enddo^D&\} ! end loop for every points

!        filename = trim(base_filename)//".Q_loc"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*) "time = ", global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), Q_loc(ix^D,1:3)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667) 
!    endif
!
!
!        filename = trim(base_filename)//".ksi_NN"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*) "time = ", global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), ksi_NN(ix^D,1:2)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667) 
!    endif
!
!!
!        filename = trim(base_filename)//".Q_eff"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*) "time = ", global_time-t_bounce*time_gf
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
!        filename = trim(base_filename)//".prim_mapped"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*) "time = ", global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)   alp     psi      rho       ye       temp     "
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), lprim(ix^D,alp_), lprim(ix^D,psi_), lprim(ix^D,rho_),&
!                                    lprim(ix^D,ye_), lprim(ix^D,temp_)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667) 
!    endif





!
!        filename = trim(base_filename)//".time_diff"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*) "time = ", global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), time_diff(ix^D,1:3,2)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667) 
!    endif
!stop

!!!!!simple leakage 
!   {do ix^D = ixO^LIM^D \}
!       !determine effective rates
!       R_eff(ix^D,:) = R_loc(ix^D,:)/(1.0d0+R_loc(ix^D,:)/R_diff(ix^D,:))
!       Q_eff(ix^D,:) = Q_loc(ix^D,:)/(1.0d0+Q_loc(ix^D,:)/Q_diff(ix^D,:))
!
!       gamma_eff(ix^D,:,1) = 1.d0/(1.0d0+R_loc(ix^D,:)/R_diff(ix^D,:))
!       gamma_eff(ix^D,:,2) = 1.d0/(1.0d0+Q_loc(ix^D,:)/Q_diff(ix^D,:))
!
!       do neu_num=1,3
!          if (Q_diff(ix^D,neu_num).eq.0.0d0) then
!             Q_eff(ix^D,neu_num) = 0.0d0
!             gamma_eff(ix^D,neu_num,2) = 0.0d0
!             ave_enr_nu(ix^D,neu_num) = 0.0d0
!          else
!             ave_enr_nu(ix^D,neu_num) = Q_eff(ix^D,neu_num)/R_eff(ix^D,neu_num)
!          endif
!
!          if (R_diff(ix^D,neu_num).eq.0.0d0) then
!             R_eff(ix^D,neu_num) = 0.0d0
!             gamma_eff(ix^D,neu_num,1) = 0.0d0
!             ave_enr_nu(ix^D,neu_num) = 0.0d0
!          else
!             ave_enr_nu(ix^D,neu_num) = Q_eff(ix^D,neu_num)/R_eff(ix^D,neu_num)
!          endif
!       enddo
!
!        Q_tot(ix^D) =  (sum(-Q_eff(ix^D,:)))
!        R_tot(ix^D) =  -(R_eff(ix^D,1)-R_eff(ix^D,2))
!    {enddo^D&\} ! end loop for every points



  ! {do ix^D = ixI^LIM^D \}
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

          lepton_blocking(ix^D,1) = 1.0d0/(1.0d0 + dexp(eta_nucleons(ix^D,1) - &
               get_fermi_integral(5,eta_nu(ns_location^D(1),1)  )/ &
               get_fermi_integral(4,eta_nu(ns_location^D(1),1))))
          lepton_blocking(ix^D,2) = 1.0d0/(1.0d0 + dexp(-eta_nucleons(ix^D,1) - &
               get_fermi_integral(5,eta_nu(ns_location^D(2),2) )/ &
               get_fermi_integral(4,eta_nu(ns_location^D(2),2))))

          F(1:2) = (4.275d0*lleak_tau(ix^D,1:2)+1.15d0)*dexp(-2.0d0*lleak_tau(ix^D,1:2))*lepton_blocking(ix^D,1:2)


             heat_rad(ix^D,1) = heat_fac*heat_const * lprim(ix^D, rho_) * mass_fraction(ix^D,1)* lum(ix_prev^D,1) * &
                  lIarea(ix^D) * heat_erms(1)**2 * F(1) * lvol(ix^D) /(lprim(ix^D, alp_)**2)
             heat_rad(ix^D,2) = heat_fac*heat_const * lprim(ix^D, rho_) * mass_fraction(ix^D,2)* lum(ix_prev^D,2) * &
                  lIarea(ix^D) * heat_erms(2)**2 * F(2) * lvol(ix^D) /(lprim(ix^D, alp_)**2)

              heat_rad(ix^D,:) =min(lum(ix_prev^D,:)/lprim(ix^D, alp_)**2,heat_rad(ix^D,:))
          !in MeV/cm^3/s (same units as Q)
          heat_rad(ix^D,1) = heat_rad(ix^D,1) / lvol(ix^D) * erg_to_mev
          heat_rad(ix^D,2) = heat_rad(ix^D,2) / lvol(ix^D) * erg_to_mev
                   ! only nue, nua
                   do neu = 1, 2
                        if (gamma_eff(ix^D,neu,1) .eq. 0.0d0) then
                              heat_rad(ix^D,neu) = 0.0d0
                        endif

                        if (gamma_eff(ix^D,neu,2) .eq. 0.0d0) then
                              heat_rad(ix^D,neu) = 0.0d0
                        endif
                   enddo
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
!  as leakage only estimate the local emission rate, not consider shocked or not.
!               if we want physical result, should apply shocked region
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

!  for emit for region, but not only for shocked region
!          !just cooling
!              lum(ix^D,:) = lum(ix_prev^D,:) + Q_eff(ix^D,:)* &
!                   mev_to_erg*lvol(ix^D)*lprim(ix^D, alp_)**2*lprim(ix^D, psi_)**6
!
!              depsdt(ix^D) = -sum(Q_eff(ix^D,1:3))*mev_to_erg/lprim(ix^D, rho_)
!              dyedt(ix^D) = -(R_eff(ix^D,1)-R_eff(ix^D,2))*massn_cgs/lprim(ix^D, rho_)
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


! put into imrpoeved ILEAS

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

!      {do ix^D = ixO^LIM^D \}
!        write(*,*)  lcoolingsource(ix^D,:),"coolingsource"
!      {enddo^D&\} ! end loop for every points



   

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
        

   {do ix^D = ixO^LIM^D \}
        if (lx(ix^D,r_) .le. 5.0d7) then
                lum_500(1) = max(lum(ix^D,1), lum_500(1))
                lum_500(2) = max(lum(ix^D,2), lum_500(2))
                lum_500(3) = max(lum(ix^D,3), lum_500(3))
        endif
    {enddo^D&\} ! end loop for every points

!    write(*,*) lum_500(:), "lum_500"

        filename = trim(base_filename)//".lum_nu"
!          filename = trim(/users/ho-yin.ng/patrick_gmunu/gmunu/tests/collapse/1D_collapse_ILEAS/output)//"/lum_nu.dat"
    if (mype == 0) then
          open(667,file=filename,status='unknown',position='append')
!          write(667,"(1P10E15.6)") global_time-t_bounce*time_gf,lum(ixOmax1,1),lum(ixOmax1,2), &
!          if (.not. bounce) then
             write(667,*) global_time,lum(ixOmax1,1),lum(ixOmax1,2), &
                  lum(ixOmax1,3), shock_radius/1.0d5, pns_radius/1.0d5
!          else
!             write(667,*) global_time,lum(ixOmax1,1),lum(ixOmax1,2), &
!                  lum(ixOmax1,3), shock_radius/1.0d5, pns_radius/1.0d5, "bounce!"
!          endif
          close(667) 
    endif

!        filename = trim(base_filename)//".lum_500"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!!          write(667,"(1P10E15.6)") global_time-t_bounce*time_gf,lum_500(1:3)
!!         if (.not. bounce) then
!            write(667,*) global_time,lum_500(1:3)
!!         else
!!            write(667,*) global_time,lum_500(1:3), "bounce"
!!         endif
!          close(667) 
!    endif


!        filename = trim(base_filename)//".lum_everygrid"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!             write(667,*) "time = ", global_time-t_bounce*time_gf
!             write(667,*)  " ix1     lx(r_)     nue    nua   nux"
!                {do ix^D = ixO^LIM^D \}
!                      write(667,*) ix^D, lx(ix^D,r_), lum(ix^D,1:3)
!                {enddo^D&\} ! end loop for every points
!             write(667,*) "#############"
!             write(667,*) "#############"
!          close(667)
!    endif



        filename = trim(base_filename)//".leak_tau"
    if (mype == 0) then
          open(667,file=filename,status='unknown',position='append')
!          write(667,*) global_time-t_bounce*time_gf, lleak_tau(5,:), ns_location1(1:3)
!         if (.not. bounce) then
             write(667,*) global_time, lleak_tau(5,:), ns_location1(1:3)
!         else
!             write(667,*) global_time, lleak_tau(5,:), ns_location1(1:3), "bounce"
!         endif
          close(667)
    endif


!        filename = trim(base_filename)//".eta_nu"
!    if (mype == 0) then
!          open(667,file=filename,status='unknown',position='append')
!!          write(667,*) global_time-t_bounce*time_gf, lleak_tau(5,:), ns_location1(1:3)
!!         if (.not. bounce) then
!             write(667,*) global_time, eta_nu(5,1:3)
!!         else
!!             write(667,*) global_time, eta_nu(5,1:3), "bounce"
!!         endif
!          close(667)
!    endif
!
!          filename = trim(adjustl(outdir))//"/nu_spheres.dat"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,"(1P10E15.6)") time-t_bounce,lx(ns_location1(1),1),lx(ns_location1(2),1), &
!               lx(ns_location1(3),1)
!          close(667)
!
!stop
!
!       endif    !correspond to  nt100 
!    endif
  end subroutine compute_update

!######################################################################


!!######################################################################
!  subroutine compute_emission(ixI^L, ixO^L)
!
!    use mod_global_parameters
!
!    implicit none
!
!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision :: beta
!    double precision :: pair_const, R_pair, Q_pair
!    double precision :: gamma, gamma_const, R_gamma
!    double precision :: block_factor_e, block_factor_a, block_factor_x
!    double precision :: enr_m, enr_p, enr_tilde_m, enr_tilde_p
!    double precision get_fermi_integral
!
!    integer :: ix^D
!
!    !electron & positron capture
!
!    beta = pi*clight*(1.0d0+3.0d0*alpha**2)*sigma_0/(hc_mevcm**3*me_mev**2)
!
!    R_loc = 0.0d0
!    Q_loc = 0.0d0
!
!    {do ix^D = ixO^LIM^D \}
!
!    !electron & positron capture
!       R_loc(ix^D,1) = beta*eta_pn(ix^D)*lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D,1))
!       Q_loc(ix^D,1) = beta*eta_pn(ix^D)*lprim(ix^D, temp_)**6*get_fermi_integral(5,eta_nucleons(ix^D,1))
!       R_loc(ix^D,2) = beta*eta_np(ix^D)*lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D,1))
!       Q_loc(ix^D,2) = beta*eta_np(ix^D)*lprim(ix^D, temp_)**6*get_fermi_integral(5,-eta_nucleons(ix^D,1))
!
!       !e-e+ pair processes from Ruffert et al.
!       block_factor_e = 1.0d0+dexp(eta_nu(ix^D,1)-0.5d0*( &
!            get_fermi_integral(4,eta_nucleons(ix^D,1))/get_fermi_integral(3,eta_nucleons(ix^D,1)) + &
!            get_fermi_integral(4,-eta_nucleons(ix^D,1))/get_fermi_integral(3,-eta_nucleons(ix^D,1)) &
!            ))
!       block_factor_a = 1.0d0+dexp(eta_nu(ix^D,2)-0.5d0*( &
!            get_fermi_integral(4,eta_nucleons(ix^D,1))/get_fermi_integral(3,eta_nucleons(ix^D,1)) + &
!            get_fermi_integral(4,-eta_nucleons(ix^D,1))/get_fermi_integral(3,-eta_nucleons(ix^D,1)) &
!            ))
!       block_factor_x = 1.0d0+dexp(eta_nu(ix^D,3)-0.5d0*( &
!            get_fermi_integral(4,eta_nucleons(ix^D,1))/get_fermi_integral(3,eta_nucleons(ix^D,1)) + &
!            get_fermi_integral(4,-eta_nucleons(ix^D,1))/get_fermi_integral(3,-eta_nucleons(ix^D,1)) &
!            ))
!
!       enr_m = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D,1))
!       enr_p = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D,1))
!
!       enr_tilde_m = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D,1))
!       enr_tilde_p = 8.0d0*pi/hc_mevcm**3*lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D,1))
!
!       pair_const = sigma_0*clight/me_mev**2*enr_m*enr_p
!
!       R_pair =  pair_const/(36.0d0*block_factor_e*block_factor_a)* &
!            ((Cv-Ca)**2+(Cv+Ca)**2)
!
!!      M1 compare, no ee for nue nua
!!       R_loc(ix^D,1) = R_loc(ix^D,1) + R_pair
!!       Q_loc(ix^D,1) = Q_loc(ix^D,1) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
!!       R_loc(ix^D,2) = R_loc(ix^D,2) + R_pair
!!       Q_loc(ix^D,2) = Q_loc(ix^D,2) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
!
!       R_pair =  pair_const/(9.0d0*block_factor_x**2)*((Cv-Ca)**2+(Cv+Ca-2.0d0)**2)
!
!       R_loc(ix^D,3) = R_loc(ix^D,3) + R_pair
!       Q_loc(ix^D,3) = Q_loc(ix^D,3) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
!
!!      M1 compare , no plasmon
!!       !plasmon decay from Ruffert et al.
!!       gamma = gamma_0*sqrt((pi**2+3.0d0*eta_nucleons(ix^D,1)**2)/3.0d0)
!!       block_factor_e = 1.0d0 + dexp(eta_nu(ix^D,1)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
!!       block_factor_a = 1.0d0 + dexp(eta_nu(ix^D,2)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
!!       block_factor_x = 1.0d0 + dexp(eta_nu(ix^D,3)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
!!
!!       gamma_const = pi**3*sigma_0*clight*lprim(ix^D, temp_)**8/(me_mev**2*3.0d0*fsc*hc_mevcm**6)* &
!!            gamma**6*dexp(-gamma)*(1.0d0+gamma)
!!
!!       R_gamma = Cv**2*gamma_const/(block_factor_e*block_factor_a)
!!       R_loc(ix^D,1) = R_loc(ix^D,1) + R_gamma
!!       Q_loc(ix^D,1) = Q_loc(ix^D,1) + R_gamma*0.5d0*lprim(ix^D, temp_)*(2.0d0+gamma**2/(1.0d0+gamma))
!!       R_loc(ix^D,2) = R_loc(ix^D,2) + R_gamma
!!       Q_loc(ix^D,2) = Q_loc(ix^D,2) + R_gamma*0.5d0*lprim(ix^D, temp_)*(2.0d0+gamma**2/(1.0d0+gamma))
!!
!!       R_gamma = (Cv-1.0d0)**2*4.0d0*gamma_const/block_factor_x**2
!!       R_loc(ix^D,3) = R_loc(ix^D,3) + R_gamma
!!       Q_loc(ix^D,3) = Q_loc(ix^D,3) + R_gamma*0.5d0*lprim(ix^D, temp_)*(2.0d0+gamma**2/(1.0d0+gamma))
!
!       !NN Bremsstrahlung (non degenerate limit, BRT06 (with fix to constant out front 1.04 -> 2.0778, Burrows)
!
!       ! In M1 comparisons,  no need Brem for nue nua
!!       if (do_NNBrem) then
!!          R_pair = 0.231d0*(2.0778d2*erg_to_mev)*0.5d0*&
!!                  (mass_fraction(ix^D,1)**2+mass_fraction(ix^D,2)**2+28.0d0/3.0d0*mass_fraction(ix^D,1)*mass_fraction(ix^D,2))* &
!!               lprim(ix^D, rho_)**2*lprim(ix^D, temp_)**(4.5d0)
!!          Q_pair = R_pair*lprim(ix^D, temp_)/0.231d0*0.504d0
!
!       !   R_loc(ix^D,1) = R_loc(ix^D,1) + R_pair
!       !   Q_loc(ix^D,1) = Q_loc(ix^D,1) + Q_pair
!
!       !   R_loc(ix^D,2) = R_loc(ix^D,2) + R_pair
!       !   Q_loc(ix^D,2) = Q_loc(ix^D,2) + Q_pair
!!          R_loc(ix^D,3) = R_loc(ix^D,3) + 4.0d0*R_pair
!!          Q_loc(ix^D,3) = Q_loc(ix^D,3) + 4.0d0*Q_pair
!!       endif
!
!!     sum  4 species tgt
!        Q_pair = (2.0778d2*erg_to_mev)*0.5d0*&
!                  (mass_fraction(ix^D,1)**2+mass_fraction(ix^D,2)**2+28.0d0/3.0d0*mass_fraction(ix^D,1)*mass_fraction(ix^D,2))* &
!               lprim(ix^D, rho_)**2*lprim(ix^D, temp_)**(5.5d0)
!        R_pair = Q_pair/3.0d0/lprim(ix^D, temp_) 
!
!        R_loc(ix^D,3) = R_loc(ix^D,3) + R_pair 
!        Q_loc(ix^D,3) = Q_loc(ix^D,3) + Q_pair 
!    {enddo^D&\} ! end loop for every points
!
!
!!   write(*,*)  Q_loc(5,:), R_loc(5,:), "Q_lco  R_lxco"
!!stop
!  end subroutine compute_emission




!#################################################
!!   ILEAS  
subroutine compute_emission(ixI^L, ixO^L)
   implicit none
  double precision :: const_beta, const_ee_nuea, const_ee_nux, const_gamma_nuea, const_gamma_nux, const_brems
  integer :: ix^D
  integer, intent(in)             :: ixI^L, ixO^L

  double precision get_fermi_integral
  double precision eps_beta_nue                 !B22
  double precision, dimension(ixI^S,2)   :: beta_block ! Beta process neu phase blocking  of mean energy of nue and nua  only 
  double precision, dimension(ixI^S)     :: gamma_plasmon
  double precision, dimension(ixI^S,3)   :: ee_block ! e- e+ annihilation phase blocking of mean energy  for 3 neu
  double precision, dimension(ixI^S,3)   :: gamma_block ! plasmon process phase blocking of mean energy for 3 neu
  double precision, dimension(ixI^S)     :: brems_energy !  brems process   mean energy for nux only
  double precision, dimension(ixI^S,2,2) :: Q_beta 
  double precision, dimension(ixI^S,3,2) :: Q_ee  
  double precision, dimension(ixI^S,3,2) :: Q_plasmon
  double precision, dimension(ixI^S,2)   :: Q_brems

  beta_block = 0.0d0
  gamma_block = 0.0d0
  ee_block = 0.0d0
  brems_energy = 0.0d0
  gamma_plasmon = 0.0d0
  Q_loc = 0.0d0
  R_loc = 0.0d0
  Q_ee = 0.0d0
  Q_beta = 0.0d0
  Q_plasmon = 0.0d0
  Q_brems = 0.0d0

  const_beta = (1.d0+3.0d0*g_A**2)*sigma_0*clight*pi/(me_mev**2)/(hc_mevcm**3)         !const of B20, 21

!  const_ee_nuea = (8.0d0 *pi / hc_mevcm**3 )**2*((Cv - Ca)**2 + (Cv + Ca)**2)*sigma_0*clight/(me_mev**2)/72.0d0   ! const of B24
   const_ee_nuea = 0.0d0
  const_ee_nux = (8.0d0*pi/hc_mevcm**3 )**2*((Cv - Ca)**2 + (Cv + Ca - 2.d0)**2)*sigma_0*clight/(me_mev**2)/18.0d0   ! B26
!   const_ee_nux = 0.0d0

! code test, no plasmon decay
!  const_gamma_nuea = (pi**3/3.d0/fsc)*Cv**2*sigma_0*clight/me_mev**2/hc_mevcm**6    !B27
!  const_gamma_nuea = 0.0d0
!  const_gamma_nux  = (4.d0*pi**3/3.d0/fsc)*(Cv - 1.d0)**2*sigma_0*clight/me_mev**2/hc_mevcm**6   !B28
!  const_gamma_nux = 0.0d0

! code test , no brems
!  const_brems      = 2.08d2*ksi_brems    !B30   is fking wrong
  const_brems      = 2.08d2*ksi_brems*erg_to_mev        !B30     wtf erg_to_mev, paper doesnt mention
!  const_brems = 0.0d0

   {do ix^D = ixI^LIM^D \}
  ! {do ix^D = ixO^LIM^D \}

!!    simple leakage in diffusion to redo eta_np
       eta_pn(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,1)-mass_fraction(ix^D,2))/(dexp(eta_hat(ix^D))-1.0d0)
       eta_pn(ix^D) = max(eta_pn(ix^D),0.0d0)
       eta_np(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,2)-mass_fraction(ix^D,1))/(dexp(-eta_hat(ix^D))-1.0d0)
       eta_np(ix^D) = max(eta_np(ix^D),0.0d0)
!       if (eta_pn(ix^D) .le. 0.0d0  .or.  eta_np(ix^D) .le. 0.0d0) then
!           stop "eta_pn, np   -ve"
!       endif


       if (lprim(ix^D, rho_).lt.1.0d11) then
          !non degenerate here, use mass fractions as chemical potentials fail at low densities
          eta_pn(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,2)
          eta_np(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,1)
       endif

!   ILEAS 
        ! nothing


      eps_beta_nue  = max( (lprim(ix^D, temp_)* get_fermi_integral(5,eta_nucleons(ix^D, 1))/get_fermi_integral(4,eta_nucleons(ix^D, 1)) - Qnp) , 0.0d0)

      beta_block(ix^D,1) = 1.d0/( 1.d0 +  dexp(-( eps_beta_nue/lprim(ix^D, temp_) - eta_nu(ix^D, 1))))


      beta_block(ix^D,2) = 1.d0/( 1.d0 +  dexp(-( get_fermi_integral(5,-eta_nucleons(ix^D, 1))/&
                            get_fermi_integral(4,-eta_nucleons(ix^D, 1))  -  eta_nu(ix^D, 2))))

      ee_block(ix^D,1)   = 1.d0/( 1.d0 +  dexp(-( (0.5d0*get_fermi_integral(4,eta_nucleons(ix^D, 1))/&
                            get_fermi_integral(3,eta_nucleons(ix^D, 1)) +&
                            0.5d0*get_fermi_integral(4,-eta_nucleons(ix^D, 1))/get_fermi_integral(3,-eta_nucleons(ix^D, 1)))- eta_nu(ix^D, 1))))

      ee_block(ix^D,2)   = 1.d0/( 1.d0 +  dexp(-( (0.5d0*get_fermi_integral(4,eta_nucleons(ix^D, 1))/&
                            get_fermi_integral(3,eta_nucleons(ix^D, 1)) +&
                            0.5d0*get_fermi_integral(4,-eta_nucleons(ix^D, 1))/get_fermi_integral(3,-eta_nucleons(ix^D, 1)))- eta_nu(ix^D, 2))))

      ee_block(ix^D,3)   = 1.d0/( 1.d0 +  dexp(-( (0.5d0*get_fermi_integral(4,eta_nucleons(ix^D, 1))/&
                            get_fermi_integral(3,eta_nucleons(ix^D, 1)) + 0.5d0*get_fermi_integral(4,-eta_nucleons(ix^D, 1))&
                           /get_fermi_integral(3,-eta_nucleons(ix^D, 1))  )      - eta_nu(ix^D, 3))))

      gamma_plasmon(ix^D) = (5.565d-2)*sqrt((1.d0/3.d0)*(pi**2 + 3.d0*eta_nucleons(ix^D, 1)**2))

      gamma_block(ix^D,1)  = 1.d0/( 1.d0 +  dexp(-( (0.5d0*(2.d0 + gamma_plasmon(ix^D)**2/&
                             (1.d0 + gamma_plasmon(ix^D)))) - eta_nu(ix^D, 1))))

      gamma_block(ix^D,2)  = 1.d0/( 1.d0 +  dexp(-( (0.5d0*(2.d0 + gamma_plasmon(ix^D)**2/&
                             (1.d0 + gamma_plasmon(ix^D)))) - eta_nu(ix^D, 2))))

      gamma_block(ix^D,3)  = 1.d0/( 1.d0 +  dexp(-( (0.5d0*(2.d0 + gamma_plasmon(ix^D)**2/&
                             (1.d0 + gamma_plasmon(ix^D)))) - eta_nu(ix^D, 3))))

      brems_energy(ix^D)  = 3.d0*lprim(ix^D, temp_)
    {enddo^D&\} ! end loop for every points
!pass
!!!!!Beta_process

    {do ix^D = ixI^LIM^D \}
   !{do ix^D = ixO^LIM^D \}

   !!! j = 0  nue
     Q_beta(ix^D,1,1) = const_beta*eta_pn(ix^D)*beta_block(ix^D,1)*&
                        (lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D, 1) - Qnp/lprim(ix^D, temp_)) +&
                         2.d0*Qnp*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1) - Qnp/lprim(ix^D, temp_)) +&
                         Qnp**2*lprim(ix^D, temp_)**3*get_fermi_integral(2,eta_nucleons(ix^D, 1) - Qnp/lprim(ix^D, temp_)))
   !!! j = 1  nue
     Q_beta(ix^D,1,2) = const_beta*eta_pn(ix^D)*beta_block(ix^D,1)*&
                        (lprim(ix^D, temp_)**6*get_fermi_integral(5,eta_nucleons(ix^D, 1) - Qnp/lprim(ix^D, temp_)) +&
                         2.d0*Qnp*lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D, 1) - Qnp/lprim(ix^D, temp_)) +&
                         Qnp**2*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1) - Qnp/lprim(ix^D, temp_)))

   !!! j = 0 nua
     Q_beta(ix^D,2,1) = const_beta*eta_np(ix^D)*beta_block(ix^D,2)*&
                        (lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D, 1)) +&
                         2.d0*Qnp*lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)) +&
                         Qnp**2*lprim(ix^D, temp_)**3*get_fermi_integral(2,-eta_nucleons(ix^D, 1)))

   !!! j = 1 nua
     Q_beta(ix^D,2,2) = const_beta*eta_np(ix^D)*beta_block(ix^D,2)*&
                        (lprim(ix^D, temp_)**6*get_fermi_integral(5,-eta_nucleons(ix^D, 1)) +&
                         3.d0*Qnp*lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D, 1)) +&
                         3.d0*Qnp**2*lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)) +&
                         Qnp**3*lprim(ix^D, temp_)**3*get_fermi_integral(2,-eta_nucleons(ix^D, 1)))

!!!!!! ee process
   !!! j = 0 nue and nua
     Q_ee(ix^D,1,1) = const_ee_nuea*ee_block(ix^D,1)*ee_block(ix^D,2) *&
                      (lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1))*&
                       lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)) +&
                       lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1))*lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)))

     Q_ee(ix^D,2,1) = Q_ee(ix^D,1,1)

   !!! j = 1 nue and nua

     Q_ee(ix^D,1,2) = const_ee_nuea*ee_block(ix^D,1)*ee_block(ix^D,2) *&
                      (lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D, 1))*&
                       lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)) +&
                       lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1))*lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D, 1)))

     Q_ee(ix^D,2,2) = Q_ee(ix^D,1,2)

   !!! j = 0 nux

     Q_ee(ix^D,3,1) = const_ee_nux*ee_block(ix^D,3)**2 *&
                      (lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1))*&
                       lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)) +&
                       lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1))*lprim(ix^D, temp_)**4*get_fermi_integral(3,-eta_nucleons(ix^D, 1)))

   !!! j = 1 nux

     Q_ee(ix^D,3,2) = const_ee_nux*ee_block(ix^D,3)**2 *&
                       (lprim(ix^D, temp_)**5*get_fermi_integral(4,eta_nucleons(ix^D, 1))*lprim(ix^D, temp_)**4*&
                       get_fermi_integral(3,-eta_nucleons(ix^D, 1)) + lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nucleons(ix^D, 1))*&
                       lprim(ix^D, temp_)**5*get_fermi_integral(4,-eta_nucleons(ix^D, 1)))

!!!!!!!Plasmon process   for all 3 neu

   !!! j = 0 nue and nua
     Q_plasmon(ix^D,1,1) = const_gamma_nuea* lprim(ix^D, temp_)**8 *&
                            gamma_plasmon(ix^D)**6*dexp(-gamma_plasmon(ix^D))*&
                           (1.d0+gamma_plasmon(ix^D))*gamma_block(ix^D,1)*gamma_block(ix^D,2)

     Q_plasmon(ix^D,2,1) = Q_plasmon(ix^D,1,1)
   !!! j = 1 nue and nua

     Q_plasmon(ix^D,1,2) = const_gamma_nuea* lprim(ix^D, temp_)**8 *&
                            gamma_plasmon(ix^D)**6*dexp(-gamma_plasmon(ix^D))*&
                           (1.d0+gamma_plasmon(ix^D))*gamma_block(ix^D,1)*&
                           gamma_block(ix^D,2)* (0.5d0*lprim(ix^D, temp_)*&
                           (2.d0 + gamma_plasmon(ix^D)**2/(1.d0 + gamma_plasmon(ix^D))))

     Q_plasmon(ix^D,2,2) = Q_plasmon(ix^D,1,2)

   !!! j = 0 nux

     Q_plasmon(ix^D,3,1) = const_gamma_nux* lprim(ix^D, temp_)**8 *&
                            gamma_plasmon(ix^D)**6*dexp(-gamma_plasmon(ix^D))*&
                           (1.d0+gamma_plasmon(ix^D))*gamma_block(ix^D,3)**2

   !!! j = 1 nux

     Q_plasmon(ix^D,3,2) = const_gamma_nux* lprim(ix^D, temp_)**8 *&
                            gamma_plasmon(ix^D)**6*dexp(-gamma_plasmon(ix^D))*&
                           (1.d0+gamma_plasmon(ix^D))*gamma_block(ix^D,3)**2*&
                           (0.5d0*lprim(ix^D, temp_)*(2.d0 + gamma_plasmon(ix^D)**2/(1.d0 + gamma_plasmon(ix^D))))

!  rubbish paper ILEAS form
!!!!!!!!Brems process for nux only
! !!! j = 1    do   j = 1 first!!    in the paper ,  it times only rho. it actually times rho**2
!     Q_brems(ix^D,2) = const_brems * ( mass_fraction(ix^D,1)**2 + mass_fraction(ix^D,2)**2 +&
!                        (28.d0/3.d0)*mass_fraction(ix^D,1)*mass_fraction(ix^D,2))*lprim(ix^D, rho_)**2*lprim(ix^D, temp_)**5.5d0
! !!! j = 0
!
!     Q_brems(ix^D,1) = Q_brems(ix^D,2)/brems_energy(ix^D)
!  times 4 's form
!     j = 0  first
      Q_brems(ix^D,1) = 0.3333d0*const_brems * ( mass_fraction(ix^D,1)**2 + mass_fraction(ix^D,2)**2 +&
                       (28.d0/3.d0)*mass_fraction(ix^D,1)*mass_fraction(ix^D,2))*lprim(ix^D, rho_)**2*lprim(ix^D, temp_)**4.5d0
      Q_brems(ix^D,2) = Q_brems(ix^D,1) * lprim(ix^D, temp_)*3.0d0 
      
      Q_brems(ix^D,1) = Q_brems(ix^D,1) *4.0d0
      Q_brems(ix^D,2) = Q_brems(ix^D,2) *4.0d0

    {enddo^D&\} ! end loop for every points

    {do ix^D = ixI^LIM^D \}
   !{do ix^D = ixO^LIM^D \}
!!!!!!!!! total production leptons number rate
    R_loc(ix^D,1)  = Q_beta(ix^D,1,1) + Q_ee(ix^D,1,1) + Q_plasmon(ix^D,1,1)

    R_loc(ix^D,2)  = Q_beta(ix^D,2,1) + Q_ee(ix^D,2,1) + Q_plasmon(ix^D,2,1)

    R_loc(ix^D,3)  = Q_ee(ix^D,3,1)  + Q_plasmon(ix^D,3,1) + Q_brems(ix^D,1)

!!!!!!!!! total production energy rate
    Q_loc(ix^D,1)  = Q_beta(ix^D,1,2) + Q_ee(ix^D,1,2) + Q_plasmon(ix^D,1,2)

    Q_loc(ix^D,2)  = Q_beta(ix^D,2,2) + Q_ee(ix^D,2,2) + Q_plasmon(ix^D,2,2)

    Q_loc(ix^D,3)  = Q_ee(ix^D,3,2) + Q_plasmon(ix^D,3,2) + Q_brems(ix^D,2)
    {enddo^D&\} ! end loop for every points


!pass!!!!!  0.0001% difference

!   write(*,*)  Q_loc(5,1,1,:)   ,"Q_loc"
!stop "R_loc"
!   write(*,*)  Q_beta(5,1,1,:,1), Q_beta(5,1,1,:,2)
   !write(*,*)  Q_loc(:,1,1,:)
!   write(*,*)  Q_ee(5,1,1,:,1), Q_ee(5,1,1,:,2)
!   write(*,*)  Q_plasmon(5,1,1,:,1), Q_plasmon(5,1,1,:,2)
!  stop  "at production_rate"
! passed !!!!!!!!  good values!!
end subroutine compute_emission
!###############################################




!! simple leakage diffusion
!  subroutine compute_diffusion(ixI^L, ixO^L)
!
!    use mod_global_parameters
!    implicit none
!
!    double precision :: scattering_kappa,abs_kappa,rate_const,dr
!    double precision :: block_factor
!    double precision :: get_fermi_integral
!    integer          :: ix^D, ix_next^D
!    integer, intent(in)             :: ixI^L, ixO^L
!
!    kappa_tilde_nu_scat = 0.0d0
!    kappa_tilde_nu_abs = 0.0d0
!!    eta_pn = 0.0d0
!!    eta_np = 0.0d0
!    zeta = 0.0d0
!    chi  = 0.0d0
!    R_diff = 0.0d0
!    Q_diff = 0.0d0
!!stop "start of diffusion"
!    {do ix^D = ixO^LIM^D \}
!
!       !scattering
!       scattering_kappa = lprim(ix^D, rho_)*avo*0.25d0*sigma_0/me_mev**2
!       kappa_tilde_nu_scat(ix^D,1,1) = mass_fraction(ix^D,1)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,1,2) = mass_fraction(ix^D,2)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,2,1) = mass_fraction(ix^D,1)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,2,2) = mass_fraction(ix^D,2)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,3,1) = mass_fraction(ix^D,1)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,3,2) = mass_fraction(ix^D,2)*scattering_kappa
!
!       scattering_kappa = lprim(ix^D, rho_)*avo*0.0625d0*sigma_0/me_mev**2* &
!            mass_fraction(ix^D,5)*(1.0d0-mass_fraction(ix^D,6)/mass_fraction(ix^D,5))**2 ! only have 1 factor of A because kappa multiples the number fraction, not mass fractions
!       kappa_tilde_nu_scat(ix^D,1,3) = mass_fraction(ix^D,4)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,2,3) = mass_fraction(ix^D,4)*scattering_kappa
!       kappa_tilde_nu_scat(ix^D,3,3) = mass_fraction(ix^D,4)*scattering_kappa
!
!! simple leakage redo the eta_np
!       eta_pn(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,1)-mass_fraction(ix^D,2))/(dexp(eta_hat(ix^D))-1.0d0)
!       eta_pn(ix^D) = max(eta_pn(ix^D),0.0d0)
!       eta_np(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,2)-mass_fraction(ix^D,1))/(dexp(-eta_hat(ix^D))-1.0d0)
!       eta_np(ix^D) = max(eta_np(ix^D),0.0d0)
!!       if (eta_pn(ix^D) .le. 0.0d0  .or.  eta_np(ix^D) .le. 0.0d0) then
!!           stop "eta_pn, np   -ve"
!!       endif
!
!
!       if (lprim(ix^D, rho_).lt.1.0d11) then
!          !non degenerate here, use mass fractions as chemical potentials fail at low densities
!          eta_pn(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,2)
!          eta_np(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,1)
!       endif
!
!       !absorption
!!       abs_kappa = lrho(i)*avo*
!       abs_kappa = (1.0d0+3.0d0*alpha**2)*0.25d0*sigma_0/me_mev**2
!       block_factor = 1.0d0 + dexp(eta_nucleons(ix^D,1)-get_fermi_integral(5,eta_nu(ix^D,1))/ &
!            get_fermi_integral(4,eta_nu(ix^D,1)))
!       kappa_tilde_nu_abs(ix^D,1,1) = eta_np(ix^D)*abs_kappa/block_factor
!       kappa_tilde_nu_abs(ix^D,2,1) = 0.0d0 !no absorption of a-type on neutrons
!       kappa_tilde_nu_abs(ix^D,3,1) = 0.0d0 !no absorption of x-type neutrinos
!       kappa_tilde_nu_abs(ix^D,1,2) = 0.0d0 !no absorption of e-type on protons
!       block_factor = 1.0d0 + dexp(-eta_nucleons(ix^D,1)-get_fermi_integral(5,eta_nu(ix^D,2))/ &
!            get_fermi_integral(4,eta_nu(ix^D,2)))
!       kappa_tilde_nu_abs(ix^D,2,2) = eta_pn(ix^D)*abs_kappa/block_factor
!       kappa_tilde_nu_abs(ix^D,3,2) = 0.0d0 !no absorption of x-type neutrinos
!       kappa_tilde_nu_abs(ix^D,1,3) = 0.0d0 !no absorption on nuclei
!       kappa_tilde_nu_abs(ix^D,2,3) = 0.0d0 !no absorption on nuclei
!       kappa_tilde_nu_abs(ix^D,3,3) = 0.0d0 !no absorption on nuclei
!
!       !sum up opacities to get zeta (again, factoring out energy dependence)
!       zeta(ix^D,1) = kappa_tilde_nu_scat(ix^D,1,1) + kappa_tilde_nu_scat(ix^D,1,2) + &
!            kappa_tilde_nu_scat(ix^D,1,3) + kappa_tilde_nu_abs(ix^D,1,1) + &
!            kappa_tilde_nu_abs(ix^D,1,2) + kappa_tilde_nu_abs(ix^D,1,3)
!
!       zeta(ix^D,2) = kappa_tilde_nu_scat(ix^D,2,1) + kappa_tilde_nu_scat(ix^D,2,2) + &
!            kappa_tilde_nu_scat(ix^D,2,3) + kappa_tilde_nu_abs(ix^D,2,1) + &
!            kappa_tilde_nu_abs(ix^D,2,2) + kappa_tilde_nu_abs(ix^D,2,3)
!
!       zeta(ix^D,3) = kappa_tilde_nu_scat(ix^D,3,1) + kappa_tilde_nu_scat(ix^D,3,2) + &
!            kappa_tilde_nu_scat(ix^D,3,3) + kappa_tilde_nu_abs(ix^D,3,1) + &
!            kappa_tilde_nu_abs(ix^D,3,2) + kappa_tilde_nu_abs(ix^D,3,3)
!    {enddo^D&\} ! end loop for every points
!
!!write(*,*)  eta_np(5),eta_pn(5), "eta_np"
!!stop
!!  write(*,*) kappa_tilde_nu_abs(:,1,1), "kappa_tilde_nu_abs_diffusiobn"
!!  write(*,*) kappa_tilde_nu_scat(:,1,1), "kappa_tilde_nu_scat_diffusiobn"
!!stop
!
!    {do ix^D = ixOmax^D, ixOmin^D, -1 \}
!!    do k=ghosts3+1,n3-ghosts3
!!    do j=ghosts2+1,n2-ghosts2
!!    do i=n1-ghosts1,ghosts1+1,-1
!        ix_next^D = ix^D + kr(1,^D)
!        dr = lxi(ix1+1, r_) - lxi(ix1, r_)
!       !integrate zeta to get chi, tau with energy dependence factored out
!       chi(ix^D,1) = chi(ix_next^D,1)  + zeta(ix^D,1)* dr
!       chi(ix^D,2) = chi(ix_next^D,2)  + zeta(ix^D,2)* dr
!       chi(ix^D,3) = chi(ix_next^D,3)  + zeta(ix^D,3)* dr
!    {enddo^D&\} ! end loop for every points
!
!    {do ix^D = ixO^LIM^D \}
!       !now we can determine diffusion rates
!       rate_const = 4.0d0*pi*clight*zeta(ix^D,1)/(hc_mevcm**3*6.0d0*chi(ix^D,1)**2)
!       R_diff(ix^D,1) = rate_const*lprim(ix^D, temp_)*get_fermi_integral(0,eta_nu(ix^D,1))
!       Q_diff(ix^D,1) = rate_const*lprim(ix^D, temp_)**2*get_fermi_integral(1,eta_nu(ix^D,1))
!
!       rate_const = 4.0d0*pi*clight*zeta(ix^D,2)/(hc_mevcm**3*6.0d0*chi(ix^D,2)**2)
!       R_diff(ix^D,2) = rate_const*lprim(ix^D, temp_)*get_fermi_integral(0,eta_nu(ix^D,2))
!       Q_diff(ix^D,2) = rate_const*lprim(ix^D, temp_)**2*get_fermi_integral(1,eta_nu(ix^D,2))
!
!       rate_const = 16.0d0*pi*clight*zeta(ix^D,3)/(hc_mevcm**3*6.0d0*chi(ix^D,3)**2)
!       R_diff(ix^D,3) = rate_const*lprim(ix^D, temp_)*get_fermi_integral(0,eta_nu(ix^D,3))
!       Q_diff(ix^D,3) = rate_const*lprim(ix^D, temp_)**2*get_fermi_integral(1,eta_nu(ix^D,3))
!
!    {enddo^D&\} ! end loop for every points
!
!  end subroutine compute_diffusion





!!  ILEAS
!subroutine compute_diffusion(ixI^L, ixO^L)
!    implicit none
!    integer, intent(in)             :: ixI^L, ixO^L
!    integer :: i,j,k,m,neu, ix^D, index, small_region, location_i, location_f
!    integer :: idim
!    double precision :: deps                      !  eps_bins(k+1) - eps_bins(k)
!    double precision, dimension(ixI^S,3,2) :: sum_div
!    double precision, dimension(ixI^S,3,2,n_eps) :: E_j ! energy-dependent neu number and energy < rad: ang: neu species: j: energy space>  (eq.21)
!!    double precision, dimension(16) :: eps_bins
!
!    double precision, dimension(ixI^S,ndim,n_eps) :: dEdx     !originially,   is dEdr (ixI^S, 3, 2, neps)
!    double precision, dimension(ixI^S,n_eps) :: Norm_Grad_E ! Norm of gradient Ej
!    double precision, dimension(ixI^S,n_eps) :: flux_lim ! flux-limiter   <ixI^S>
!    double precision, dimension(ixI^S,ndim,n_eps) :: Deno != flux_lim * (-c/3kappa) * grad_Ej  with ndim  Deno_r, Deno_theta
!    double precision, dimension(ixI^S,3,2,ndim) :: integral ! energy integral of Deno, include integral_r, integral_theta
!    double precision, dimension(ixI^S,ndim) :: div_dEdx 
!    double precision, dimension(ixI^S,3,n_eps) :: neu_blocking
!
!   neu_blocking = 0.0d0
!   dEdx = 0.0d0
!   Norm_Grad_E = 0.0d0
!   flux_lim = 0.0d0
!   Deno = 0.0d0
!   integral = 0.0d0
!   sum_div = 0.0d0
!   integral = 0.0d0
!   div_dEdx = 0.0d0
!   time_diff = 0.0d0
!   E_j = 0.0d0
!
!
!    {do ix^D = ixI^LIM^D \}
!     do m = epsbins_init, epsbins_max
!         E_j(ix^D,1,1,m) = (4.d0*pi/hc_mevcm**3)*eps_bins(m)**2/(1.d0 + dexp(eps_bins(m)/lprim(ix^D, temp_) - eta_nu(ix^D,1)))
!         E_j(ix^D,2,1,m) = (4.d0*pi/hc_mevcm**3)*eps_bins(m)**2/(1.d0 + dexp(eps_bins(m)/lprim(ix^D, temp_) - eta_nu(ix^D,2)))
!         E_j(ix^D,3,1,m) = 4.0d0*(4.d0*pi/hc_mevcm**3)*eps_bins(m)**2/(1.d0 + dexp(eps_bins(m)/lprim(ix^D, temp_) - eta_nu(ix^D,3)))
!
!
!         E_j(ix^D,1,2,m) = (4.d0*pi/hc_mevcm**3)*eps_bins(m)**3/(1.d0 + dexp(eps_bins(m)/lprim(ix^D, temp_) - eta_nu(ix^D,1)))
!         E_j(ix^D,2,2,m) = (4.d0*pi/hc_mevcm**3)*eps_bins(m)**3/(1.d0 + dexp(eps_bins(m)/lprim(ix^D, temp_) - eta_nu(ix^D,2)))
!         E_j(ix^D,3,2,m) = 4.0d0*(4.d0*pi/hc_mevcm**3)*eps_bins(m)**3/(1.d0 + dexp(eps_bins(m)/lprim(ix^D, temp_) - eta_nu(ix^D,3)))
!     enddo
!    {enddo^D&\}
!
!
!    {do ix^D = ixI^LIM^D \}
!     do m = epsbins_init, epsbins_max
!
!      neu_blocking(ix^D,1,m)  = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)/lprim(ix^D, temp_) - eta_nu_eq(ix^D,1)  )))
!      neu_blocking(ix^D,2,m)  = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)/lprim(ix^D, temp_) - eta_nu_eq(ix^D,2)  )))
!      neu_blocking(ix^D,3,m)  = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)/lprim(ix^D, temp_) - eta_nu_eq(ix^D,3)  )))
!
!     enddo
!    {enddo^D&\}
!
!!    do ix1 = 1,  20
!!      write(*,*) ix^D, E_j(ix1,1,1,2), " E_j"
!!    enddo
!!      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!!   write(*,*) E_j(3:7,3,2,2), "3-7 E_j"
!!stop
!
!!############################1D case
!     integral = 0.0d0    
! 
!!E_j(ghostzone) pass!!
!
!      do j = 1, 2     ! j = 1 is numebr,   j = 2 is energy
!       do neu = 1, 3   
!!code test
!!       do neu = 1, 2
!
!        dEdx = 0.0d0
!        flux_lim = 0.0d0
!        Norm_grad_E = 0.0d0
!        Deno = 0.0d0
!!         include  dEdr,  dEdtheta
!         do m = epsbins_init, epsbins_max
!          do idim = 1, ndim
!           call gradient_any_coordinate_five_pt(ixI^L, ixO^L, E_j(:^D&,neu,j,m), dEdx(:^D&,idim,m), idim)  ! 
!!           call gradient_any_coordinate_two_pt(ixI^L, ixO^L, E_j(:^D&,neu,j,m), dEdx(:^D&,idim,m), idim)  ! 
!          enddo
!         enddo ! eps do loop
!
!!  ch3ecked
!!       do ix1 = 1 , 20 
!!   write(*,*)       E_j(ix1,1,1,2) , dEdx(ix1,1,2), lxi(ix1,r_), "E_j    dEdx     r"
!!        enddo
!!stop
!
!! passed..
!!        {do ix^D = ixI^LIM^D \}
!!         write(*,*) E_j(ix^D,neu,j,2), dEdx(ix^D,1,2),  lxi(ix^D,r_) ,"E_j   dE_j/dr   r"
!!        {enddo^D&\}
!!stop
!
!!!       code test
!!        dEdx = 0.0d0
!!         do m = epsbins_init, epsbins_max
!!          do idim = 1, ndim
!!           call gradient_any_coordinate_five_pt(ixI^L, ixO^L, E_j(:^D&,3,2,m), dEdx(:^D&,idim,m), idim)  ! 
!!!           call gradient_any_coordinate_two_pt(ixI^L, ixO^L, E_j(:^D&,neu,j,m), dEdx(:^D&,idim,m), idim)  ! 
!!          enddo
!!         enddo ! eps do loop
!!       write(*,*)  E_j(5,3,2,:), "E_J"
!!    !   write(*,*)  lxi(5:6,r_) , "lxi"
!!       write(*,*)  dEdx(5,1,:), "dEdx"
!!!       write(*,*) kappa_diff_total(5,3,:), "kappa"
!!        norm_grad_E = 0.0d0
!!                do idim = 1, ndim
!!                 Norm_grad_E(5,:) = Norm_grad_E(5,:)  + dEdx(5,idim,:)**2
!!                enddo
!!                 Norm_grad_E(5,:) = sqrt(Norm_grad_E(5,:))
!!                 flux_lim(5,:) = 1.d0/( 1.0d0 + (1.d0/3.d0/kappa_diff_total(5,3,:))*(Norm_grad_E(5,:)/E_j(5,3,2,:)))
!!               write(*,*)  flux_lim(5,:), "flux_lim"
!!
!!stop
!
!        {do ix^D = ixI^LIM^D \}
!!        {do ix^D = ixO^LIM^D \}
!       
!             do m = epsbins_init, epsbins_max
!
!                do idim = 1, ndim
!                 Norm_grad_E(ix^D,m) = Norm_grad_E(ix^D,m)  + dEdx(ix^D,idim,m)**2
!                enddo
!                 Norm_grad_E(ix^D,m) = sqrt(Norm_grad_E(ix^D,m))
!   
!               if (Norm_grad_E(ix^D,m) .lt. 0.0d0) then
!                    
!                   write(*,*)ix^D,Norm_grad_E(ix^D,m), "Norm_grad_E  < 0.0d0"
!                    stop "Norm grad < 0"
!               endif
!
!
!                 flux_lim(ix^D,m) = 1.d0/( 1.0d0 + (1.d0/3.d0/kappa_diff_total(ix^D,neu,m))*(Norm_grad_E(ix^D,m)/E_j(ix^D,neu,j,m)))
!!       code test
!!                 flux_lim(ix^D,m) = 1.d0/( 1.0d0 + (neu_blocking(ix^D,neu,m)/3.d0/kappa_diff_total(ix^D,neu,m))*(Norm_grad_E(ix^D,m)/E_j(ix^D,neu,j,m)))
!
!                !    write(*,*) flux_lim(ix^D), "flux_)lim"
!                 
!                     if (E_j(ix^D,neu,j,m) .eq. 0.0d0 .or. kappa_diff_total(ix^D,neu,m) .eq. 0.0d0) then           !avoid  E_j = 0 --> flux_limit --> NA --> all NAs shit!
!                        flux_lim(ix^D,m) = 0.0d0                                                                    !if set flux_lim = 0.0d0 --> integral for big eps would not contribute
!!                        flux_lim(ix^D) = 1.0d0
!
!!                       write(*,*) ix^D,E_j(ix^D,neu,j,m), kappa_diff_total(ix^D,neu,m), "E_j = 0.0d0 or kappa_diff_total = 0.0d0"
!                     endif
!
!                ! include Deno_x1, Deno_x2
!                do idim = 1, ndim
!                 Deno(ix^D,idim,m) = (-clight/3.d0/kappa_diff_total(ix^D,neu,m))*dEdx(ix^D,idim,m)*(flux_lim(ix^D,m))
!                enddo
!             enddo ! eps do loop
!
!             ! for integral for all driections against  eps
!           do m = epsbins_init, epsbins_max
!!              deps = eps_bins(m) - eps_bins(m-1)    ! should not use the fking  0.0Mev  --> deps = 5 - 0= 5 Mev, 5 times larger
!!              deps = eps_bins(m+1) - eps_bins(m)
!              deps = eps_binsi(m+1) - eps_binsi(m)
!              do idim = 1, ndim
!               integral(ix^D,neu,j,idim) = integral(ix^D,neu,j,idim)  +  Deno(ix^D,idim,m) * deps  !retangle rule
!              enddo
!           enddo ! eps do loop
!          
!         {enddo^D&\}
!
!         enddo
!        enddo
!
!
!
!!        do ix1 =4, 10
!!  write(*,*)  flux_lim(ix1,2), "flux"
!!        enddo
!!write(*,*) '#############'
!!     numerical writing test passed!!!!!
!
!!write(*,*) dEdx(5,1,2), "dEdx"
!        !  you must complete the integration then do other things
!!   end the energy bins integral
!!stop
!
!       
!      do j = 1, 2     ! j = 1 is numebr,   j = 2 is energy
!       do neu = 1, 3   
!         do idim = 1, ndim
!        {do ix^D = ixI^LIM^D \}
!!        {do ix^D = ixO^LIM^D \}
!            integral(ix^D,neu,j,idim) = integral(ix^D,neu,j,idim) * lprim(ix^D, alp_) * lprim(ix^D, psi_)**2
!        {enddo^D&\}
!         enddo
!        enddo
!      enddo
!
!!   write(*,*)  integral(3:7,3,2,1) , "integral"
!!   write(*,*)  lx(3:7,r_), "lx"
!!   write(*,*)  lxi(5:6, r_), "lxi"
!!        do ix1 =4, 10
!!            write(*,*)  integral(ix1,2,2,1),ix1, "integral(1,1)"
!!        enddo
!!write(*,*) '#############'
!!       stop
!
!          sum_div = 0.0d0
!     ! take the shift vector (beta) to be negligible 
!      do j = 1, 2     ! j = 1 is numebr,   j = 2 is energy
!       do neu = 1, 3   
!          div_dEdx = 0.0d0
!
!          do idim = 1, ndim
!!             call div_any_coordinate_five_pt(ixI^L, ixO^L, integral(:^D&,neu,j,idim), div_dEdx(:^D&, idim), idim) 
!             call div_any_coordinate_two_pt(ixI^L, ixO^L, integral(:^D&,neu,j,idim), div_dEdx(:^D&, idim), idim) 
!          enddo
!
!            {do ix^D = ixO^LIM^D \}
!               do idim = 1, ndim
!                  sum_div(ix^D,neu,j) = sum_div(ix^D,neu,j) + div_dEdx(ix^D,idim)
!               enddo
!            {enddo^D&\}
!          
!            {do ix^D = ixO^LIM^D \}
!           ! {do ix^D = ixO^LIM^D \}
!                time_diff(ix^D,neu,j) = lprim(ix^D, psi_)**2 * bar_E_j(ix^D,neu,j)/sum_div(ix^D,neu,j)
!                if (sum_div(ix^D,neu,j) .eq. 0.0d0) then   !avoid infinity
!                     time_diff(ix^D,neu,j) = 9.9d99
!                write(*,*)  "sim div zero"            
!                endif
!                
!
!            {enddo^D&\}
!       enddo
!      enddo
!
!
!!        do ix1 =4, 10
!!          write(*,*) time_diff(ix1,2,2),"time_diff"   
!!        enddo
!!write(*,*) "################"
!! passs  numerical writing check!! 
!!          write(*,*) bar_E_j(5,3,2) , "bar_E_j"
!!         write(*,*) lprim(5,psi_), "psi"
!!          write(*,*) time_diff(5,3,2)
!!stop 
!!            {do ix^D = ixO^LIM^D \}
!!                write(*,*) time_diff(ix^D,neu,j),div_dEdx(ix^D,1), "time_diff  sum_div"
!!            {enddo^D&\}
!
!!code tets
!!          time_diff(5,neu,2) = -10.d0
!!          time_diff(7,neu,2) = -10.d0
!!          time_diff(10,neu,2) = -10.d0
!!        
!!          time_diff(18,neu,2) = -10.0d0
!!          time_diff(24,neu,2) = -10.0d0
!
!
!
!
!!  region bounded by two sinks --> sink as well
!
!      do j = 1, 2     ! j = 1 is numebr,   j = 2 is energy
!        neu = 1       !  for nue
!            {do ix^D = ixO^LIM^D \}
!                if (time_diff(ix^D,neu,j) .lt. 0.0d0) then
!                     if (lleak_tau(ix^D,neu) .gt. 0.666666d0) then
!                        if (time_diff(ix^D+2,neu,j) .lt. 0.0d0) then
!                              if (lleak_tau(ix^D+1,neu) .gt. 0.666666d0) then
!                                   time_diff(ix^D+1,neu,j) = 9.9d90
!                              endif
!                        else if (time_diff(ix^D+3,neu,j) .lt. 0.0d0) then
!                             do small_region = ix^D, ix^D+2
!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                                   time_diff(small_region,neu,j) = 9.9d90
!                              endif
!                             enddo
!!                        else if (time_diff(ix^D+4,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+3
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+5,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+4
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+6,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+5
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+7,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+7
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+8,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+7
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+9,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+8
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+10,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+9
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+11,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+10
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+12,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+11
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+13,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+12
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+14,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+13
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+15,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+14
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+16,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+15
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+17,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+16
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+18,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+17
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+19,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+18
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+20,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+19
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!!                        else if (time_diff(ix^D+21,neu,j) .lt. 0.0d0) then
!!                             do small_region = ix^D, ix^D+20
!!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!!                                   time_diff(small_region,neu,j) = 9.9d90
!!                              endif
!!                             enddo
!
!
!                        endif
!                      ! optical thin, no need, as if -ve, still to be absolute velus   
!                     endif
!                endif
!            {enddo^D&\}
!
!          neu = 2
!            {do ix^D = ixO^LIM^D \}
!                if (time_diff(ix^D,neu,j) .lt. 0.0d0) then
!                     if (lleak_tau(ix^D,neu) .gt. 0.666666d0) then
!                        if (time_diff(ix^D+2,neu,j) .lt. 0.0d0) then
!                              if (lleak_tau(ix^D+1,neu) .gt. 0.666666d0) then
!                                   time_diff(ix^D+1,neu,j) = 9.9d90
!                              endif
!                        else if (time_diff(ix^D+3,neu,j) .lt. 0.0d0) then
!                             do small_region = ix^D, ix^D+2
!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                                   time_diff(small_region,neu,j) = 9.9d90
!                              endif
!                             enddo
!                   !     else if (time_diff(ix^D+4,neu,j) .lt. 0.0d0) then
!                   !          do small_region = ix^D, ix^D+3
!                   !           if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                   !                time_diff(small_region,neu,j) = 9.9d90
!                   !           endif
!                   !          enddo
!                   !     else if (time_diff(ix^D+5,neu,j) .lt. 0.0d0) then
!                   !          do small_region = ix^D, ix^D+4
!                   !           if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                   !                time_diff(small_region,neu,j) = 9.9d90
!                   !           endif
!                   !          enddo
!                   !     else if (time_diff(ix^D+6,neu,j) .lt. 0.0d0) then
!                   !          do small_region = ix^D, ix^D+5
!                   !           if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                   !                time_diff(small_region,neu,j) = 9.9d90
!                   !           endif
!                   !          enddo
!                   !     else if (time_diff(ix^D+7,neu,j) .lt. 0.0d0) then
!                   !          do small_region = ix^D, ix^D+6
!                   !           if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                   !                time_diff(small_region,neu,j) = 9.9d90
!                   !           endif
!                   !          enddo
!                   !     else if (time_diff(ix^D+8,neu,j) .lt. 0.0d0) then
!                   !          do small_region = ix^D, ix^D+7
!                   !           if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                   !                time_diff(small_region,neu,j) = 9.9d90
!                   !           endif
!                   !          enddo
!                   !     else if (time_diff(ix^D+9,neu,j) .lt. 0.0d0) then
!                   !          do small_region = ix^D, ix^D+8
!                   !           if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                   !                time_diff(small_region,neu,j) = 9.9d90
!                   !           endif
!                   !          enddo
!                        endif
!                     endif
!                endif
!            {enddo^D&\}
!
!
!
!
!           neu = 3
!            {do ix^D = ixO^LIM^D \}
!                if (time_diff(ix^D,neu,j) .lt. 0.0d0) then
!                     if (lleak_tau(ix^D,neu) .gt. 0.666666d0) then
!                        if (time_diff(ix^D+2,neu,j) .lt. 0.0d0) then
!                              if (lleak_tau(ix^D+1,neu) .gt. 0.666666d0) then
!                                   time_diff(ix^D+1,neu,j) = 9.9d90
!                              endif
!                        else if (time_diff(ix^D+3,neu,j) .lt. 0.0d0) then
!                             do small_region = ix^D, ix^D+2
!                              if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                                   time_diff(small_region,neu,j) = 9.9d90
!                              endif
!                             enddo
!                  !      else if (time_diff(ix^D+4,neu,j) .lt. 0.0d0) then
!                  !           do small_region = ix^D, ix^D+3
!                  !            if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                  !                 time_diff(small_region,neu,j) = 9.9d90
!                  !            endif
!                  !           enddo
!                  !      else if (time_diff(ix^D+5,neu,j) .lt. 0.0d0) then
!                  !           do small_region = ix^D, ix^D+4
!                  !            if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                  !                 time_diff(small_region,neu,j) = 9.9d90
!                  !            endif
!                  !           enddo
!                  !      else if (time_diff(ix^D+6,neu,j) .lt. 0.0d0) then
!                  !           do small_region = ix^D, ix^D+5
!                  !            if (lleak_tau(small_region,neu) .gt. 0.666666d0) then 
!                  !                 time_diff(small_region,neu,j) = 9.9d90
!                  !            endif
!                  !           enddo
!                        endif
!                     endif
!                endif
!            {enddo^D&\}
!      enddo
!           
!
!      do j = 1, 2     ! j = 1 is numebr,   j = 2 is energy
!       do neu = 1, 3   
!         ! check if time_diff is -ve   (non-physical)
!            {do ix^D = ixO^LIM^D \}
!                if (time_diff(ix^D,neu,j) .lt. 0.d0) then
!
!                  !         write(*,*)  time_diff(ix^D,neu,j),"time_diff -ve"
!                       if (lleak_tau(ix^D,neu) .gt. 0.666666d0) then
!                           time_diff(ix^D,neu,j) = 9.9d90   !infinity
!               !            write(*,*) time_diff(ix^D,neu,j), ix^D, ".gt. than 0.66666 to be infinity time_diff"
!                       else
!                           time_diff(ix^D,neu,j) = time_diff(ix^D,neu,j) * (-1.0d0)
!                !           write(*,*) time_diff(ix^D,neu,j), ix^D, ".lt. than 0.666666 to be absoulte value outside"
!                       endif
!                endif
!            {enddo^D&\}
!       enddo
!      enddo
!
!                
!      do j = 1, 2     ! j = 1 is numebr,   j = 2 is energy
!       do neu = 1, 3   
!            {do ix^D = ixO^LIM^D \}
!                !  optical thin --> recover free streaming limit
!                       if (lleak_tau(ix^D,neu) .lt. 0.666666d0) then
!                           if (time_diff(ix^D,neu,j) .ge. time_prod(ix^D,neu,j)) then
!                                time_diff(ix^D,neu,j) = 9.9d99   !stop the leakage at this grid
!                           endif
!                       endif
!            {enddo^D&\}
!       enddo
!      enddo
!! -4 ns +8 region, Q_eff = 0.0d0
!!     do j = 1,2
!!        do neu = 1,3
!!                location_i = ns_location1(neu) - 4
!!                location_f = ns_location1(neu) 
!!            time_diff(location_i:location_f,neu,j) = 9.9d99
!!                location_i = ns_location1(neu)
!!                location_f = ns_location1(neu) + 8
!!            time_diff(location_i:location_f,neu,j) = 9.9d99
!!        enddo
!!     enddo
!
!
!!  write(*,*) div_dEdx(5,1), "div_dEdx"
!!code tets
!!            do ix1 = 1, 50
!!                    write(*,*)  time_diff(ix1,1,2), "time_diff" 
!!            enddo
!
!!        {do ix^D = ixI^LIM^D \}
!!                write(*,*) time_diff(ix^D,3,2),  zeta_N(ix^D,1), integral(ix^D,3,2,1), sum_div(ix^D,3,2),"time_diff(nux)  zeta_NN integral "
!!        {enddo^D&\}
!
!                !  pass:  time_prod, lleak_tau, kappa_diff_total, integral, bar_j_E, E_J, dEdx
!                ! not pass: time_diff, sum_div
!                !! alll passss    because   div_E is d(r**2 E)/dr  =  d(r(ix+1)**2E(ix+1) -  d(r(ix)**2 E(ix))
! 
!
!end subroutine compute_diffusion

! need calculate the ghost zones for later integral
subroutine find_kappa_opa_diffusion(ixI^L, ixO^L)
   implicit none
  real*8 get_fermi_integral
  integer, intent(in) :: ixI^L, ixO^L
  !real*8 :: mass_heavy, charge_heavy    ! average ,mass number and charge number  of heavy nuclei
  real*8 :: const
  integer ::  i,j,k,ix^D, m

  double precision, dimension(ixI^S,3,n_eps) :: neu_blocking
  double precision, dimension(ixI^S,2,n_eps) :: e_ae_blocking
! for kappa_opa_diff_total
  double precision, dimension(ixI^S,3,2,n_eps) :: kappa_diff_scat_nucleons
  double precision, dimension(ixI^S,3,2,n_eps) :: kappa_diff_scat_nuclei
  double precision, dimension(ixI^S,2,n_eps) :: kappa_diff_abs
!  double precision, dimension(16) :: eps_bins 
!  double precision, dimension(ixI^S) :: 
!  double precision, dimension(ixI^S) :: 
!  double precision, dimension(ixI^S) :: 

    neu_blocking = 0.0d0
    e_ae_blocking = 0.0d0

    kappa_diff_scat_nucleons = 0.0d0
    kappa_diff_scat_nuclei = 0.0d0
    kappa_diff_abs = 0.0d0
    kappa_diff_total = 0.0d0


   const = ((1.0d0 + 3*g_a**2)/4)*sigma_0


  do m = epsbins_init, epsbins_max
    {do ix^D = ixI^LIM^D \}
!    {do ix^D = ixO^LIM^D \}

!  do m = epsbins_init+1, epsbins_max

!      eta_nu_eq(ix^D,1) = eta_nucleons(ix^D, 1) + eta_p(ix^D) - eta_n(ix^D) - Qnp/lprim(ix^D, temp_) !!!put to tau calculation
!      eta_nu_eq(ix^D,2) =  -eta_nu_eq(ix^D,1)

      neu_blocking(ix^D,1,m)  = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)/lprim(ix^D, temp_) - eta_nu_eq(ix^D,1)  )))
      neu_blocking(ix^D,2,m)  = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)/lprim(ix^D, temp_) - eta_nu_eq(ix^D,2)  )))
  !    neu_blocking(i,j,3,k)  = 1.0d0/(1.0d0 + dexp(-(eps_bins(k)/ltemp(i,j))))



      e_ae_blocking(ix^D,1,m) = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)+Qnp)/lprim(ix^D, temp_) + eta_nucleons(ix^D, 1)))
      e_ae_blocking(ix^D,2,m) = 1.0d0/(1.0d0 + dexp(-(eps_bins(m)-Qnp)/lprim(ix^D, temp_) - eta_nucleons(ix^D, 1)))


!!!!!!!!!!!!!
!   simple leakage   will redo the eta_pn  here  !if not use simple_leakage, pls turn if off
!       eta_pn(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,1)-mass_fraction(ix^D,2))/(dexp(eta_hat(ix^D))-1.0d0)
!       eta_pn(ix^D) = max(eta_pn(ix^D),0.0d0)
!       eta_np(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,2)-mass_fraction(ix^D,1))/(dexp(-eta_hat(ix^D))-1.0d0)
!       eta_np(ix^D) = max(eta_np(ix^D),0.0d0)
!!       if (eta_pn(ix^D) .le. 0.0d0  .or.  eta_np(ix^D) .le. 0.0d0) then
!!           stop "eta_pn, np   -ve"
!!       endif
!
!       if (lprim(ix^D, rho_).lt.1.0d11) then
!          !non degenerate here, use mass fractions as chemical potentials fail at low densities
!          eta_pn(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,2)
!          eta_np(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,1)
!       endif
!!!!!!!!!simple leakage


         {enddo^D&\}
  enddo

!!!!start calculating Diffusion opacity
  do m = epsbins_init, epsbins_max
    {do ix^D = ixI^LIM^D \}
!    {do ix^D = ixO^LIM^D \}


   !!!!!!3neu - absorbed by n and p
          kappa_diff_abs(ix^D,1,m) = const*eta_np(ix^D)*e_ae_blocking(ix^D,1,m)*(eps_bins(m) + Qnp)**2*&
                              (1.0d0 - ((me_mev)/(eps_bins(m) + Qnp))**2)**(1.d0/2.d0)/(me_mev)**2/neu_blocking(ix^D,1,m)

       if (eps_bins(m) .gt. (Qnp + me_mev)) then   

          kappa_diff_abs(ix^D,2,m) = const*eta_pn(ix^D)*e_ae_blocking(ix^D,2,m)*(eps_bins(m) - Qnp)**2*&
                                  (1.0d0 - ((me_mev)/(eps_bins(m) - Qnp))**2)**(1.d0/2.d0)/(me_mev)**2/neu_blocking(ix^D,2,m)

       else
          kappa_diff_abs(ix^D,2,m) = 0.0d0

       endif

   !!!!!!!3neu - nucleons scattering  p,n
      kappa_diff_scat_nucleons(ix^D,1,1,m) = (Cp*sigma_0*ksi_NN(ix^D,1))*(eps_bins(m)/me_mev)**2
      kappa_diff_scat_nucleons(ix^D,2,1,m) = kappa_diff_scat_nucleons(ix^D,1,1,m)
      kappa_diff_scat_nucleons(ix^D,3,1,m) = kappa_diff_scat_nucleons(ix^D,1,1,m)

      kappa_diff_scat_nucleons(ix^D,1,2,m) = (Cn*sigma_0*ksi_NN(ix^D,2))*(eps_bins(m)/me_mev)**2
      kappa_diff_scat_nucleons(ix^D,2,2,m) = kappa_diff_scat_nucleons(ix^D,1,2,m)
      kappa_diff_scat_nucleons(ix^D,3,2,m) = kappa_diff_scat_nucleons(ix^D,1,2,m)

   !!!!!!3neu - nucleis scattering   alpha and heavy

      kappa_diff_scat_nuclei(ix^D,1,1,m) = (1.d0/6.d0)*(16.d0)*((1.d0/2.d0) -1.d0 + (1.d0/2.d0)*(2.d0 - (1.d0/2.d0) - Cv))**2*&
                                            sigma_0*(eps_bins(m)/me_mev)**2*n_N(ix^D,3)

      kappa_diff_scat_nuclei(ix^D,2,1,m) = kappa_diff_scat_nuclei(ix^D,1,1,m)
      kappa_diff_scat_nuclei(ix^D,3,1,m) = kappa_diff_scat_nuclei(ix^D,1,1,m)

                                        !!!!!!!!!charge_heavy and mass_heavy need read table!!!
      kappa_diff_scat_nuclei(ix^D,1,2,m) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*((1.d0/2.d0) -1.d0 +(mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*&
                                            (2.d0 - (1.d0/2.d0) - Cv))**2*sigma_0*(eps_bins(m)/me_mev)**2*n_N(ix^D,4)

      kappa_diff_scat_nuclei(ix^D,2,2,m) = kappa_diff_scat_nuclei(ix^D,1,2,m)
      kappa_diff_scat_nuclei(ix^D,3,2,m) = kappa_diff_scat_nuclei(ix^D,1,2,m)


    !!!!!!do total diff opacity

    !  nue
!code test
      kappa_diff_total(ix^D,1,m) = kappa_diff_abs(ix^D,1,m) + kappa_diff_scat_nucleons(ix^D,1,1,m) +&
!      kappa_diff_total(ix^D,1,m) = kappa_diff_scat_nucleons(ix^D,1,1,m) +&
                                   kappa_diff_scat_nucleons(ix^D,1,2,m) + kappa_diff_scat_nuclei(ix^D,1,1,m) +&
                                   kappa_diff_scat_nuclei(ix^D,1,2,m)
    !  nua
!code test
      kappa_diff_total(ix^D,2,m) = kappa_diff_abs(ix^D,2,m) + kappa_diff_scat_nucleons(ix^D,2,1,m) +&
!      kappa_diff_total(ix^D,2,m) = kappa_diff_scat_nucleons(ix^D,2,1,m) +&
                                   kappa_diff_scat_nucleons(ix^D,2,2,m) + kappa_diff_scat_nuclei(ix^D,2,1,m) +&
                                   kappa_diff_scat_nuclei(ix^D,2,2,m)
    !  nux
      kappa_diff_total(ix^D,3,m) = kappa_diff_scat_nucleons(ix^D,3,1,m) + kappa_diff_scat_nucleons(ix^D,3,2,m) +&
                                   kappa_diff_scat_nuclei(ix^D,3,1,m) + kappa_diff_scat_nuclei(ix^D,3,2,m)
         {enddo^D&\}
      enddo


  ! write(*,*)  kappa_diff_total(5,:,2), "kappa_total_diff"
!checked
!write(*,*) kappa_diff_total(5,1,:),  "kappa_diff_total"
!  write(*,*)  eps_bins(:)
!stop
!passed!
end subroutine find_kappa_opa_diffusion



  subroutine finding_array_size(n^D)
   use mod_global_parameters
   implicit none
!   integer, intent(in)             :: ixI^L, ixO^L
!      ixG^LL --> ixI^L   ;   ixM^LL --> ixO^L
   integer :: nM_local^D, iigrid, igrid, nM^D
   integer :: nG_local^D, nnG^D
   integer, intent(inout) :: n^D
      character*1024 filename


{^IFONED
   nM_local1 = 0
   nG_local1 = 0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       !without ghostzones
!         write(*,*)  ps(igrid)%x(ixMlo1:ixMhi1,r_)
!        write(*,*) '###############'
        nM_local1 = SIZE(ps(igrid)%x(ixMlo1:ixMhi1,1))+nM_local1

 !       write(*,*) size(ps(igrid)%x(ixMlo1:ixMhi1,r_))
!        nG_local1 = SIZE(ps(igrid)%x(ixGlo1:ixGhi1,1))+nG_local1
  
 !       write(*,*) nM_local1 , "nM_local1"
    enddo
        call MPI_ALLREDUCE(nM_local1, nM1, 1, mpi_integer, &
                           MPI_SUM, icomm, ierrmpi)
! write(*,*) '################################ffkkfkfkfkfk#########################'
! write(*,*) '################################ffkkfkfkfkfk#########################'
! write(*,*) '################################ffkkfkfkfkfk#########################'
! write(*,*) '################################ffkkfkfkfkfk#########################'
! write(*,*) '################################ffkkfkfkfkfk#########################'

!        call MPI_ALLREDUCE(nG_local1, nnG1, 1, mpi_integer, &
!                           MPI_SUM, icomm, ierrmpi)
!  write(*,*) nM1,  "nM1"
!   4 inner ghostzones (1:ghostzone1) ,  4 outer ghostzones(n1-ghostzone1 : n1)
!        n1 = nM1 + ghostzones1*2
        n1 = nM1 + ghostzones1

        filename = trim(base_filename)//".n1_number"
    if (mype == 0) then
          open(667,file=filename,status='unknown',position='append')
                      write(667,*) global_time, n1
          close(667) 
    endif

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

    allocate(zeta(n^D,3)) 
    allocate(chi(n^D,3)) 

    allocate(lprim_local(n^D, nprim))
    allocate(lleak_tau_local(n^D, 3))
      
!   ILEAS
    allocate(kappa_opa_total_j(n^D,3,2))
    allocate(kappa_opa_abs_j(n^D,2,2))

    allocate(n_N(n^D,4))
    allocate(zeta_N(n^D, 2))
    allocate(ksi_NN(n^D, 2))
    allocate(n_b(n^D))
    allocate(eta_pn(n^D)) 
    allocate(eta_np(n^D)) 

    allocate(bar_E_j(n^D,3,2))
    allocate(time_prod(n^D,3,2))
    allocate(time_diff(n^D,3,2))
    allocate(kappa_diff_total(n^D,3,n_eps))

    allocate(eta_nu_eq(n^D,3))


!   ILEAS
    kappa_opa_total_j = 0.0d0
    kappa_opa_abs_j = 0.0d0
    n_N = 0.0d0
    zeta_N = 0.0d0 
    ksi_NN = 0.0d0 
    eta_pn = 0.0d0
    eta_np = 0.0d0
    n_b = 0.0d0
    bar_E_j = 0.0d0
    time_prod = 0.0d0
    time_diff = 0.0d0
    kappa_diff_total = 0.0d0
    
    eta_nu_eq = 0.0d0


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

    chi = 0.0d0
    zeta = 0.0d0
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
!   stop "all local x(r_) inside ILEAS "
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
!     stop "all local xi(r_) inside ILEAS "
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
        ! dont map outer ghostzones, as we need the reflective outer boundary
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
!     stop "all local lprim inside ILEAS "

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


! mapping inner ghostzones values |1 |2 |3 |4 | copy real zones |4 |3 |2 |1 | first 4
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
!     stop "all local lprim inside ILEAS "




! mapping outer ghostzones values |93 |94 |95 |96 | copy real zones |96 |95 |94 |93 | last 4
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
!pass
!   write(*,*)  lprim(n1-ghostzones1*2:n1,rho_), "rho_"
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
        if (lprim(i, index_prim) .eq. 0.0d0) then
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

!maybe n1-ghostzones1 to n1  cannot be mapped,  but is okay, as coolingsource in atmosphere region do not matter
!  coolingsource(idir) = momentum  ,  idir+1  =  energy  ,   idir+2  = ye
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
    !deallocate ILEAS variables
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

    deallocate(chi)
    deallocate(zeta)
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


!   ILEAS
    deallocate(kappa_opa_total_j)
    deallocate(kappa_opa_abs_j)

    deallocate(n_N)
    deallocate(n_b)
    deallocate(zeta_N)
    deallocate(ksi_NN)

    deallocate(eta_pn)
    deallocate(eta_np)

    deallocate(eta_nu_eq)

    deallocate(bar_E_j)
    deallocate(time_prod)
    deallocate(time_diff)
    deallocate(kappa_diff_total)
  end subroutine deallocate_arrays

!######################################################################
!!!  simple leakage
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

   

!   write(*,*)  kappa_tot(5,:,1,1)!pass
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



       jm1=1 !energy optical depth, switch to 0 for number
       {do ix^D = ixO^LIM^D \}
          local_eta_nux(ix^D) = 0.0d0   ! (A2)
          ! (A5) equilibrium eta, we have rest masses in our chemical potentials
          ! no need to include mass difference in eta
          eta_nue_eq(ix^D) = eta_nu(ix^D,1)
          ! (A3); note that the ^0 etas are set to 0.0d0
          local_eta_nue(ix^D) =  eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,1))) 
          ! (A4)
          local_eta_nua(ix^D) = -eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,2))) 

          !assuming completely dissociated, valid in side shock
          xlye = lprim(ix^D, ye_)
          ! (A8)
          xyn = (1.0d0-xlye) / (1.0d0 + 2.0d0/3.0d0 * max(eta_nucleons(ix^D,3),0.0d0))
          xyp = xlye / (1.0d0 + 2.0d0/3.0d0*max(eta_nucleons(ix^D,2),0.0d0))
          t1 = dexp(-eta_hat(ix^D))
          ! (A13)
          xynp = max((2.0d0*xlye-1.0d0)/ (t1-1.0d0),0.0d0)
          ! (A14)
          xypn = max(xynp * t1,0.0d0)
 
          ! electron neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nue(ix^D)) / & 
               get_fermi_integral(2+jm1,local_eta_nue(ix^D))
          ! (A15)
          t2 = 1.0d0 + dexp(eta_nucleons(ix^D,1)-get_fermi_integral(5,local_eta_nue(ix^D)) / &
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
          t2 = 1.0d0 + dexp(-eta_nucleons(ix^D,1)-get_fermi_integral(5,local_eta_nua(ix^D)) / &
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
    
       ! compute relative change xerr
       xerr = 0.0d0
       {do ix^D = ixO^LIM^D \}
          xerr = max(xerr,abs(kappa_tot(ix^D,1)/kappa_tot_p(ix^D,1)-1.0d0))
          xerr = max(xerr,abs(kappa_tot(ix^D,2)/kappa_tot_p(ix^D,2)-1.0d0))
          xerr = max(xerr,abs(kappa_tot(ix^D,3)/kappa_tot_p(ix^D,3)-1.0d0))
       {enddo^D&\} ! end loop for every points

       icount = icount + 1

    enddo


    
    if(icount.ge.icount_max) then
       write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
       stop "icount > icount_max in ILEAS; ILEAS.F90"
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

!   write(*,*) eta_nu(5,1), "eta_nu1"
   call compute_zeta_n_N(ixI^L, ixO^L)
!  write(*,*) lprim(5,temp_), lprim(
!  write(*,*)  eta_nu(5,:), lleak_tau(5,:)
!stop "after leaktau"
 
      !  write(*,*) lleak_tau(5,:), lprim(5,temp_),"leak_tau after ruf_tau"
!        stop
  end subroutine find_ruf_tau
!!######################################################################

!ILEAS  
!######################################################################
!  subroutine find_ruf_tau(ixI^L, ixO^L)
!
!    use mod_global_parameters
!
!    implicit none
!    integer, intent(in)             :: ixI^L, ixO^L
!    integer jm1,ix^D, ix_next^D
!    integer icount
!    integer,parameter :: icount_max = 2000
!    ! **** EOS ****
!    double precision, dimension(ixI^S) :: kappa_const_scat_n
!    double precision, dimension(ixI^S) :: kappa_const_scat_p
!!    double precision, dimension(ixI^S) :: kappa_const_scat_h
!    double precision, dimension(ixI^S) :: kappa_const_abs
!    double precision :: csn_0,csp_0,t1,t2
!    double precision :: xerr
!    double precision, parameter :: xerr_out = 1.0d-12
!!    double precision, dimension(ixI^S,3) :: kappa_tot   
!    double precision, dimension(ixI^S,3,2) :: kappa_tot_p 
!    double precision, dimension(ixI^S,3) :: kappa_scat_n ! 1/cm
!    double precision, dimension(ixI^S,3) :: kappa_scat_p ! 1/cm
!!    double precision, dimension(ixI^S) :: kappa_scat_h ! 1/cm
!    double precision, dimension(ixI^S) :: kappa_abs_n  ! 1/cm
!    double precision, dimension(ixI^S) :: kappa_abs_p  ! 1/cm
!    double precision, dimension(ixI^S) :: local_eta_nue
!    double precision, dimension(ixI^S) :: local_eta_nua
!    double precision, dimension(ixI^S) :: local_eta_nux
!    double precision, dimension(ixI^S) :: eta_nue_eq
!
!
!    character*1024:: filename
!
!    double precision :: dr
!    double precision :: xlye,xyn,xynp,xyp,xypn
!    double precision :: get_fermi_integral
!
!    ! initialize some suff:
!    kappa_opa_total_j = 1.0d0
!    kappa_tot_p = 1.0d0
!    kappa_scat_n = 1.0d-5 ! 1/cm
!    kappa_scat_p = 1.0d-5 ! 1/cm
!    kappa_abs_n  = 1.0d-5 ! 1/cm
!    kappa_abs_p  = 1.0d-5 ! 1/cm
!
!    local_eta_nux = 0.0d0
!    local_eta_nue = 0.0d0
!    local_eta_nua = 0.0d0
!
!    ! Neutrino-nucleon scattering transport cross-section
!    
!    ! C_s,N in equation (A1) of Ruffert et al.
!    !
!    csn_0 = (1.0d0 + 5.0D0*alpha**2) / 24.0d0
!    csp_0 = (4.0d0*(Cv-1.0d0)**2 + 5.0d0*alpha**2) / 24.0d0
!    !
!!write(*,*) ixO^LIM1
!!stop 
!
!    {do ix^D = ixI^LIM^D \}
!   ! {do ix^D = ixO^LIM^D \}
!       ! constant parts of kappa (A6)
!       t1 = sigma_0 * avo * lprim(ix^D, rho_) * (lprim(ix^D, temp_)/me_mev)**2
!       kappa_const_scat_n(ix^D) = csn_0 * t1
!       kappa_const_scat_p(ix^D) = csp_0 * t1
!       ! (A11) constant part
!       kappa_const_abs(ix^D) = (1.0d0+3.0d0*alpha**2)/4.0d0 * t1
!    {enddo^D&\} ! end loop for every points
!
!    ! Loop to get converged result for tau.
!    ! This is discussed in the text between equations (A5) and
!    ! (A6). Note that for the initial iteration the kappas are set to 1.0d-5
!    ! unless we have tau from previous time, then use it as starting point
!    icount = 1
!    xerr = 1.0d0
!
!! initial values of kappa_opa_total_j (j=2)
!       kappa_opa_total_j(ixI^S,1,2) = &
!              kappa_scat_p(ixI^S,1) &
!            + kappa_scat_n(ixI^S,1) &
!            + kappa_abs_n(ixI^S)
!       ! antis; (A18)
!       kappa_opa_total_j(ixI^S,2,2) = &
!              kappa_scat_p(ixI^S,2) &
!            + kappa_scat_n(ixI^S,2) &
!            + kappa_abs_p(ixI^S)
!
!       ! nu_xl (A19)
!       kappa_opa_total_j(ixI^S,3,2) = &
!            + kappa_scat_p(ixI^S,3) &
!            + kappa_scat_n(ixI^S,3)
!
!!code test     tried    same result..
!!   {do ix^D = ixO^LIM^D \}
!!       lleak_tau(ix^D,:) =  (10**(11.d0)/lprim(ix^D,rho_))**2
!!   {enddo^D&\} ! end loop for every points
!
!
!    do while(xerr.gt.xerr_out .and. icount.lt.icount_max)
!       ! copy over current into previous kappa
!       ! they have same array size
!     write(*,*) icount,'iconut'
!       ! compute relative change xerr
!       xerr = 0.0d0
!       {do ix^D = ixI^LIM^D \}
!       !{do ix^D = ixO^LIM^D \}
!          xerr = max(xerr,abs(kappa_opa_total_j(ix^D,1,2)/kappa_tot_p(ix^D,1,2)-1.0d0))
!          xerr = max(xerr,abs(kappa_opa_total_j(ix^D,2,2)/kappa_tot_p(ix^D,2,2)-1.0d0))
!          xerr = max(xerr,abs(kappa_opa_total_j(ix^D,3,2)/kappa_tot_p(ix^D,3,2)-1.0d0))
!       {enddo^D&\} ! end loop for every points
!
!       icount = icount + 1
!       kappa_tot_p = kappa_opa_total_j
!
!!fixme when 2D
!
!    if(icount.gt.2) then
!
!         lleak_tau = 0.0d0
!         {do ix^D = ixImax^D-1, ixImin^D, -1 \}
!        ! {do ix^D = ixOmax^D, ixOmin^D, -1 \}
!               dr = lxi(ix1+1,r_) - lxi(ix1,r_)
!               ix_next^D = ix^D + kr(1,^D)
!               lleak_tau(ix^D, 1:3) = lleak_tau(ix_next^D,1:3) + &
!                     kappa_opa_total_j(ix^D,1:3,2) * dr * lprim(ix^D, psi_)**2
!         {enddo^D&\} ! end loop for every points
!    endif
!
!
!
!       {do ix^D = ixI^LIM^D \}
!      ! {do ix^D = ixO^LIM^D \}
!          local_eta_nux(ix^D) = 0.0d0   ! (A2)
!
!          eta_nue_eq(ix^D) = eta_nu(ix^D,1)
!
!          local_eta_nue(ix^D) =  eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,1))) 
!          local_eta_nua(ix^D) = -eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,2))) 
!       {enddo^D&\} ! end loop for every points
!
!    
!
!       call find_kappa_opa_total(local_eta_nue, local_eta_nua, ixI^L, ixO^L)
!
!    enddo
!
!
!    
!    if(icount.ge.icount_max) then
!       write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
!       stop "icount > icount_max in ILEAS; ILEAS.F90"
!    endif
!
!!fixme when 2D
!    ! Recompute tau based on the most recent kappa_tot
!         lleak_tau = 0.0d0
!         {do ix^D = ixImax^D-1, ixImin^D, -1 \}
!         !{do ix^D = ixOmax^D, ixOmin^D, -1 \}
!               dr = lxi(ix1+1,r_) - lxi(ix1,r_)
!               ix_next^D = ix^D + kr(1,^D)
!               lleak_tau(ix^D, 1:3) = lleak_tau(ix_next^D,1:3) + &
!                     kappa_opa_total_j(ix^D,1:3,2) * dr * lprim(ix^D, psi_)**2
!         {enddo^D&\} ! end loop for every points
!!         call cal_leak_tau(ixI^L, ixO^L, kappa_tot, lx, lxi, leak_tau, lprim)
!
!
! !   first_iteration = .false.    
!!###################  i took equatorial  ns_lcoation  for 2D
!
!
!  ns_location^D(1:3) = ixOmin^D
!
!  {do ix^D = ixI^LIM^D \}
!  !{do ix^D = ixO^LIM^D \}
!       !nu-spheres for heating          ns  assumes  shperical symmetric for all
!       !                                dimensions  (ns_location_i, 1, 1)
!       if (lleak_tau(ix^D,1).gt.0.66666d0) then
!          ns_location1(1) = ix1
!       endif
!
!       if (lleak_tau(ix^D,2).gt.0.66666d0) then
!          ns_location1(2) = ix1
!       endif
!
!       if (lleak_tau(ix^D,3).gt.0.66666d0) then
!          ns_location1(3) = ix1
!       endif
!  {enddo^D&\} 
!
!    heat_erms(1) = lprim(ns_location^D(1), temp_)*&
!                   sqrt(get_fermi_integral(5,local_eta_nue(ns_location^D(1)))/&
!                   get_fermi_integral(3,local_eta_nue(ns_location^D(1))))
!
!    heat_erms(2) = lprim(ns_location^D(2), temp_)*&
!                   sqrt(get_fermi_integral(5,local_eta_nua(ns_location^D(2)))/&
!                   get_fermi_integral(3,local_eta_nua(ns_location^D(2))))
!
!    heat_erms(3) = lprim(ns_location^D(3), temp_)*&
!                   sqrt(get_fermi_integral(5,0.0d0)/ &
!                   get_fermi_integral(3,0.0d0))
!
!    heat_em(1) = lprim(ns_location^D(1), temp_)*get_fermi_integral(5,local_eta_nue(ns_location^D(1)))/&
!         get_fermi_integral(4,local_eta_nue(ns_location^D(1)))
!    heat_em(2) = lprim(ns_location^D(2), temp_)*get_fermi_integral(5,local_eta_nua(ns_location^D(2)))/&
!         get_fermi_integral(4,local_eta_nua(ns_location^D(2)))
!    heat_em(3) = lprim(ns_location^D(3), temp_)*get_fermi_integral(5,0.0d0)/ &
!         get_fermi_integral(4,0.0d0) !not used
!
!    !set degeneracy factors to interpolated values
!    eta_nu(ixI^S,1) = local_eta_nue(ixI^S)
!    eta_nu(ixI^S,2) = local_eta_nua(ixI^S)
!    eta_nu(ixI^S,3) = local_eta_nux(ixI^S)
!    
!!  write(*,*) lprim(5,temp_), lprim(
!!  write(*,*)  eta_nu(5,:), lleak_tau(5,:)
!!stop "after leaktau"
! 
! !  call compute_zeta_n_N(ixI^L, ixO^L)
!
!!  lleak_tau(4) and lleak_tau(5)   are the same
!!passed
!!        write(*,*) lleak_tau(1:8,1),"leak_tau after ruf_tau"
!!        write(*,*) eta_nu(1:8,1),"eta_nu after ruf_tau"
!
!!  notpass!!! 
!   
!  end subroutine find_ruf_tau


!ILEAS   !density function 
!!######################################################################
!  subroutine find_ruf_tau(ixI^L, ixO^L)
!
!    use mod_global_parameters
!
!    implicit none
!    integer, intent(in)             :: ixI^L, ixO^L
!    integer jm1,ix^D, ix_next^D
!    integer icount
!    integer,parameter :: icount_max = 2000
!    ! **** EOS ****
!    double precision, dimension(ixI^S) :: kappa_const_scat_n
!    double precision, dimension(ixI^S) :: kappa_const_scat_p
!!    double precision, dimension(ixI^S) :: kappa_const_scat_h
!    double precision, dimension(ixI^S) :: kappa_const_abs
!    double precision :: csn_0,csp_0,t1,t2
!    double precision :: xerr
!    double precision, parameter :: xerr_out = 1.0d-12
!!    double precision, dimension(ixI^S,3) :: kappa_tot   
!    double precision, dimension(ixI^S,3,2) :: kappa_tot_p 
!    double precision, dimension(ixI^S,3) :: kappa_scat_n ! 1/cm
!    double precision, dimension(ixI^S,3) :: kappa_scat_p ! 1/cm
!!    double precision, dimension(ixI^S) :: kappa_scat_h ! 1/cm
!    double precision, dimension(ixI^S) :: kappa_abs_n  ! 1/cm
!    double precision, dimension(ixI^S) :: kappa_abs_p  ! 1/cm
!    double precision, dimension(ixI^S) :: local_eta_nue
!    double precision, dimension(ixI^S) :: local_eta_nua
!    double precision, dimension(ixI^S) :: local_eta_nux
!    double precision, dimension(ixI^S) :: eta_nue_eq
!
!
!    character*1024:: filename
!
!    double precision :: dr
!    double precision :: xlye,xyn,xynp,xyp,xypn
!    double precision :: get_fermi_integral
!
!    ! initialize some suff:
!!    kappa_opa_total_j = 1.0d0
!    kappa_tot_p = 1.0d0
!!    kappa_scat_n = 1.0d-5 ! 1/cm
!!    kappa_scat_p = 1.0d-5 ! 1/cm
!!    kappa_abs_n  = 1.0d-5 ! 1/cm
!!    kappa_abs_p  = 1.0d-5 ! 1/cm
!
!    local_eta_nux = 0.0d0
!    local_eta_nue = 0.0d0
!    local_eta_nua = 0.0d0
!
!    ! Neutrino-nucleon scattering transport cross-section
!    
!    ! C_s,N in equation (A1) of Ruffert et al.
!    !
!    csn_0 = (1.0d0 + 5.0D0*alpha**2) / 24.0d0
!    csp_0 = (4.0d0*(Cv-1.0d0)**2 + 5.0d0*alpha**2) / 24.0d0
!    !
!!write(*,*) ixO^LIM1
!!stop 
!
!       {do ix^D = ixO^LIM^D \}
!           lleak_tau(ix^D,:) =  (10**(11.d0)/lprim(ix^D,rho_))**2
!       {enddo^D&\} ! end loop for every points
!
!       {do ix^D = ixI^LIM^D \}
!      ! {do ix^D = ixO^LIM^D \}
!          local_eta_nux(ix^D) = 0.0d0   ! (A2)
!          eta_nue_eq(ix^D) = eta_nu(ix^D,1)
!          local_eta_nue(ix^D) =  eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,1))) 
!          local_eta_nua(ix^D) = -eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,2))) 
!       {enddo^D&\} ! end loop for every points
!
!       call find_kappa_opa_total(local_eta_nue, local_eta_nua, ixI^L, ixO^L)
!
!    icount = 1
!    xerr = 1.0d0
!
!
!
!    do while(xerr.gt.xerr_out .and. icount.lt.icount_max)
!       ! copy over current into previous kappa
!       ! they have same array size
!       !write(*,*) kappa_opa_total_j(5,1,2), kappa_tot_p(5,1,2)
!
!
!       ! write(*,*)  icount, 'incount'
!
!!fixme when 2D
!    if(icount.ge.2) then
!         lleak_tau = 0.0d0
!         {do ix^D = ixImax^D-1, ixImin^D, -1 \}
!        ! {do ix^D = ixOmax^D, ixOmin^D, -1 \}
!               dr = lxi(ix1+1,r_) - lxi(ix1,r_)
!               ix_next^D = ix^D + kr(1,^D)
!               lleak_tau(ix^D, 1:3) = lleak_tau(ix_next^D,1:3) + &
!                     kappa_opa_total_j(ix^D,1:3,2) * dr * lprim(ix^D, psi_)**2
!         {enddo^D&\} ! end loop for every points
!
!
!       {do ix^D = ixI^LIM^D \}
!      ! {do ix^D = ixO^LIM^D \}
!          local_eta_nux(ix^D) = 0.0d0   ! (A2)
!          eta_nue_eq(ix^D) = eta_nu(ix^D,1)
!          local_eta_nue(ix^D) =  eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,1))) 
!          local_eta_nua(ix^D) = -eta_nue_eq(ix^D) * (1.0d0-dexp(-lleak_tau(ix^D,2))) 
!       {enddo^D&\} ! end loop for every points
!
!       call find_kappa_opa_total(local_eta_nue, local_eta_nua, ixI^L, ixO^L)
!    endif
!
!
!
!       ! compute relative change xerr
!       xerr = 0.0d0
!!       {do ix^D = ixI^LIM^D \}
!       {do ix^D = ixO^LIM^D \}
!          xerr = max(xerr,abs(kappa_opa_total_j(ix^D,1,2)/kappa_tot_p(ix^D,1,2)-1.0d0))
!          xerr = max(xerr,abs(kappa_opa_total_j(ix^D,2,2)/kappa_tot_p(ix^D,2,2)-1.0d0))
!          xerr = max(xerr,abs(kappa_opa_total_j(ix^D,3,2)/kappa_tot_p(ix^D,3,2)-1.0d0))
!       {enddo^D&\} ! end loop for every points
!
!       icount = icount + 1
!       kappa_tot_p = kappa_opa_total_j
!
!    enddo
!
!
!    
!    if(icount.ge.icount_max) then
!       write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
!       stop "icount > icount_max in ILEAS; ILEAS.F90"
!    endif
!
!!fixme when 2D
!    ! Recompute tau based on the most recent kappa_tot
!         lleak_tau = 0.0d0
!         {do ix^D = ixImax^D-1, ixImin^D, -1 \}
!         !{do ix^D = ixOmax^D, ixOmin^D, -1 \}
!               dr = lxi(ix1+1,r_) - lxi(ix1,r_)
!               ix_next^D = ix^D + kr(1,^D)
!               lleak_tau(ix^D, 1:3) = lleak_tau(ix_next^D,1:3) + &
!                     kappa_opa_total_j(ix^D,1:3,2) * dr * lprim(ix^D, psi_)**2
!         {enddo^D&\} ! end loop for every points
!!         call cal_leak_tau(ixI^L, ixO^L, kappa_tot, lx, lxi, leak_tau, lprim)
!
!
! !   first_iteration = .false.    
!!###################  i took equatorial  ns_lcoation  for 2D
!
!
!  ns_location^D(1:3) = 1+ghostzones1
!
!  {do ix^D = ixI^LIM^D \}
!  !{do ix^D = ixO^LIM^D \}
!       !nu-spheres for heating          ns  assumes  shperical symmetric for all
!       !                                dimensions  (ns_location_i, 1, 1)
!       if (lleak_tau(ix^D,1).gt.0.66666d0) then
!          ns_location1(1) = ix1
!       endif
!
!       if (lleak_tau(ix^D,2).gt.0.66666d0) then
!          ns_location1(2) = ix1
!       endif
!
!       if (lleak_tau(ix^D,3).gt.0.66666d0) then
!          ns_location1(3) = ix1
!       endif
!  {enddo^D&\} 
!
!    heat_erms(1) = lprim(ns_location^D(1), temp_)*&
!                   sqrt(get_fermi_integral(5,local_eta_nue(ns_location^D(1)))/&
!                   get_fermi_integral(3,local_eta_nue(ns_location^D(1))))
!
!    heat_erms(2) = lprim(ns_location^D(2), temp_)*&
!                   sqrt(get_fermi_integral(5,local_eta_nua(ns_location^D(2)))/&
!                   get_fermi_integral(3,local_eta_nua(ns_location^D(2))))
!
!    heat_erms(3) = lprim(ns_location^D(3), temp_)*&
!                   sqrt(get_fermi_integral(5,0.0d0)/ &
!                   get_fermi_integral(3,0.0d0))
!
!    heat_em(1) = lprim(ns_location^D(1), temp_)*get_fermi_integral(5,local_eta_nue(ns_location^D(1)))/&
!         get_fermi_integral(4,local_eta_nue(ns_location^D(1)))
!    heat_em(2) = lprim(ns_location^D(2), temp_)*get_fermi_integral(5,local_eta_nua(ns_location^D(2)))/&
!         get_fermi_integral(4,local_eta_nua(ns_location^D(2)))
!    heat_em(3) = lprim(ns_location^D(3), temp_)*get_fermi_integral(5,0.0d0)/ &
!         get_fermi_integral(4,0.0d0) !not used
!
!    !set degeneracy factors to interpolated values
!    eta_nu(ixI^S,1) = local_eta_nue(ixI^S)
!    eta_nu(ixI^S,2) = local_eta_nua(ixI^S)
!    eta_nu(ixI^S,3) = local_eta_nux(ixI^S)
!    
!!  write(*,*) lprim(5,temp_), lprim(
!!  write(*,*)  eta_nu(5,:), lleak_tau(5,:)
!!stop "after leaktau"
! 
! !  call compute_zeta_n_N(ixI^L, ixO^L)
!
!!  lleak_tau(4) and lleak_tau(5)   are the same
!!passed
!!        write(*,*) lleak_tau(1:8,1),"leak_tau after ruf_tau"
!!        write(*,*) eta_nu(1:8,1),"eta_nu after ruf_tau"
!
!!  notpass!!! 
!   
!  end subroutine find_ruf_tau


  
!!!!!!!!!!!Find average absorption opa and scattering opa for absorption and optical depth!!!!!!!
subroutine find_kappa_opa_total(local_eta_nue, local_eta_nua, ixI^L, ixO^L)
  implicit none

    integer, intent(in)             :: ixI^L, ixO^L
    integer jm1,ix^D, ix_next^D
  character(len=100) filename
   !integer,parameter :: icount_max = 200, jcount_max = 200


   ! real*8 :: kappa_total(n1,n2,3,2)    ! <rad:ang:neu:j>
   ! real*8 :: kappa_nu_abs(n1,n2,2,2)   !   <rad:ang:neu (nu_e and anu_e):j>
   ! real*8 :: kappa_scat(n1,n2,3,4,2) !scattering for nucleons and nucleis for all types neu   <rad:ang:neu:nucleons or nucleis:j>

   ! real*8 :: local_T
    double precision, intent(in) :: local_eta_nue(ixI^S)     ! eta_nu for all 3 types
    double precision, intent(in) :: local_eta_nua(ixI^S)           !!!!eta_ae = - eta_e
   ! real*8,parameter :: tor = 10**(-7)
    real*8 get_fermi_integral
    real*8 :: const, xyn, xyp
!    real*8 :: n_N(n1,n2,4)               !number den of given nucleon type or nuclei type 1=p,2=n,3=alpha,4=heavy
!    real*8 :: E_fermi(n1,n2,2)           !Fermi energy for p,n      1=p, 2=n
   ! kappa_total(:,:,:,:) = 1.0d0
   ! kappa_nu_abs(:,:,:,:) = 1.0d-5
   ! kappa_scat(:,:,:,:) = 1.0d-5

   ! local_eta_nu(:,:,:) = 0.0d0
   ! local_eta_e(:,:) = 0.0d0
    double precision, dimension(ixI^S,3,2,2) :: kappa_opa_scat_nucleons_j
    double precision, dimension(ixI^S,3,2,2) :: kappa_opa_scat_nuclei_j
    double precision, dimension(ixI^S,2) :: lep_blocking
    double precision, dimension(ixI^S,2) :: deno
    double precision, dimension(ixI^S,2) :: deno2



   kappa_opa_total_j = 0.0d0
   kappa_opa_scat_nucleons_j = 0.0d0
   kappa_opa_scat_nuclei_j = 0.0d0
!  needed to use in other sub
   kappa_opa_abs_j = 0.0d0
   lep_blocking = 0.0d0
   deno = 0.0d0
   deno2 = 0.0d0

   call compute_zeta_n_N(ixI^L, ixO^L)


   const =  ((1.d0+3.d0*g_A**2)/(4.d0*(me_mev)**2))*sigma_0

  !pass

  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
!     if (mu_no_rest_mass) then
!      eta_np(ix^D) = avo*lrho(ix^D)*(lyp(ix^D) - lyn(ix^D))/(dexp(eta_p(ix^D) - eta_n(ix^D) ) - 1.0d0)     !B4
!      eta_pn(ix^D) = avo*lrho(ix^D)*(lyp(ix^D) - lyn(ix^D))/(-dexp(eta_n(ix^D) - eta_p(ix^D) ) + 1.0d0)    !B5
!     else
!      eta_np(ix^D) = avo*lrho(ix^D)*(lyp(ix^D) - lyn(ix^D))/&
!                      (dexp(eta_p(ix^D) - eta_n(ix^D) + Qnp/lprim(ix^D, temp_)) - 1.0d0)     !B4
!      eta_pn(ix^D) = avo*lrho(ix^D)*(lyp(ix^D) - lyn(ix^D))/&
!                      (-dexp(eta_n(ix^D) - eta_p(ix^D) - Qnp/lprim(ix^D, temp_)) + 1.0d0)    !B5
!     endif
!
!
!      eta_np(ix^D) = max(eta_np(ix^D), 0.0d0)
!      eta_pn(ix^D) = max(eta_pn(ix^D), 0.0d0)
!
!      if (lrho(ix^D).lt.1.0d11) then
!               !non degenerate here,    chem potential fails at low density
!          eta_pn(ix^D) = avo*lrho(ix^D)* lyp(ix^D)
!          eta_np(ix^D) = avo*lrho(ix^D)* lyn(ix^D)
!      endif

!pass
      !has rest_mass diff correction      should i erase the Qnp for eta_nue???????/////////
      lep_blocking(ix^D,1) = 1.0d0/( 1.d0 + dexp(-((lprim(ix^D, temp_)*get_fermi_integral(5,local_eta_nue(ix^D))/&
                              get_fermi_integral(4,local_eta_nue(ix^D)) + Qnp )/lprim(ix^D, temp_) - eta_nucleons(ix^D,1))))      !!!!! for  e-  meanenergy
      lep_blocking(ix^D,2) = 1.0d0/( 1.d0 + dexp(-(get_fermi_integral(5,local_eta_nua(ix^D))/&
                              get_fermi_integral(4,local_eta_nua(ix^D)) + eta_nucleons(ix^D,1)))) !!!! for e+  mean energy

            if (lep_blocking(ix^D,1) .le. 0.0d0 .or. lep_blocking(ix^D,2) .le. 0.0d0) then
                 write(*,*)  lep_blocking(ix^D,1), lep_blocking(ix^D,2), ix^D
                write(*,*) "-ve lepton blocking"
                !stop  "lep_block  -ve "
            endif

  {enddo^D&\} 

!   write(*,*) 1/lep_blocking(5,2) , local_eta_nua(5), eta_nucleons(5,1),"lep_blocking  local_eta_nua eta_nucleosn"

!!!!!!!! bar_opa for absorption  for nue and nua only
  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
         !!!!!nue
           deno(ix^D,1) = (lprim(ix^D, temp_)**2*get_fermi_integral(4,local_eta_nue(ix^D)) +&
                            2.d0*Qnp*lprim(ix^D, temp_)*get_fermi_integral(3,local_eta_nue(ix^D)) +&
                            Qnp**2*get_fermi_integral(2,local_eta_nue(ix^D)))/get_fermi_integral(2,local_eta_nue(ix^D))
!!j = 0
           deno(ix^D,2) = (lprim(ix^D, temp_)**2*get_fermi_integral(5,local_eta_nue(ix^D)) +&
                            2.d0*Qnp*lprim(ix^D, temp_)*get_fermi_integral(4,local_eta_nue(ix^D)) +&
                            Qnp**2*get_fermi_integral(3,local_eta_nue(ix^D)))/get_fermi_integral(3,local_eta_nue(ix^D))
            if (deno(ix^D,1) .le. 0.0d0 .or. deno(ix^D,2) .le. 0.0d0) then
              write(*,*)    deno(ix^D,1), deno(ix^D,2)
            stop  "deno1  is -ve "
            endif
         !!j = 1
           kappa_opa_abs_j(ix^D,1,1) = const*eta_np(ix^D)*lep_blocking(ix^D,1)*deno(ix^D,1)   !!!(B13)  j = 0

           kappa_opa_abs_j(ix^D,1,2) = const*eta_np(ix^D)*lep_blocking(ix^D,1)*deno(ix^D,2)   !!!(B13)  j = 1

         !!!!!nua
                 !!! j = 0
           deno2(ix^D,1) = ((lprim(ix^D, temp_)**2*get_fermi_integral(4,local_eta_nua(ix^D) -&
                            Qnp/lprim(ix^D, temp_)) + (2.d0)*Qnp*lprim(ix^D, temp_)*get_fermi_integral(3,local_eta_nua(ix^D)-&
                            Qnp/lprim(ix^D, temp_)))/get_fermi_integral(2,local_eta_nua(ix^D))) +&
                           (1.d0*Qnp**2*get_fermi_integral(2,local_eta_nua(ix^D) -&
                            Qnp/lprim(ix^D, temp_))/get_fermi_integral(2,local_eta_nua(ix^D)))

                  !!! j = 1
           deno2(ix^D,2) = ((lprim(ix^D, temp_)**2*get_fermi_integral(5,local_eta_nua(ix^D) -&
                            Qnp/lprim(ix^D, temp_)) + (2.d0 + 1.d0)*Qnp*lprim(ix^D, temp_)*&
                            get_fermi_integral(4,local_eta_nua(ix^D)-&
                            Qnp/lprim(ix^D, temp_)))/get_fermi_integral(3,local_eta_nua(ix^D))) +&
                           (3.d0*Qnp**2*get_fermi_integral(3,local_eta_nua(ix^D) -&
                            Qnp/lprim(ix^D, temp_))/get_fermi_integral(3,local_eta_nua(ix^D))) +&
                           ((Qnp**3/lprim(ix^D, temp_))*get_fermi_integral(2,local_eta_nua(ix^D) -&
                            Qnp/lprim(ix^D, temp_))/get_fermi_integral(3,local_eta_nua(ix^D)))
!pass
            if (deno2(ix^D,1) .le. 0.0d0 .or. deno2(ix^D,2) .le. 0.0d0) then
             write(*,*)    deno2(ix^D,1), deno2(ix^D,2)
            stop  "deno2  is -ve "
            endif

           kappa_opa_abs_j(ix^D,2,1) = const*eta_pn(ix^D)*lep_blocking(ix^D,2)*deno2(ix^D,1)   !!!B14  j = 0


           kappa_opa_abs_j(ix^D,2,2) = const*eta_pn(ix^D)*lep_blocking(ix^D,2)*deno2(ix^D,2)   !!!B14  j = 1


  {enddo^D&\} 
!passs!!!!
!write(*,*)  kappa_opa_abs_j(5,:,2) , "kappa_opa_abs_j   for pront and neu"
!      write(*,*)  "kappa_opa_abs_j"
!       write(*,*)  kappa_opa_abs_j(5,1,1,:,1), "kappa_opa_abs"
!       write(*,*)  kappa_opa_abs_j(5,1,1,:,2)






  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
!!!!j = 0  Nucleons 3 neu
   !!   with p
         !!!j = 0
         kappa_opa_scat_nucleons_j(ix^D,1,1,1) =  Cp*sigma_0*ksi_NN(ix^D,1)*&
                                                   (lprim(ix^D, temp_)/me_mev)**2*get_fermi_integral(4,local_eta_nue(ix^D))/&
                                                   get_fermi_integral(2,local_eta_nue(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,2,1,1) =  Cp*sigma_0*ksi_NN(ix^D,1)*& 
                                                   (lprim(ix^D, temp_)/me_mev)**2*get_fermi_integral(4,local_eta_nua(ix^D))/&
                                                   get_fermi_integral(2,local_eta_nua(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,3,1,1) =  Cp*sigma_0*ksi_NN(ix^D,1)*&
                                                   (lprim(ix^D, temp_)/me_mev)**2*get_fermi_integral(4,0.0d0)/&
                                                   get_fermi_integral(2,0.0d0)

   !!   with n

         kappa_opa_scat_nucleons_j(ix^D,1,2,1) =  Cn*sigma_0*ksi_NN(ix^D,2)*&
                                                   (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(4,local_eta_nue(ix^D))/get_fermi_integral(2,local_eta_nue(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,2,2,1) =  Cn*sigma_0*ksi_NN(ix^D,2)*&
                                                   (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(4,local_eta_nua(ix^D))/get_fermi_integral(2,local_eta_nua(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,3,2,1) =  Cn*sigma_0*ksi_NN(ix^D,2)*&
                                                  (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(4,0.0d0)/get_fermi_integral(2,0.0d0)


!!!!!!!!!!!!!!!j = 1 Nucleons 3 neu
  !!   with p    !!! j = 1
         kappa_opa_scat_nucleons_j(ix^D,1,1,2) =  Cp*sigma_0*ksi_NN(ix^D,1)*&
                                                   (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(5,local_eta_nue(ix^D))/get_fermi_integral(3,local_eta_nue(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,2,1,2) =  Cp*sigma_0*ksi_NN(ix^D,1)*&
                                                  (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(5,local_eta_nua(ix^D))/get_fermi_integral(3,local_eta_nua(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,3,1,2) =  Cp*sigma_0*ksi_NN(ix^D,1)*&
                                                  (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(5,0.0d0)/get_fermi_integral(3,0.0d0)

  !!  with n

         kappa_opa_scat_nucleons_j(ix^D,1,2,2) =  Cn*sigma_0*ksi_NN(ix^D,2)*&
                                                  (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(5,local_eta_nue(ix^D))/get_fermi_integral(3,local_eta_nue(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,2,2,2) =  Cn*sigma_0*ksi_NN(ix^D,2)*&
                                                  (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(5,local_eta_nua(ix^D))/get_fermi_integral(3,local_eta_nua(ix^D))

         kappa_opa_scat_nucleons_j(ix^D,3,2,2) =  Cn*sigma_0*ksi_NN(ix^D,2)*&
                                                  (lprim(ix^D, temp_)/me_mev)**2*&
                                                   get_fermi_integral(5,0.0d0)/get_fermi_integral(3,0.0d0)

!!!!!!!!!!!!!!!!!!!!! j = 0  Nuclei 3 neu
  !!  with alpha

         kappa_opa_scat_nuclei_j(ix^D,1,1,1) = (1.d0/6.d0)*16.d0*(Ca - 1.0d0 + (0.5d0)*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,3)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(4,local_eta_nue(ix^D))/get_fermi_integral(2,local_eta_nue(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,2,1,1) = (1.d0/6.d0)*16.d0*(Ca - 1.0d0 + (0.5d0)*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,3)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(4,local_eta_nua(ix^D))/get_fermi_integral(2,local_eta_nua(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,3,1,1) = (1.d0/6.d0)*16.d0*(Ca - 1.0d0 + (0.5d0)*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,3)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(4,0.0d0)/get_fermi_integral(2,0.0d0)


  !!  with heavy

         kappa_opa_scat_nuclei_j(ix^D,1,2,1) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*(Ca - 1.0d0 + &
                                                (mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*&
                                                (2.d0 - Ca - Cv))**2*sigma_0*n_N(ix^D,4)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(4,local_eta_nue(ix^D))/get_fermi_integral(2,local_eta_nue(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,2,2,1) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*(Ca - 1.0d0 +  &
                                                (mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*&
                                                (2.d0 - Ca - Cv))**2*sigma_0*n_N(ix^D,4)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(4,local_eta_nua(ix^D))/get_fermi_integral(2,local_eta_nua(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,3,2,1) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*&
                                                (Ca - 1.0d0 + (mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,4)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(4,0.0d0)/get_fermi_integral(2,0.0d0)


!!! j = 1 Nuclei 3 neu
  !!  with alpha

         kappa_opa_scat_nuclei_j(ix^D,1,1,2) = (1.d0/6.d0)*16.d0*(Ca - 1.0d0 + (0.5d0)*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,3)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(5,local_eta_nue(ix^D))/&
                                                get_fermi_integral(3,local_eta_nue(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,2,1,2) = (1.d0/6.d0)*16.d0*(Ca - 1.0d0 + (0.5d0)*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,3)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(5,local_eta_nua(ix^D))/&
                                                get_fermi_integral(3,local_eta_nua(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,3,1,2) = (1.d0/6.d0)*16.d0*(Ca - 1.0d0 + (0.5d0)*(2.d0 - Ca - Cv))**2*&
                                                sigma_0*n_N(ix^D,3)*(lprim(ix^D, temp_)/me_mev)**2*&
                                                get_fermi_integral(5,0.0d0)/&
                                                get_fermi_integral(3,0.0d0)


  !!  with heavy

         kappa_opa_scat_nuclei_j(ix^D,1,2,2) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*&
                                                (Ca - 1.0d0 + (mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*&
                                                (2.d0 - Ca - Cv))**2*sigma_0*n_N(ix^D,4)*&
                                                (lprim(ix^D, temp_)/me_mev)**2*get_fermi_integral(5,local_eta_nue(ix^D))/&
                                                get_fermi_integral(3,local_eta_nue(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,2,2,2) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*&
                                                (Ca - 1.0d0 + (mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*&
                                                (2.d0 - Ca - Cv))**2*sigma_0*n_N(ix^D,4)*&
                                                (lprim(ix^D, temp_)/me_mev)**2*get_fermi_integral(5,local_eta_nua(ix^D))/&
                                                get_fermi_integral(3,local_eta_nua(ix^D))

         kappa_opa_scat_nuclei_j(ix^D,3,2,2) = (1.d0/6.d0)*(mass_fraction(ix^D,5)**2)*&
                                                (Ca - 1.0d0 + (mass_fraction(ix^D,6)/mass_fraction(ix^D,5))*&
                                                (2.d0 - Ca - Cv))**2*sigma_0*n_N(ix^D,4)*&
                                                (lprim(ix^D, temp_)/me_mev)**2*get_fermi_integral(5,0.0d0)/&
                                                get_fermi_integral(3,0.0d0)


  {enddo^D&\} 


!    write(*,*) kappa_opa_scat_nucleons_j(5,:,2,2),  " kappa_opa_scat_ neutrobns"
!    write(*,*) kappa_opa_scat_nucleons_j(5,:,1,2),  " kappa_opa_scat_   protons"
!    write(*,*) kappa_opa_scat_nuclei_j(5,:,1,2)  , "kappa_opa_scatting alpha"
!    write(*,*) kappa_opa_scat_nuclei_j(5,:,2,2)  , "kappa_opa_ scat  heavy"
!stop


!!!!!!!!!!! do total bar_opa for 3neu   !here is improved_leakage defined kappa_total
  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
 !!!! j = 0          for 3 neu
       kappa_opa_total_j(ix^D,1,1) = kappa_opa_abs_j(ix^D,1,1) + kappa_opa_scat_nucleons_j(ix^D,1,1,1) +&
                                      kappa_opa_scat_nucleons_j(ix^D,1,2,1) + kappa_opa_scat_nuclei_j(ix^D,1,1,1) +&
                                      kappa_opa_scat_nuclei_j(ix^D,1,2,1)

       kappa_opa_total_j(ix^D,2,1) = kappa_opa_abs_j(ix^D,2,1) + kappa_opa_scat_nucleons_j(ix^D,2,1,1) +&
                                      kappa_opa_scat_nucleons_j(ix^D,2,2,1) + kappa_opa_scat_nuclei_j(ix^D,2,1,1) +&
                                      kappa_opa_scat_nuclei_j(ix^D,2,2,1)

       kappa_opa_total_j(ix^D,3,1) = kappa_opa_scat_nucleons_j(ix^D,3,1,1) + kappa_opa_scat_nucleons_j(ix^D,3,2,1) +&
                                      kappa_opa_scat_nuclei_j(ix^D,3,1,1) + kappa_opa_scat_nuclei_j(ix^D,3,2,1)

 !!!! j = 1     for 3 neu

       kappa_opa_total_j(ix^D,1,2) = kappa_opa_abs_j(ix^D,1,2) + kappa_opa_scat_nucleons_j(ix^D,1,1,2) +&
                                      kappa_opa_scat_nucleons_j(ix^D,1,2,2) + kappa_opa_scat_nuclei_j(ix^D,1,1,2) +&
                                      kappa_opa_scat_nuclei_j(ix^D,1,2,2)

       kappa_opa_total_j(ix^D,2,2) = kappa_opa_abs_j(ix^D,2,2) + kappa_opa_scat_nucleons_j(ix^D,2,1,2) +&
                                      kappa_opa_scat_nucleons_j(ix^D,2,2,2) + kappa_opa_scat_nuclei_j(ix^D,2,1,2) +&
                                      kappa_opa_scat_nuclei_j(ix^D,2,2,2)

       kappa_opa_total_j(ix^D,3,2) = kappa_opa_scat_nucleons_j(ix^D,3,1,2) + kappa_opa_scat_nucleons_j(ix^D,3,2,2) +&
                                      kappa_opa_scat_nuclei_j(ix^D,3,1,2) + kappa_opa_scat_nuclei_j(ix^D,3,2,2)
  {enddo^D&\} 
!write(*,*)  kappa_opa_total(5,1,1,1,1), "kappa_opa_total"
!!!!!!!!!!!!! do total bar_opa for 3neu    !here is simple leakage kappa_total
!   {do ix^D = ixO^LIM^D \}
! !!!! j = 0          for 3 neu
!       kappa_opa_total_j(ix^D,1,1) = kappa_opa_abs_j(ix^D,1,1) + kappa_opa_scat_nucleons_j(ix^D,1,1,1) +&
!                                      kappa_opa_scat_nucleons_j(ix^D,1,2,1) 
!                                   !+ kappa_opa_scat_nuclei_j(ix^D,1,1,1) +&
!                                   !   kappa_opa_scat_nuclei_j(ix^D,1,2,1)
!
!       kappa_opa_total_j(ix^D,2,1) = kappa_opa_abs_j(ix^D,2,1) + kappa_opa_scat_nucleons_j(ix^D,2,1,1) +&
!                                      kappa_opa_scat_nucleons_j(ix^D,2,2,1)
!                                    ! + kappa_opa_scat_nuclei_j(ix^D,2,1,1) +&
!                                    !  kappa_opa_scat_nuclei_j(ix^D,2,2,1)
!       kappa_opa_total_j(ix^D,3,1) = kappa_opa_scat_nucleons_j(ix^D,3,1,1) + kappa_opa_scat_nucleons_j(ix^D,3,2,1) 
!                                        !+&
!                                     ! kappa_opa_scat_nuclei_j(ix^D,3,1,1) + kappa_opa_scat_nuclei_j(ix^D,3,2,1)
!
! !!!! j = 1     for 3 neu
!
!       kappa_opa_total_j(ix^D,1,2) = kappa_opa_abs_j(ix^D,1,2) + kappa_opa_scat_nucleons_j(ix^D,1,1,2) +&
!                                      kappa_opa_scat_nucleons_j(ix^D,1,2,2)
!                                        ! + kappa_opa_scat_nuclei_j(ix^D,1,1,2) +&
!                                    !  kappa_opa_scat_nuclei_j(ix^D,1,2,2)
!
!       kappa_opa_total_j(ix^D,2,2) = kappa_opa_abs_j(ix^D,2,2) + kappa_opa_scat_nucleons_j(ix^D,2,1,2) +&
!                                      kappa_opa_scat_nucleons_j(ix^D,2,2,2)
!                                        ! + kappa_opa_scat_nuclei_j(ix^D,2,1,2) +&
!                                    !  kappa_opa_scat_nuclei_j(ix^D,2,2,2)
!
!       kappa_opa_total_j(ix^D,3,2) = kappa_opa_scat_nucleons_j(ix^D,3,1,2) + kappa_opa_scat_nucleons_j(ix^D,3,2,2)   
!                                                !+&
!                                     ! kappa_opa_scat_nuclei_j(ix^D,3,1,2) + kappa_opa_scat_nuclei_j(ix^D,3,2,2)
!  {enddo^D&\} 
!       write(*,*)  kappa_opa_total_j(5,1,2), "kappa_opa_total_j energy"


!
!       write(*,*)  "kappa_opa_totala_j"
!       write(*,*)  kappa_opa_total_j(5,1,1,:,1)
!       write(*,*)  kappa_opa_total_j(5,1,1,:,2)
!          filename = trim(adjustl(outdir))//"/lepton_blocking.xg"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,*) '"Time = ',time-t_bounce, "ishock =", ishock(1)
!          do i=ghosts1+1,600
!             write(667,"(i8,1P10E15.6)") i,lep_blocking(i,1,1,1)
!          enddo
!          write(667,*) " "
!          write(667,*) " "
!          close(667)
!
!
!          filename = trim(adjustl(outdir))//"/deno1and2.xg"
!          open(667,file=filename,status='unknown',position='append')
!          write(667,*) '"Time = ',time-t_bounce, "ishock =", ishock(1)
!          do i=ghosts1+1,600
!             write(667,"(i8,1P10E15.6)") i,deno(i,1,1,1), deno2(i,1,1,1)
!          enddo
!          write(667,*) " "
!          write(667,*) " "
!          close(667)




end subroutine find_kappa_opa_total



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_zeta_n_N(ixI^L, ixO^L)

   implicit none
  integer, intent(in)             :: ixI^L, ixO^L
  integer jm1,ix^D, ix_next^D
  double precision, dimension(ixI^S,2) :: E_fermi  !Fermi energy for p,n      1=p, 2=n
  double precision, dimension(ixI^S) :: leta_n
  double precision, dimension(ixI^S) :: leta_p
  real*8  get_fermi_integral_half
  real*8  get_fermi_integral
  real*8  xlye, t1, xynp, xypn, xyp, xyn


   n_N(:,1) = 0.0d0
   n_N(:,2) = 0.0d0
   eta_pn = 0.0d0
   eta_np = 0.0d0
   zeta_N = 0.0d0
   E_fermi = 0.0d0
   ksi_NN = 0.0d0
   leta_p = 0.0d0
   leta_n = 0.0d0

  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
   leta_p(ix^D) =  eta_nucleons(ix^D,2) ! - 938/ltemp(:,:,:)         xmu_p  range from  -5 to -30   onyl!!
   leta_n(ix^D) =  eta_nucleons(ix^D,3) ! - 938/ltemp(:,:,:)
  {enddo^D&\} 


    !already  checked ans,   ~ 1d-5 error
!   write(*,*)     get_fermi_integral_half(-5.0d0)/(sqrt(pi)/2.0d0), "-5"
!   write(*,*)     get_fermi_integral_half(0.0d0)/(sqrt(pi)/2.0d0), "0.0"
!   write(*,*)     get_fermi_integral_half(1.0d0)/(sqrt(pi)/2.0d0), "1"
!   write(*,*)     get_fermi_integral_half(3.0d0)/(sqrt(pi)/2.0d0), "3"
!   write(*,*)     get_fermi_integral_half(4.0d0)/(sqrt(pi)/2.0d0), "4"
!stop


!  ILEAS version
  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}


! ! code test   simple leakage to calculate leak_tau
        t1 = dexp(-eta_hat(ix^D))
       xlye = lprim(ix^D, ye_)

!   for code tets    simple leakage to calculate kappa_abs_j
       xynp = max((2.0d0*xlye-1.0d0)/ (t1-1.0d0),0.0d0)
       eta_np(ix^D) = avo*lprim(ix^D, rho_)*xynp
       xypn = max(xynp * t1,0.0d0)
       eta_pn(ix^D) = avo*lprim(ix^D, rho_) *xypn



!!!!!   eta_pn = ksi_pn    !move to compute_zeta_N_n     !simple leakage version when doing diffusion (not completely disscoication)
!       eta_pn(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,1)-mass_fraction(ix^D,2))/(dexp(eta_hat(ix^D))-1.0d0)
!       eta_pn(ix^D) = max(eta_pn(ix^D),0.0d0)
!       eta_np(ix^D) = avo*lprim(ix^D, rho_)*(mass_fraction(ix^D,2)-mass_fraction(ix^D,1))/(dexp(-eta_hat(ix^D))-1.0d0)
!       eta_np(ix^D) = max(eta_np(ix^D),0.0d0)
!       !if (eta_pn(ix^D) .le. 0.0d0  .or.  eta_np(ix^D) .le. 0.0d0) then
!!           stop "eta_pn, np   -ve"
!!       endif
!
!       if (lprim(ix^D, rho_).lt.1.0d11) then
!          !non degenerate here, use mass fractions as chemical potentials fail at low densities
!          eta_pn(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,2)
!          eta_np(ix^D) = avo*lprim(ix^D, rho_)*mass_fraction(ix^D,1)
!       endif

!!!!!!!!!!!!ILEAS
     !  eta_p , eta_n  here  include rest mass of not!!!!
       !Now just follow gr1d    only  eta_pn(ILEAS)  ksi_NN(ILEAS) need these !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      n_N(ix^D,1)  = (4.d0*pi/(hc_mevcm)**3)*(2.d0*mp_mev*lprim(ix^D, temp_))**(3.d0/2.d0)*&
!                      get_fermi_integral_half(leta_p(ix^D))
!
!
!      n_N(ix^D,2)  = (4.d0*pi/(hc_mevcm)**3)*(2.d0*mn_mev*lprim(ix^D, temp_))**(3.d0/2.d0)*&
!                      get_fermi_integral_half(leta_n(ix^D))

!        write(*,*) n_N(ix^D,2)
             !here   we do the New  nucleon blocking factors for
             !degenerate/non-degenerate Fermi gas ,    inverting the relation by
             !using the num_den of proton and neutron

!  ILEAS eta_pn and eta_np    !Rampp. Radiation hydrodynamics with neutrinos: Stellar core collapse and dexplosion mechanism of type
!  II supernovae. PhD thesis, Technische Universit?at M?unchen, 2000. (A40)
!      eta_pn(ix^D) = ( n_N(ix^D,2) - n_N(ix^D,1) )/ &
!                      (dexp(leta_n(ix^D) - leta_p(ix^D) - Qnp/lprim(ix^D, temp_)) - 1.0d0)
!
!      eta_np(ix^D) = ( n_N(ix^D,1) - n_N(ix^D,2) )/ &
!                      (dexp(leta_p(ix^D) - leta_n(ix^D) + Qnp/lprim(ix^D, temp_)) - 1.0d0)
!
!!           if (eta_np(ix^D) .le. 0.0d0 .or. eta_pn(ix^D) .le. 0.0d0) then       !would be -ve in optical thin region
!!            write(*,*)  eta_np(ix^D), eta_pn(ix^D), ix1
!!            write(*,*)  "eta  is -ve "
!!           endif
!
!      eta_np(ix^D) = max(eta_np(ix^D), 0.0d0)    !would be -ve in optical thin region
!      eta_pn(ix^D) = max(eta_pn(ix^D), 0.0d0)
!        write(*,*) eta_np(ix^D), lx(ix^D,r_),"eta_np"
!     unit checked


!      E_fermi(ix^D,1) = (planck_mevs**2/8.0d0/pi**2/(m_b_g))*(3.d0*pi**2*n_N(ix^D,1))**(2.d0/3.d0)
!      E_fermi(ix^D,2) = (planck_mevs**2/8.0d0/pi**2/(m_b_g))*(3.d0*pi**2*n_N(ix^D,2))**(2.d0/3.d0)   !m_b_g in g
!      E_fermi(ix^D,1) = (planck_mevs**2/8.0d0/pi**2/(m_b_mev))*(3.d0*pi**2*n_N(ix^D,1))**(2.d0/3.d0)
!      E_fermi(ix^D,2) = (planck_mevs**2/8.0d0/pi**2/(m_b_mev))*(3.d0*pi**2*n_N(ix^D,2))**(2.d0/3.d0)   ! m_b_mev  of coz needs MeV
!
!      zeta_N(ix^D,1) = 3.d0*lprim(ix^D, temp_)/2.d0/E_fermi(ix^D,1)
!      zeta_N(ix^D,2) = 3.d0*lprim(ix^D, temp_)/2.d0/E_fermi(ix^D,2)


  {enddo^D&\} 

!pass  !pass even for ILEAS (small difference)
!write(*,*) eta_np(5), eta_pn(5), "eta_np"
!write(*,*) E_fermi(5,:), "E_fermi"
!write(*,*) zeta_N(5,:), "zeta_N"
!stop


!      write(*,*)  "compute zeta"
!      write(*,*)   zeta_N(5,1,1,:), "zeta"
!      write(*,*)   n_N(:,1,1,:)  , "n_N"
!      write(*,*)   E_fermi(:,1,1,:),   "E_fermi"


  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
!!!! simple leakage leawk_tau

!  increase so much the leak_tau
          xyn = (1.0d0-lprim(ix^D, ye_)) / (1.0d0 + 2.0d0/3.0d0 * max(eta_nucleons(ix^D,3),0.0d0))
          xyp = lprim(ix^D, ye_) / (1.0d0 + 2.0d0/3.0d0*max(eta_nucleons(ix^D,2),0.0d0))

        ksi_NN(ix^D,1) =  avo * lprim(ix^D, rho_) * xyp
        ksi_NN(ix^D,2) =  avo * lprim(ix^D, rho_) * xyn

!   ILEAS
!! ksi _p
!       ksi_NN(ix^D,1) = (zeta_N(ix^D,1)/sqrt(1.d0 + zeta_N(ix^D,1)**2)) * avo * lprim(ix^D, rho_) *&
!                         mass_fraction(ix^D, 2) 
!
!! ksi_n
!       ksi_NN(ix^D,2) = (zeta_N(ix^D,2)/sqrt(1.d0 + zeta_N(ix^D,2)**2)) * avo * lprim(ix^D, rho_) *&
!                         mass_fraction(ix^D, 1) 

        !write(*,*) ksi_NN(ix^D,1) , "ksi_NN"
  {enddo^D&\} 


end subroutine compute_zeta_n_N


subroutine time_prod_and_bar_E(ixI^L, ixO^L)

  implicit none
    integer, intent(in)             :: ixI^L, ixO^L
   real*8 get_fermi_integral
   integer :: i,j,k, ix^D
   time_prod = 0.0d0
   bar_E_j = 0.0d0

!  get integrated over neutrino spectrum neutrino energy and number density
  {do ix^D = ixI^LIM^D \}
  !{do ix^D = ixO^LIM^D \}
!!!!!! j = 0
       bar_E_j(ix^D,1,1) = (4.d0*pi/(hc_mevcm**3))*lprim(ix^D, temp_)**3*get_fermi_integral(2,eta_nu(ix^D,1))
       bar_E_j(ix^D,2,1) = (4.d0*pi/(hc_mevcm**3))*lprim(ix^D, temp_)**3*get_fermi_integral(2,eta_nu(ix^D,2))
       bar_E_j(ix^D,3,1) = 4.0d0*(4.d0*pi/(hc_mevcm**3))*lprim(ix^D, temp_)**3*get_fermi_integral(2,eta_nu(ix^D,3))
!!!!!! j = 1
       bar_E_j(ix^D,1,2) = (4.d0*pi/(hc_mevcm**3))*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nu(ix^D,1))
       bar_E_j(ix^D,2,2) = (4.d0*pi/(hc_mevcm**3))*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nu(ix^D,2))
       bar_E_j(ix^D,3,2) = 4.0d0*(4.d0*pi/(hc_mevcm**3))*lprim(ix^D, temp_)**4*get_fermi_integral(3,eta_nu(ix^D,3))
  {enddo^D&\}


!!!!!!#######################!production time scale
    {do ix^D = ixI^LIM^D \}
    !{do ix^D = ixO^LIM^D \}
         time_prod(ix^D,1,1) = bar_E_j(ix^D,1,1)/R_loc(ix^D,1)
         time_prod(ix^D,2,1) = bar_E_j(ix^D,2,1)/R_loc(ix^D,2)
         time_prod(ix^D,3,1) = bar_E_j(ix^D,3,1)/R_loc(ix^D,3)


         time_prod(ix^D,1,2) = bar_E_j(ix^D,1,2)/Q_loc(ix^D,1)
         time_prod(ix^D,2,2) = bar_E_j(ix^D,2,2)/Q_loc(ix^D,2)
         time_prod(ix^D,3,2) = bar_E_j(ix^D,3,2)/Q_loc(ix^D,3)
     {enddo^D&\}

end subroutine 





  subroutine finding_intermediate_variables(ixI^L, ixO^L)
   ! use mod_global_parameters
    implicit none

    integer, intent(in) :: ixI^L, ixO^L
    double precision    :: W2v2(ixI^S)
    double precision    :: lfac2(ixI^S)    
    integer :: idir, ix^D

     ! calculate lfac  (spherical coor)
         {do ix^D = ixI^LIM^D \}
         !{do ix^D = ixO^LIM^D \}
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


! output  dEdx  array
 subroutine gradient_any_coordinate_five_pt(ixI^L, ixO^L, array, dEdx, input_d)
    use mod_geometry
   implicit none

    double precision, intent(in), dimension(ixI^S) :: array
    double precision, intent(inout), dimension(ixI^S) :: dEdx       !  input:  dEdx(:, r_) or dEdx(:, theta_)
    double precision :: five_pt_stencils
    integer, intent(in)                            :: ixI^L, ixO^L
     integer :: ix^D, ix_2plus^D, ix_1plus^D, ix_2minus^D, ix_1minus^D, input_d
 
! cartesian coordinate
  ! Del E = dEdx x^ + dEdy y^ + dEdz z^

  if (coordinate == cartesian) then 
    !  5 pt stencil needs  +2 , -2 grids
    {do ix^D = ixImin^D, ixImax^D-1 \}
    !  for later second derivative  need 2 more ghostzone first derivative
!    {do ix^D = ixOmin^D-2, ixOmax^D \}
    !{do ix^D = ixO^LIM^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_2plus^D = ix^D + kr(input_d,^D) + kr(input_d,^D)
    ix_1minus^D = ix^D - kr(input_d,^D)   
    ix_2minus^D = ix^D - kr(input_d,^D) - kr(input_d,^D)  
 
        dEdx(ix^D) =  five_pt_stencils(array(ix_2plus^D), array(ix_1plus^D), array(ix_1minus^D), array(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
   
    {enddo^D&\} ! end loop for every points

  elseif (coordinate == spherical) then
  ! Del E = dEdr r^ + 1/r dEdtheta theta^ + 1/r/sin(theta) dEdphi phi^
    !  5 pt stencil needs  +2 , -2 grids

!    {do ix^D = ixImin^D, ixImax^D-1 \}
!    {do ix^D = ixI^LIM^D \}
  !  {do ix^D = ixOmin^D-1, ixOmax^D \}
    !  for later second derivative  need 2 more ghostzone first derivative
    {do ix^D = ixOmin^D-2, ixOmax^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_2plus^D = ix^D + kr(input_d,^D) + kr(input_d,^D)
    ix_1minus^D = ix^D - kr(input_d,^D)   
    ix_2minus^D = ix^D - kr(input_d,^D) - kr(input_d,^D)  

        if (input_d == 1) then 
        dEdx(ix^D) =  five_pt_stencils(array(ix_2plus^D), array(ix_1plus^D), array(ix_1minus^D), array(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))

        else if (input_d == 2) then
{^NOONED
        dEdx(ix^D) =  five_pt_stencils(array(ix_2plus^D), array(ix_1plus^D), array(ix_1minus^D), array(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
        dEdx(ix^D) = dEdx(ix^D) * 1.0d0 / lx(ix^D, r_)
}
        else if (input_d == 3) then
{^NOONED
        dEdx(ix^D) =  five_pt_stencils(array(ix_2plus^D), array(ix_1plus^D), array(ix_1minus^D), array(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
        dEdx(ix^D) = dEdx(ix^D) * 1.0d0/ lx(ix^D, r_) / dsin(lx(ix^D, theta)) 
}
        endif
    {enddo^D&\} ! end loop for every points
  endif
! spherical coordinate
 end subroutine


 subroutine div_any_coordinate_five_pt(ixI^L, ixO^L, array, div_dEdx, input_d)
    use mod_geometry
   implicit none

    double precision, intent(in), dimension(ixI^S) :: array
    double precision, dimension(ixI^S) :: array_local
    double precision, intent(inout), dimension(ixI^S) :: div_dEdx       !  input:  div_dEdx(:, r_) or div_dEdx(:, theta_)
    double precision :: five_pt_stencils
    integer, intent(in)                            :: ixI^L, ixO^L
    integer :: ix^D, ix_2plus^D, ix_1plus^D, ix_2minus^D, ix_1minus^D, input_d
! cartesian coordinate
  !  Div_E = (dE_x/dx) + (dE_y/dy) + (dE_z/dz)
!     five_pt_stencils = 0.0d0
     array_local = array

  if (coordinate == cartesian) then 
    !  5 pt stencil needs  +2 , -2 grids

    {do ix^D = ixO^LIM^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_2plus^D = ix^D + kr(input_d,^D) + kr(input_d,^D)
    ix_1minus^D = ix^D - kr(input_d,^D)   
    ix_2minus^D = ix^D - kr(input_d,^D) - kr(input_d,^D)  
 
        div_dEdx(ix^D) =  five_pt_stencils(array_local(ix_2plus^D), array_local(ix_1plus^D), array_local(ix_1minus^D), array_local(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
   
    {enddo^D&\} ! end loop for every points

  elseif (coordinate == spherical) then
  ! Div_E = (d(E_r r**2)/dr) / r**2 + (d(E_theta sin(theta))/dtheta) / (r sin(theta)) + (dE_phi/dphi) / (r sin(theta))
    !  5 pt stencil needs  +2 , -2 grids
        if (input_d == 1) then    ! as  d(r**2 E)/dr  =   d r(ix+1)**2 E(ix+1) - d r(ix)**2 E(ix)
           {do ix^D = ixO^LIM^D \}
                array_local(ix^D) = array_local(ix^D) * lx(ix^D, r_)**2    !  (E_r * r**2)
           {enddo^D&\} ! end loop for every points
        else if (input_d == 2) then
           {do ix^D = ixO^LIM^D \}
                array_local(ix^D) = array_local(ix^D) * dsin(lx(ix^D, theta_))  ! (E_theta * sin(theta))
           {enddo^D&\} ! end loop for every points
        else if (input_d == 3) then
                !nothing
        endif

    {do ix^D = ixO^LIM^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_2plus^D = ix^D + kr(input_d,^D) + kr(input_d,^D)
    ix_1minus^D = ix^D - kr(input_d,^D)   
    ix_2minus^D = ix^D - kr(input_d,^D) - kr(input_d,^D)  

        if (input_d == 1) then
        div_dEdx(ix^D) =  five_pt_stencils(array_local(ix_2plus^D), array_local(ix_1plus^D), array_local(ix_1minus^D), array_local(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))

        else if (input_d == 2) then
{^NOONED
        div_dEdx(ix^D) =  five_pt_stencils(array_local(ix_2plus^D), array_local(ix_1plus^D), array_local(ix_1minus^D), array_local(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
}
        else if (input_d == 3) then
{^NOONED
        div_dEdx(ix^D) =  five_pt_stencils(array_local(ix_2plus^D), array_local(ix_1plus^D), array_local(ix_1minus^D), array_local(ix_2minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
}
        endif
    {enddo^D&\} ! end loop for every points


        if (input_d == 1) then    ! as  d(r**2 E)/dr  =   d r(ix+1)**2 E(ix+1) - d r(ix)**2 E(ix)
           {do ix^D = ixO^LIM^D \}
                div_dEdx(ix^D) = div_dEdx(ix^D) / (lx(ix^D, r_)**2)
           {enddo^D&\} ! end loop for every points
        else if (input_d == 2) then
           {do ix^D = ixO^LIM^D \}
                div_dEdx(ix^D) = div_dEdx(ix^D) * 1.0d0 / lx(ix^D, r_) / dsin(lx(ix^D, theta_))
           {enddo^D&\} ! end loop for every points
        else if (input_d == 3) then
           {do ix^D = ixO^LIM^D \}
                div_dEdx(ix^D) = div_dEdx(ix^D) * 1.0d0/ lx(ix^D, r_) / dsin(lx(ix^D, theta_)) 
           {enddo^D&\} ! end loop for every points
        endif
   endif
 end subroutine


!################################## 3 pt stencils
! output  dEdx  array
 subroutine gradient_any_coordinate_two_pt(ixI^L, ixO^L, array, dEdx, input_d)
    use mod_geometry
   implicit none

    double precision, intent(in), dimension(ixI^S) :: array
    double precision, intent(inout), dimension(ixI^S) :: dEdx       !  input:  dEdx(:, r_) or dEdx(:, theta_)
    double precision :: two_pt_stencils
    integer, intent(in)                            :: ixI^L, ixO^L
     integer :: ix^D, ix_1plus^D, ix_1minus^D, input_d

!    two_pt_stencils = 0.0d0 
! cartesian coordinate
  ! Del E = dEdx x^ + dEdy y^ + dEdz z^

  if (coordinate == cartesian) then 
    !  5 pt stencil needs  +2 , -2 grids
    {do ix^D = ixImin^D, ixImax^D-1 \}
    !  for later second derivative  need 2 more ghostzone first derivative
!    {do ix^D = ixOmin^D-2, ixOmax^D \}
    !{do ix^D = ixO^LIM^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_1minus^D = ix^D - kr(input_d,^D) 
 
        dEdx(ix^D) =  two_pt_stencils(array(ix_1plus^D), array(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
   
    {enddo^D&\} ! end loop for every points

  elseif (coordinate == spherical) then
  ! Del E = dEdr r^ + 1/r dEdtheta theta^ + 1/r/sin(theta) dEdphi phi^
    !  5 pt stencil needs  +2 , -2 grids

!    {do ix^D = ixImin^D, ixImax^D-1 \}
!    {do ix^D = ixI^LIM^D \}
  !  {do ix^D = ixOmin^D-1, ixOmax^D \}
    !  for later second derivative  need 2 more ghostzone first derivative
    {do ix^D = ixOmin^D-2, ixOmax^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_1minus^D = ix^D - kr(input_d,^D) 

        if (input_d == 1) then 
        dEdx(ix^D) =  two_pt_stencils(array(ix_1plus^D), array(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))

        else if (input_d == 2) then
{^NOONED
        dEdx(ix^D) =  two_pt_stencils(array(ix_1plus^D), array(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
        dEdx(ix^D) = dEdx(ix^D) * 1.0d0 / lx(ix^D, r_)
}
        else if (input_d == 3) then
{^NOONED
        dEdx(ix^D) =  two_pt_stencils(array(ix_1plus^D), array(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
        dEdx(ix^D) = dEdx(ix^D) * 1.0d0/ lx(ix^D, r_) / dsin(lx(ix^D, theta)) 
}
        endif
    {enddo^D&\} ! end loop for every points
  endif
! spherical coordinate
 end subroutine


 subroutine div_any_coordinate_two_pt(ixI^L, ixO^L, array, div_dEdx, input_d)
    use mod_geometry
   implicit none

    double precision, intent(in), dimension(ixI^S) :: array
    double precision, dimension(ixI^S) :: array_local
    double precision, intent(inout), dimension(ixI^S) :: div_dEdx       !  input:  div_dEdx(:, r_) or div_dEdx(:, theta_)
    double precision :: two_pt_stencils
    integer, intent(in)                            :: ixI^L, ixO^L
    integer :: ix^D, ix_1plus^D, ix_1minus^D, input_d
! cartesian coordinate
  !  Div_E = (dE_x/dx) + (dE_y/dy) + (dE_z/dz)
     array_local = array
!     two_pt_stencils = 0.0d0

  if (coordinate == cartesian) then 
    !  5 pt stencil needs  +2 , -2 grids

    {do ix^D = ixO^LIM^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_1minus^D = ix^D - kr(input_d,^D) 
 
        div_dEdx(ix^D) =  two_pt_stencils(array_local(ix_1plus^D), array_local(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
   
    {enddo^D&\} ! end loop for every points

  elseif (coordinate == spherical) then
  ! Div_E = (d(E_r r**2)/dr) / r**2 + (d(E_theta sin(theta))/dtheta) / (r sin(theta)) + (dE_phi/dphi) / (r sin(theta))
    !  5 pt stencil needs  +2 , -2 grids
        if (input_d == 1) then    ! as  d(r**2 E)/dr  =   d r(ix+1)**2 E(ix+1) - d r(ix)**2 E(ix)
           {do ix^D = ixO^LIM^D \}
                array_local(ix^D) = array_local(ix^D) * lx(ix^D, r_)**2    !  (E_r * r**2)
           {enddo^D&\} ! end loop for every points
        else if (input_d == 2) then
           {do ix^D = ixO^LIM^D \}
                array_local(ix^D) = array_local(ix^D) * dsin(lx(ix^D, theta_))  ! (E_theta * sin(theta))
           {enddo^D&\} ! end loop for every points
        else if (input_d == 3) then
                !nothing
        endif

    {do ix^D = ixO^LIM^D \}
    ix_1plus^D = ix^D + kr(input_d,^D)   
    ix_1minus^D = ix^D - kr(input_d,^D) 

        if (input_d == 1) then
        div_dEdx(ix^D) =  two_pt_stencils(array_local(ix_1plus^D), array_local(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))

        else if (input_d == 2) then
{^NOONED
        div_dEdx(ix^D) =  two_pt_stencils(array_local(ix_1plus^D), array_local(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
}
        else if (input_d == 3) then
{^NOONED
        div_dEdx(ix^D) =  two_pt_stencils(array_local(ix_1plus^D), array_local(ix_1minus^D),&
                lxi(ix_1plus^D,input_d), lxi(ix^D,input_d))
}
        endif
    {enddo^D&\} ! end loop for every points


        if (input_d == 1) then    ! as  d(r**2 E)/dr  =   d r(ix+1)**2 E(ix+1) - d r(ix)**2 E(ix)
           {do ix^D = ixO^LIM^D \}
                div_dEdx(ix^D) = div_dEdx(ix^D) / (lx(ix^D, r_)**2)
           {enddo^D&\} ! end loop for every points
        else if (input_d == 2) then
           {do ix^D = ixO^LIM^D \}
                div_dEdx(ix^D) = div_dEdx(ix^D) * 1.0d0 / lx(ix^D, r_) / dsin(lx(ix^D, theta_))
           {enddo^D&\} ! end loop for every points
        else if (input_d == 3) then
           {do ix^D = ixO^LIM^D \}
                div_dEdx(ix^D) = div_dEdx(ix^D) * 1.0d0/ lx(ix^D, r_) / dsin(lx(ix^D, theta_)) 
           {enddo^D&\} ! end loop for every points
        endif
   endif
 end subroutine
end module mod_grhd_ccsn_ILEAS_ILEAS





function five_pt_stencils(var_plus2, var_plus1, var_minus1, var_minus2, lxi_plus, lxi_minus)
   implicit none
   
   ! input grid_number
   !  y(i+2), y(i+1), y(i-1), y(i-2)
   double precision, intent(in)    :: var_plus2, var_plus1, var_minus1, var_minus2, &
                                      lxi_plus, lxi_minus
   !  inner cell coordinate  corresponding to dimension - D
   double precision :: dx, five_pt_stencils
     
   five_pt_stencils = 0.0d0
   dx = 0.0d0     

   dx = lxi_plus - lxi_minus 
   
   five_pt_stencils = ( - var_plus2 + 8.0d0*var_plus1 - 8.0d0*var_minus1 + var_minus2)/ 12.0d0/dx   
       
!  write(*,*)  var_minus2, var_minus1, var_plus1, var_plus2, "in five pt"
!        write(*,*) lxi_plus, lxi_minus, dx , five_pt_stencils

 return 
end function 

!   var_plus1 = next grid,  var_minus1 = this grid   fixme
function two_pt_stencils(var_plus1, var_minus1, lxi_plus, lxi_minus)
   implicit none
   
   ! input grid_number
   !  y(i+2), y(i+1), y(i-1), y(i-2)
   double precision, intent(in)    :: var_plus1, var_minus1, &
                                      lxi_plus, lxi_minus
   !  inner cell coordinate  corresponding to dimension - D
   double precision :: dx, two_pt_stencils
   two_pt_stencils = 0.0d0
   dx = 0.0d0     

   dx = lxi_plus - lxi_minus 
   
   two_pt_stencils = (var_plus1 - var_minus1) / 2.0d0/dx
       

 return 
end function 

!  we only consider j = 1/2   for calculate the number density for p/n P. Van Halen, and D. L. Pulfrey 1984
function get_fermi_integral_half(eta)
   implicit none
  real*8 get_fermi_integral_half
  real*8 j
  real*8 F_x
  real*8 eta, r
  real*8 z(7), a(7) , b(7), c(7)
  real*8  gamma_half
  integer :: iter
! for x < 0   !table 1      for eq(4)
  z(1) = 1.0d0
  z(2) = 0.353568d0
  z(3) = 0.192439d0
  z(4) = 0.122973d0
  z(5) = 0.077134d0
  z(6) = 0.036228d0
  z(7) = 0.008346d0

! for    x   >  4   , coefficients   eq(6)   table 2
  a(1) = 0.752253d0
  a(2) = 0.928195d0
  a(3) = 0.680839d0
  a(4) = 25.7829d0
  a(5) = -553.636d0
  a(6) = 3531.43d0
  a(7) = -3254.65d0
! for  0 < x < 2        table 3
  b(1) = 0.765147d0
  b(2) = 0.604911d0
  b(3) = 0.189885d0
  b(4) = 0.020307d0
  b(5) = -0.004380d0
  b(6) = -0.000366d0
  b(7) = 0.000133d0
! for 2  < x  < 4       table 3
  c(1) = 0.777114d0
  c(2) = 0.581307d0
  c(3) = 0.206132d0
  c(4) = 0.017680d0
  c(5) = -0.006549d0
  c(6) = 0.000784d0
  c(7) = -0.000036d0

  gamma_half = 0.8862269d0   !sqrt(pi)/2.0d0        !gamma function( 1/2 +1)
  r = 1.0d0
  j = 0.5d0        ! our fermi integral  j = 1/2

  F_x = 0.0d0

   if (eta  .le. 0.0d0) then    !eq(4)

        do iter = 1, 7

           F_x = F_x + (-1.0d0)**(r + 1.0d0) * z(iter) * dexp(r * eta)

           r = r + 1.0d0
        enddo
   else if (eta .ge. 4.0d0  ) then   !eq(6)

        do iter = 1, 7

           F_x = F_x + ( a(iter) /  eta**( 2.0d0 * ( r - 1.0d0 )) )

           r = r + 1.0d0
        enddo

        F_x = F_x * eta**( j + 1.0d0)
   else if (eta .lt. 4.0d0 .and. eta .ge. 2.0d0) then    !eq(7)

        do iter = 1, 7
           F_x = F_x + c(iter) * eta**( r - 1.0d0)

           r = r + 1.0d0
        enddo
   else  !      0.0d0 < eta < 2.0d0   !eq(7)


        do iter = 1, 7
           F_x = F_x + b(iter) * eta**( r - 1.0d0)

           r = r + 1.0d0
        enddo
   endif



   get_fermi_integral_half = F_x * gamma_half
return
end function get_fermi_integral_half



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
             log10(1.0d0+dexp(eta))
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
             log10(1.0d0+dexp(eta))
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


