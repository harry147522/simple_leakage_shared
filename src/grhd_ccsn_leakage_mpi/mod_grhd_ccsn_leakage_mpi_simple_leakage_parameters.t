!> Module containing all hydrodynamics
module mod_grhd_ccsn_leakage_mpi_simple_leakage_parameters
  implicit none
  public
      integer :: nleak = 0



      integer, allocatable  :: lcoolingsource(:)   ! ndir are mom, ndir + 1 : tau, ndir +2 : Dye
!      integer  :: lprim   ! cgs primvarialbes
! inside prim      integer  :: lleak_tau   ! read from cons(:,leak_tau(:))
!for mapping
!      integer  :: lprim_local   ! cgs primvarialbes
!      integer  :: lleak_tau_local   ! read from cons(:,leak_tau(:))
!      integer  :: lx    ! cgs coordinate local     ! only  depend 1:ndim to see r,theta,phi
!      integer  :: lxi   ! cgs cell face coordinate local
!      integer  :: lx_map    ! cgs coordinate local  for mapping, as length_gf will chane a little bit
!metric
!      integer  :: lfac    ! Lorentz factor
!      integer  :: gamma_ij
! leakage
      !number emission rate used in determining change in ye
      integer, allocatable  :: R_eff(:) !effective number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
      integer, allocatable  :: R_loc(:) !local number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
      integer, allocatable  :: R_diff(:) !diffusive number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
      integer, allocatable  :: R_tot(:)

      integer, allocatable  :: Q_eff(:) !effective energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
      !energy emission rate used in determining change in energy
      integer, allocatable  :: Q_loc(:) !local energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
      integer, allocatable  :: Q_diff(:) !diffusive energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3
      integer, allocatable  :: Q_tot(:)

!      integer  :: chi !optical depth with E^2 factored out,indices <radial zones:neutrino species> - 1 / MeV^2
!      integer  :: zeta !Eq. A21, mean free path with E^2 factored out, indices <radial zones:neutrino species> - 1 / MeV^2 / cm


      integer  :: lvol !cell volume in cm^3, indices <radial position>
      integer  :: lIarea !inverse cell area in cm^-2, indices <radial position>
      integer  :: lmass !cell mass in g, indices <radial position>

      !degeneracy parameters
         ! 1: electrons  2: p  3: n
      integer, allocatable  :: eta_nucleons(:) !e degeneracy, indices <radial position>, including rest mass - dimension(:^D&)less

      integer  :: eta_hat !n degeneracy - p dengeneracy - Qnp/temp, removes mass difference, indices <radial position> - dimensionless
         ! 1: nue  2: nua  3: nux
      integer, allocatable  :: eta_nu(:) !nu_e degeneracy, from eta_e - eta_n + eta_p, indices <radial position> - dimensionless

      integer, allocatable  :: lepton_blocking(:) !blocking terms for electrons and positrons, indices <radial position,lepton type> - dimensionless

      !mass fractions
      ! 1. xxn   2. xxp   3.  xxa   4. xxh   5. xabar   6. xzbar
      integer, allocatable  :: mass_fraction(:)

      integer  :: eta_pn !proton number density corrected for neutron final state blocking, indices <radial position> - dimensionless
      integer  :: eta_np !neutron number density corrected for proton final state blocking, indices <radial position> - dimensionless



    !  integer, allocatable  :: lum(:) 

  !> Maximum number of leakage variables
  integer, parameter :: max_nvar_leak = 500

  !> leak variable names
  character(len=16) :: leak_names(max_nvar_leak)
   
  contains   

  !> Set leakage variable
  function var_set_leakvar(name_leak, ix) result(iw)
    character(len=*), intent(in)  :: name_leak !< leakage var name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)

    integer                       :: iw

    ! total number of primitive variables
    nleak  = nleak + 1
    iw     = nleak


    if (.not. present(ix)) then
      leak_names(nleak) = name_leak
    else
      write(leak_names(nleak),"(A,I0)") name_leak, ix
    end if


  end function var_set_leakvar


end module mod_grhd_ccsn_leakage_mpi_simple_leakage_parameters
