!> Module containing all hydrodynamics
module mod_grhd_ccsn_leakage
  use mod_grhd_ccsn_leakage_phys
  use mod_grhd_ccsn_leakage_ppm
  use mod_grhd_ccsn_leakage_simple_leakage

  use mod_gmunu

  implicit none
  public

contains

  subroutine grhd_ccsn_leakage_activate()
    call grhd_ccsn_leakage_phys_init()
    call grhd_ccsn_leakage_ppm_init()
!    call grhd_simple_leakage_init()
  end subroutine grhd_ccsn_leakage_activate

end module mod_grhd_ccsn_leakage
