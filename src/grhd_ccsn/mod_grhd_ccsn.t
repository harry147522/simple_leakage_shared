!> Module containing all hydrodynamics
module mod_grhd_ccsn
  use mod_grhd_ccsn_phys
  use mod_grhd_ccsn_ppm

  use mod_gmunu

  implicit none
  public

contains

  subroutine grhd_ccsn_activate()
    call grhd_ccsn_phys_init()
    call grhd_ccsn_ppm_init()
  end subroutine grhd_ccsn_activate

end module mod_grhd_ccsn
