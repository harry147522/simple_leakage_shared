!> Module containing all hydrodynamics
module mod_grhd_ccsn_ILEAS
  use mod_grhd_ccsn_ILEAS_phys
  use mod_grhd_ccsn_ILEAS_ppm
  use mod_grhd_ccsn_ILEAS_ILEAS

  use mod_gmunu

  implicit none
  public

contains

  subroutine grhd_ccsn_ILEAS_activate()
    call grhd_ccsn_ILEAS_phys_init()
    call grhd_ccsn_ILEAS_ppm_init()
!    call grhd_ILEAS_init()
  end subroutine grhd_ccsn_ILEAS_activate

end module mod_grhd_ccsn_ILEAS
