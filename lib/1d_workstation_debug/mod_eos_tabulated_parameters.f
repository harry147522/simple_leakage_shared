!> Module for eos
module mod_eos_tabulated_parameters

  implicit none
  public

  character(len=128)                      :: eos_table_name = '' 
  double precision                        :: eos_precision = 1.0d-9

  integer                                 :: nrho,ntemp,nye

  integer                                 :: warn_from !warn from given reflevel

  double precision                        :: energy_shift = 0.0d0

  double precision                        :: t_max_hack = 240.0d0

  ! index variable mapping:
  !  1 -> logpress
  !  2 -> logenergy
  !  3 -> entropy
  !  4 -> munu
  !  5 -> cs2
  !  6 -> dedT
  !  7 -> dpdrhoe
  !  8 -> dpderho
  !  9 -> muhat
  ! 10 -> mu_e
  ! 11 -> mu_p
  ! 12 -> mu_n
  ! 13 -> xa
  ! 14 -> xh
  ! 15 -> xn
  ! 16 -> xp
  ! 17 -> abar
  ! 18 -> zbar
  ! 19 -> gamma


  ! basics
  integer                                 :: nvars = 19
  double precision, allocatable, save     :: eos_tables(:,:,:,:)
  integer, parameter                      :: i_logpress = 1
  integer, parameter                      :: i_logenergy = 2
  integer, parameter                      :: i_entropy = 3
  integer, parameter                      :: i_munu = 4
  integer, parameter                      :: i_cs2 = 5 
  integer, parameter                      :: i_dedT = 6
  integer, parameter                      :: i_dpdrhoe = 7
  integer, parameter                      :: i_dpderho = 8
  integer, parameter                      :: i_muhat = 9
  integer, parameter                      :: i_mu_e = 10
  integer, parameter                      :: i_mu_p = 11
  integer, parameter                      :: i_mu_n = 12
  integer, parameter                      :: i_xa = 13
  integer, parameter                      :: i_xh = 14
  integer, parameter                      :: i_xn = 15
  integer, parameter                      :: i_xp = 16
  integer, parameter                      :: i_abar = 17
  integer, parameter                      :: i_zbar = 18
  integer, parameter                      :: i_gamma = 19

  double precision, allocatable :: logrho_table(:)
  double precision, allocatable :: logtemp_table(:)
  double precision, allocatable :: ye_table(:)

  ! integer that stores if the speed of
  ! sound in the table has already been
  ! divided by the specific enthalpy
  integer :: have_rel_cs2
   
end module mod_eos_tabulated_parameters
