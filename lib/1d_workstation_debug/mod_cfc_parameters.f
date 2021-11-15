module mod_cfc_parameters
  use mod_global_parameters, only: bigdouble, biginteger
  implicit none
  public 

  ! last metric update time
  double precision                          ::  cfc_t_last_update = 0.0d0
  integer                                   ::  cfc_it_last_update = 0

  
  ! tolerence for initializing psi
  double precision                          ::  cfc_psi_tol_init = 1.0d-12

  ! maximun res for 1: alp, 2: psi, 3: beta/X
  double precision, dimension(1:3)          ::  cfc_tol = (/ 1.0d-6, 1.0d-6,&
      1.0d-6 /)
  ! maximun iteration for each metric solver
  integer                                   ::  cfc_it_max = 1000
  ! Print out the status of the solver at every cfc_print cfc iteration
  integer                                   ::  cfc_print = 10000000

  ! solve the metric at every N steps
  integer                                   ::  cfc_dit_update = biginteger
  ! allowed smallest time step for metric solver (for dit_update only)
  double precision                          ::  cfc_smallest_dt = 0.0d0
  ! solve the metric at every dt
  double precision                          ::  cfc_dt_update = bigdouble

  ! N smooths will be applied for 1: cycling up, 2: cycling down
  integer, dimension(1:2)                   :: cfc_n_cycle = (/2,2/)
  logical                                   :: cfc_redblack = .True.

  logical                                   :: use_cfc = .False.
  logical                                   :: update_metric = .False.

  logical                                   :: reconstruct_cfc = .True.

  ! number of points used to interpolate metric at the cell interface
  integer                                   :: cfc_n_interpolation = 4

end module mod_cfc_parameters
