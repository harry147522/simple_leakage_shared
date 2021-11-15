module mod_grhd_phys_one_d_core
  use mod_physics
  use mod_grhd_phys_parameters, only: r_core, oneDcore

  implicit none
  private

  ! Public methods
  public :: grhd_phys_one_d_core_init

contains

  !> Initialize the module
  subroutine grhd_phys_one_d_core_init()
    use mod_usr_methods

    integer :: itr, idir

    
    if ( oneDcore ) then
      usr_set_flux => grhd_set_flux
      usr_internal_bc => grhd_internal_bc
    end if
   

  end subroutine grhd_phys_one_d_core_init

  
  !> 1D core treatment: internal_bc
  subroutine grhd_internal_bc(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    ! fixme: should be averaged radial direction

    ! set vtheta and vphi = zero
    where ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_) <= r_core )
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,W_vel(2)) = 0.0d0
      ! Note: v3 is large for rapity rotating NS
      !w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,W_vel(3)) = 0.0d0
    end where
  end subroutine grhd_internal_bc

  !> 1D core treatment: flux
  subroutine grhd_set_flux(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
     ixCmax1,ixCmax2,idim,fC,xi)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixCmin1,&
       ixCmin2,ixCmax1,ixCmax2, idim
    ! face-center xi
    double precision,intent(in)  :: xi(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    ! face-center flux
    double precision,intent(inout) :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ncons,1:ndim)

    integer                      :: i_cons
    
    ! nothing to do for radial direction flux
    ! fixme: should be averaged
    if (idim == 1) return

    ! non-radial direction flux are set to be zero
    do i_cons = 1, ncons
       where (xi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,r_) <= r_core )
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,i_cons,idim) = 0.0d0
       end where
    end do
  end subroutine grhd_set_flux
 

end module mod_grhd_phys_one_d_core
