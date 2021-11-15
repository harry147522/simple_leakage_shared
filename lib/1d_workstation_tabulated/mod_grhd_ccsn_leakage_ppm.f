!> General relativistic hydrodynamics PPM module
module mod_grhd_ccsn_leakage_ppm
  use mod_physics

  implicit none
  private

  ! Public methods
  public :: grhd_ccsn_leakage_ppm_init

contains

  !> Initialize the module
  subroutine grhd_ccsn_leakage_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => grhd_ccsn_leakage_ppm_flatcd
    phys_ppm_flatsh => grhd_ccsn_leakage_ppm_flatsh
  end subroutine grhd_ccsn_leakage_ppm_init

  subroutine grhd_ccsn_leakage_ppm_flatcd(rec_from, rec_to, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, ixLmin1,ixLmax1, ixRmin1,ixRmax1, w, d2w, drho, dp)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: rec_from, rec_to, ixImin1,ixImax1,&
        ixOmin1,ixOmax1, ixLmin1,ixLmax1, ixRmin1,ixRmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nprim),&
       d2w(ixGlo1:ixGhi1,rec_from:rec_to)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1), dp(ixGlo1:ixGhi1)

    double precision                :: gamma_eos_tmp(ixImin1:ixImax1),&
        d2w_press(ixImin1:ixImax1)
    integer                         :: ix1 

    ! cs2 needed to be updated based on the updated prim vars
    call phys_get_csound2(ixImin1,ixImax1, ixOmin1,ixOmax1,w(ixImin1:ixImax1,&
        1:nprim), gamma_eos_tmp(ixOmin1:ixOmax1))
    gamma_eos_tmp(ixOmin1:ixOmax1) = gamma_eos_tmp(ixOmin1:ixOmax1) * &
       w(ixOmin1:ixOmax1,rho_)/w(ixOmin1:ixOmax1,&
       press_)* (1.0d0 + w(ixOmin1:ixOmax1,press_)/w(ixOmin1:ixOmax1,&
       rho_) + w(ixOmin1:ixOmax1,eps_))
    do ix1 = ixOmin1,ixOmax1 
    ! pressure is needed to be updated based on the d2w
       call eos_get_pressure_one_grid(d2w_press(ix1), d2w(ix1, rho_), d2w(ix1,&
           eps_), ye=d2w(ix1, ye_),temp= d2w(ix1,temp_))
    enddo

    drho(ixOmin1:ixOmax1) = gamma_eos_tmp(ixOmin1:ixOmax1) * &
       abs(d2w(ixOmin1:ixOmax1,rho_)) / min( w(ixLmin1:ixLmax1, rho_),&
        w(ixRmin1:ixRmax1, rho_) )

    dp(ixOmin1:ixOmax1) = gamma_eos_tmp(ixOmin1:ixOmax1) * &
       abs(d2w(ixOmin1:ixOmax1,press_)) / min( w(ixLmin1:ixLmax1, press_),&
        w(ixRmin1:ixRmax1, press_) )
  end subroutine grhd_ccsn_leakage_ppm_flatcd

  subroutine grhd_ccsn_leakage_ppm_flatsh(ixImin1,ixImax1, ixOmin1,ixOmax1,&
      ixLLmin1,ixLLmax1, ixLmin1,ixLmax1, ixRmin1,ixRmax1, ixRRmin1,ixRRmax1,&
      idims, w, beta, z, dv)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
        ixLLmin1,ixLLmax1, ixLmin1,ixLmax1, ixRmin1,ixRmax1, ixRRmin1,ixRRmax1
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nprim)
    double precision, intent(inout) :: beta(ixGlo1:ixGhi1), z(ixGlo1:ixGhi1),&
        dv(ixGlo1:ixGhi1)

    double precision                :: cs2(ixImin1:ixImax1)
    integer                         :: jxOmin1,jxOmax1, hxOmin1,hxOmax1,&
        ixmin1,ixmax1

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where ( abs( w(ixRRmin1:ixRRmax1, press_) - w(ixLLmin1:ixLLmax1,&
        press_) )> smalldouble )
       beta(ixOmin1:ixOmax1) = abs( w(ixRmin1:ixRmax1,&
           press_) - w(ixLmin1:ixLmax1, press_) ) /( w(ixRRmin1:ixRRmax1,&
           press_) - w(ixLLmin1:ixLLmax1, press_) )
    else where
       beta(ixOmin1:ixOmax1) = 0.0d0
    end where

    ! eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
    call phys_get_csound2(ixImin1,ixImax1, ixOmin1,ixOmax1, w(ixImin1:ixImax1,&
        1:nprim), cs2(ixOmin1:ixOmax1))
    z(ixOmin1:ixOmax1) = abs( w(ixRmin1:ixRmax1, press_) - w(ixLmin1:ixLmax1,&
        press_) ) /( w(ixOmin1:ixOmax1, rho_)*cs2(ixOmin1:ixOmax1) )

    hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
    jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
    ixmin1=ixOmin1-1;ixmax1=ixOmax1+1;

    !dv(ixO^S)=(w(jxO^S,veloc(idims))-w(hxO^S,veloc(idims)) )
    dv(ixOmin1:ixOmax1)=(w(jxOmin1:jxOmax1,W_vel(idims))-w(hxOmin1:hxOmax1,&
       W_vel(idims)) )

  end subroutine grhd_ccsn_leakage_ppm_flatsh

end module mod_grhd_ccsn_leakage_ppm
