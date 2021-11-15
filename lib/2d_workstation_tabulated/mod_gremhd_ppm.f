!> General relativistic hydrodynamics PPM module
module mod_gremhd_ppm
  use mod_physics

  implicit none
  private

  ! Public methods
  public :: gremhd_ppm_init

contains

  !> Initialize the module
  subroutine gremhd_ppm_init()
    use mod_physics_ppm
    phys_ppm_flatcd => gremhd_ppm_flatcd
    phys_ppm_flatsh => gremhd_ppm_flatsh
  end subroutine gremhd_ppm_init

  subroutine gremhd_ppm_flatcd(rec_from, rec_to, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLmin1,ixLmin2,ixLmax1,ixLmax2,&
      ixRmin1,ixRmin2,ixRmax1,ixRmax2, w, d2w, drho, dp)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: rec_from, rec_to, ixImin1,ixImin2,&
       ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim),d2w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,rec_from:rec_to)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
        dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

    double precision                :: gamma_eos_tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2), d2w_press(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: ix1,ix2 

    ! cs2 needed to be updated based on the updated prim vars
    call phys_get_csound2(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim),&
        gamma_eos_tmp(ixImin1:ixImax1,ixImin2:ixImax2))
    gamma_eos_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       gamma_eos_tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       press_)* (1.0d0 + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       press_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) + w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,eps_))
    do ix1 = ixOmin1,ixOmax1 
    do ix2 = ixOmin2,ixOmax2 
    ! pressure is needed to be updated based on the d2w
       call eos_get_pressure_one_grid(d2w_press(ix1,ix2), d2w(ix1,ix2, rho_),&
           d2w(ix1,ix2, eps_))
    enddo
    enddo

    drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = gamma_eos_tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * abs(d2w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)) / min( w(ixLmin1:ixLmax1,ixLmin2:ixLmax2, rho_),&
        w(ixRmin1:ixRmax1,ixRmin2:ixRmax2, rho_) )

    dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = gamma_eos_tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * abs(d2w_press(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)) / min( w(ixLmin1:ixLmax1,ixLmin2:ixLmax2, press_),&
        w(ixRmin1:ixRmax1,ixRmin2:ixRmax2, press_) )
  end subroutine gremhd_ppm_flatcd

  subroutine gremhd_ppm_flatsh(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2, ixLmin1,&
     ixLmin2,ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2, ixRRmin1,&
     ixRRmin2,ixRRmax1,ixRRmax2, idims, w, beta, z, dv)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,&
        ixLmin1,ixLmin2,ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
        ixRRmin1,ixRRmin2,ixRRmax1,ixRRmax2
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(inout) :: beta(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
        z(ixGlo1:ixGhi1,ixGlo2:ixGhi2), dv(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

    double precision                :: cs2(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: jxOmin1,jxOmin2,jxOmax1,jxOmax2,&
        hxOmin1,hxOmin2,hxOmax1,hxOmax2

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where ( abs( w(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
        press_) - w(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2,&
        press_) )> smalldouble )
       beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs( w(ixRmin1:ixRmax1,&
          ixRmin2:ixRmax2, press_) - w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
           press_) ) /( w(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
           press_) - w(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2, press_) )
    else where
       beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
    end where

    ! eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
    call phys_get_csound2(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim),&
        cs2(ixImin1:ixImax1,ixImin2:ixImax2))
    z(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs( w(ixRmin1:ixRmax1,&
       ixRmin2:ixRmax2, press_) - w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
        press_) ) /( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_)*cs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )

    hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
    hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
    jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
    jxOmax1=ixOmax1+kr(idims,1);jxOmax2=ixOmax2+kr(idims,2);

    !dv(ixO^S)=(w(jxO^S,veloc(idims))-w(hxO^S,veloc(idims)) )
    dv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
       W_vel(idims))-w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,W_vel(idims)) )

  end subroutine gremhd_ppm_flatsh

end module mod_gremhd_ppm
