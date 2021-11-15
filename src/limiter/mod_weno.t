module mod_weno
  ! All kinds of (W)ENO schemes
  !
  ! 2019.9.19 WENO(-JS)5 is transplant from the BHAC code by nanami;
  ! 2019.9.20 WENO3 is coded up by nanami;
  ! 2019.9.21 WENO-Z5 is coded up by nanami;
  ! 2019.9.22 WENO-Z+5 is transplant from the BHAC code by nanami;
  ! 2019.10.30 WENO(-JS)7 is coded up by nanami;
  ! 2019.10.31 MPWENO7 is coded up by nanami;
  ! 2019.11.1 exENO7 is code up by nanami;
  ! 2019.11.7 clean up the code, comment out the interpolation variation.
  !
  ! check Jiang & Shu 1996 for the basic idea of WENO;
  ! see Shu 2009 (SIAM review paper) for the implement of WENO5;
  ! see Borges et al. 2008 for the WENO-Z variation;
  ! see Acker et al. 2016 for the WENO-Z+ variation;
  ! see Balsara et al. 2000 for the MPWENO variation and basically WENO7 scheme;
  ! while the extended ENO scheme is now only used for tests, see Xu et al. 2019 for details.
   
  implicit none
  private

  public :: WENO3limiter
  public :: WENO5limiter
  public :: WENO7limiter
  public :: exENO7limiter

contains

  subroutine WENO3limiter(rec_from,rec_to,ixI^L,iL^L,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nprim)
    double precision, intent(inout) :: wRC(ixI^S,1:nprim),wLC(ixI^S,1:nprim) 
    !> local
    integer                         :: iLm^L, iLp^L, iLpp^L
    double precision                :: f_array(ixI^S,1:nprim,2), d_array(2)
    double precision                :: beta(ixI^S,1:nprim,2)
    double precision                :: u1_coeff(2), u2_coeff(2)
    double precision                :: alpha_array(ixI^S,1:nprim,2), alpha_sum(ixI^S,1:nprim), flux(ixI^S,1:nprim)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-12

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    d_array(1:2) = (/ 1.0d0/4.0d0, 3.0d0/4.0d0 /)
    u1_coeff(1:2) = (/ -1.d0/2.d0, 3.d0/2.d0 /)
    u2_coeff(1:2) = (/ 1.d0/2.d0, 1.d0/2.d0 /)
    
    !> left side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLm^S,rec_from:rec_to) + u1_coeff(2) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iL^S,rec_from:rec_to)  + u2_coeff(2) * w(iLp^S,rec_from:rec_to)
  
    beta(iL^S,rec_from:rec_to,1) = (w(iL^S,rec_from:rec_to) - w(iLm^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = (w(iLp^S,rec_from:rec_to) - w(iL^S,rec_from:rec_to))**2
  
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    do i = 1,2
       alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
       alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
    end do
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,2
       flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
  
    !> left value at right interface
    wLC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)
  
    !> right side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLpp^S,rec_from:rec_to) + u1_coeff(2) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLp^S,rec_from:rec_to)  + u2_coeff(2) * w(iL^S,rec_from:rec_to)
  
    beta(iL^S,rec_from:rec_to,1) = (w(iLpp^S,rec_from:rec_to) - w(iLp^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = (w(iLp^S,rec_from:rec_to) - w(iL^S,rec_from:rec_to))**2
  
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    do i = 1,2
       alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
       alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
    end do
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,2
       flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
  
    !> right value at right interface
    wRC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine WENO3limiter

  subroutine WENO5limiter(rec_from,rec_to,ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixI^L, iL^L, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nprim)
    double precision, intent(inout) :: wRC(ixI^S,1:nprim),wLC(ixI^S,1:nprim) 
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,1:nprim,3), d_array(3)
    double precision                :: beta(ixI^S,1:nprim,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nprim), tmp(ixI^S,1:nprim)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nprim,3), alpha_sum(ixI^S,1:nprim), flux(ixI^S,1:nprim)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-18
    double precision                :: lambda
    double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    lambda = dxdim**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> left side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLmm^S,rec_from:rec_to) + u1_coeff(2) * w(iLm^S,rec_from:rec_to) + u1_coeff(3) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLm^S,rec_from:rec_to)  + u2_coeff(2) * w(iL^S,rec_from:rec_to)  + u2_coeff(3) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * w(iL^S,rec_from:rec_to)   + u3_coeff(2) * w(iLp^S,rec_from:rec_to) + u3_coeff(3) * w(iLpp^S,rec_from:rec_to)  
  
    beta(iL^S,rec_from:rec_to,1) = beta_coeff(1) * (w(iLmm^S,rec_from:rec_to) + w(iL^S,rec_from:rec_to) - 2.0d0*w(iLm^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLmm^S,rec_from:rec_to) - 4.0d0 * w(iLm^S,rec_from:rec_to) + 3.0d0*w(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = beta_coeff(1) * (w(iLm^S,rec_from:rec_to) + w(iLp^S,rec_from:rec_to) - 2.0d0 * w(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLm^S,rec_from:rec_to) - w(iLp^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,3) = beta_coeff(1) * (w(iL^S,rec_from:rec_to) + w(iLpp^S,rec_from:rec_to) - 2.0d0 * w(iLp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, rec_from:rec_to) - 4.0d0 * w(iLp^S,rec_from:rec_to) + w(iLpp^S,rec_from:rec_to))**2
 
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
         alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(2)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3))
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.d0 + (tau(iL^S,rec_from:rec_to) / &
                                      (beta(iL^S,rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(3)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3))
      do i = 1,3
        tmp(iL^S,rec_from:rec_to) = (tau(iL^S,rec_from:rec_to) + weno_eps_machine) / (beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.0d0 + tmp(iL^S,rec_from:rec_to)**2 + lambda/tmp(iL^S,rec_from:rec_to))
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    end select

    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
  
    !> left value at right interface
    wLC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)
  
    !> right side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLppp^S,rec_from:rec_to) + u1_coeff(2) * w(iLpp^S,rec_from:rec_to) + u1_coeff(3) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLpp^S,rec_from:rec_to)  + u2_coeff(2) * w(iLp^S,rec_from:rec_to)  + u2_coeff(3) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * w(iLp^S,rec_from:rec_to)   + u3_coeff(2) * w(iL^S,rec_from:rec_to)   + u3_coeff(3) * w(iLm^S,rec_from:rec_to)  
  
    beta(iL^S,rec_from:rec_to,1) = beta_coeff(1) * (w(iLppp^S,rec_from:rec_to) + w(iLp^S,rec_from:rec_to) - 2.0d0*w(iLpp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLppp^S,rec_from:rec_to) - 4.0d0 * w(iLpp^S,rec_from:rec_to) + 3.0d0*w(iLp^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = beta_coeff(1) * (w(iLpp^S,rec_from:rec_to) + w(iL^S,rec_from:rec_to) - 2.0d0 * w(iLp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLpp^S,rec_from:rec_to) - w(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,3) = beta_coeff(1) * (w(iLp^S,rec_from:rec_to) + w(iLm^S,rec_from:rec_to) - 2.0d0 * w(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, rec_from:rec_to) - 4.0d0 * w(iL^S,rec_from:rec_to) + w(iLm^S,rec_from:rec_to))**2
  
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(2) 
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3))
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.d0 + (tau(iL^S,rec_from:rec_to) / &
                                      (beta(iL^S,rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(3)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3))
      do i = 1,3
        tmp(iL^S,rec_from:rec_to) = (tau(iL^S,rec_from:rec_to) + weno_eps_machine) / (beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.0d0 + tmp(iL^S,rec_from:rec_to)**2 + lambda/tmp(iL^S,rec_from:rec_to))
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    end select
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
  
    !> right value at right interface
    wRC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine WENO5limiter

  subroutine WENO7limiter(rec_from,rec_to,ixI^L,iL^L,idims,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: w(ixI^S,1:nprim)
    double precision, intent(inout) :: wRC(ixI^S,1:nprim),wLC(ixI^S,1:nprim) 
    !> local
    integer                         :: iLm^L, iLmm^L, iLmmm^L
    integer                         :: iLp^L, iLpp^L, iLppp^L, iLpppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L

    double precision, dimension(4)  :: d_array, u1_coeff, u2_coeff, u3_coeff, u4_coeff
    double precision, dimension(ixI^S,1:nprim,4)  :: f_array, beta, alpha_array
    double precision, dimension(ixI^S)         :: a, b, c, tmp, tmp2, tmp3
    double precision, dimension(ixI^S,1:nprim)    :: alpha_sum, d, dm4
    double precision, dimension(ixI^S,1:nprim)    :: flux, flux_min, flux_max, flux_ul, flux_md, flux_lc
    integer                         :: i, iw
    double precision, parameter     :: mpalpha = 2.d0, mpbeta = 4.d0
    double precision, parameter     :: weno_eps_machine = 1.0d-18

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLmmm^L=iLmm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iLpppp^L=iLppp^L+kr(idims,^D);

    d_array(1:4) = (/ 1.d0/35.d0, 12.d0/35.d0, 18.d0/35.d0, 4.d0/35.d0 /)
    u1_coeff(1:4) = (/ -1.d0/4.d0, 13.d0/12.d0, -23.d0/12.d0, 25.d0/12.d0 /)
    u2_coeff(1:4) = (/ 1.d0/12.d0, -5.d0/12.d0, 13.d0/12.d0, 1.d0/4.d0 /)
    u3_coeff(1:4) = (/ -1.d0/12.d0, 7.d0/12.d0, 7.d0/12.d0, -1.d0/12.d0 /)
    u4_coeff(1:4) = (/ 1.d0/4.d0, 13.d0/12.d0, -5.d0/12.d0, 1.d0/12.d0 /)
    
    !> left side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLmmm^S,rec_from:rec_to) + u1_coeff(2) * w(iLmm^S,rec_from:rec_to) + u1_coeff(3) * w(iLm^S,rec_from:rec_to) &
                             + u1_coeff(4) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLmm^S,rec_from:rec_to)  + u2_coeff(2) * w(iLm^S,rec_from:rec_to)  + u2_coeff(3) * w(iL^S,rec_from:rec_to)  &
                             + u2_coeff(4) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * w(iLm^S,rec_from:rec_to)   + u3_coeff(2) * w(iL^S,rec_from:rec_to)   + u3_coeff(3) * w(iLp^S,rec_from:rec_to)   &
                             + u3_coeff(4) * w(iLpp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,4) = u4_coeff(1) * w(iL^S,rec_from:rec_to)    + u4_coeff(2) * w(iLp^S,rec_from:rec_to)    + u4_coeff(3) * w(iLpp^S,rec_from:rec_to)  &
                             + u4_coeff(4) * w(iLppp^S,rec_from:rec_to)
  
    beta(iL^S,rec_from:rec_to,1) = w(iLmmm^S,rec_from:rec_to) * (547.d0 * w(iLmmm^S,rec_from:rec_to) - 3882.d0 * w(iLmm^S,rec_from:rec_to) + 4642.d0 * w(iLm^S,rec_from:rec_to) &
                          - 1854.d0 * w(iL^S,rec_from:rec_to)) &
                          + w(iLmm^S,rec_from:rec_to) * (7043.d0 * w(iLmm^S,rec_from:rec_to) - 17246.d0 * w(iLm^S,rec_from:rec_to) + 7042.d0 * w(iL^S,rec_from:rec_to)) &
                          + w(iLm^S,rec_from:rec_to) * (11003.d0 * w(iLm^S,rec_from:rec_to) - 9402.d0 * w(iL^S,rec_from:rec_to)) + 2107.d0 * w(iL^S,rec_from:rec_to)**2
    beta(iL^S,rec_from:rec_to,2) = w(iLmm^S,rec_from:rec_to) * (267.d0 * w(iLmm^S,rec_from:rec_to) - 1642.d0 * w(iLm^S,rec_from:rec_to) + 1602.d0 * w(iL^S,rec_from:rec_to) &
                          - 494.d0 * w(iLp^S,rec_from:rec_to))  &
                          + w(iLm^S,rec_from:rec_to) * (2843.d0 * w(iLm^S,rec_from:rec_to) - 5966.d0 * w(iL^S,rec_from:rec_to) + 1922.d0 * w(iLp^S,rec_from:rec_to)) &
                          + w(iL^S,rec_from:rec_to) * (3443.d0 * w(iL^S,rec_from:rec_to) - 2522.d0 * w(iLp^S,rec_from:rec_to)) + 547.d0 * w(iLp^S,rec_from:rec_to) ** 2
    beta(iL^S,rec_from:rec_to,3) = w(iLm^S,rec_from:rec_to) * (547.d0 * w(iLm^S,rec_from:rec_to) - 2522.d0 * w(iL^S,rec_from:rec_to) + 1922.d0 * w(iLp^S,rec_from:rec_to) &
                          - 494.d0 * w(iLpp^S,rec_from:rec_to))  &
                          + w(iL^S,rec_from:rec_to) * (3443.d0 * w(iL^S,rec_from:rec_to) - 5966.d0 * w(iLp^S,rec_from:rec_to) + 1602.d0 * w(iLpp^S,rec_from:rec_to)) &
                          + w(iLp^S,rec_from:rec_to) * (2843.d0 * w(iLp^S,rec_from:rec_to) - 1642.d0 * w(iLpp^S,rec_from:rec_to)) + 267.d0 * w(iLpp^S,rec_from:rec_to) ** 2
    beta(iL^S,rec_from:rec_to,4) = w(iL^S,rec_from:rec_to) * (2107.d0 * w(iL^S,rec_from:rec_to) - 9402.d0 * w(iLp^S,rec_from:rec_to) + 7042.d0 * w(iLpp^S,rec_from:rec_to) &
                          - 1854.d0 * w(iLppp^S,rec_from:rec_to))  &
                          + w(iLp^S,rec_from:rec_to) * (11003.d0 * w(iLp^S,rec_from:rec_to) - 17246.d0 * w(iLpp^S,rec_from:rec_to) + 4642.d0 * w(iLppp^S,rec_from:rec_to)) &
                          + w(iLpp^S,rec_from:rec_to) * (7043.d0 * w(iLpp^S,rec_from:rec_to) - 3882.d0 * w(iLppp^S,rec_from:rec_to)) &
                          + 547.d0 * w(iLppp^S,rec_from:rec_to) ** 2

    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    do i = 1,4
       alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
       alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
    end do
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,4
       flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
    
    select case(var)
    ! case1 for wenojs, case2 for mpweno
    case(1) 
      wLC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)
    case(2)
      idmax^D=iLmax^D; idmin^D=iLmin^D-kr(idims,^D);
      idm^L=id^L-kr(idims,^D);
      idp^L=id^L+kr(idims,^D);
  
      iemax^D=idmax^D+kr(idims,^D); iemin^D=idmin^D;
      iem^L=ie^L-kr(idims,^D);
      iep^L=ie^L+kr(idims,^D);
  
      d(ie^S,rec_from:rec_to) = w(iep^S,rec_from:rec_to)-2.0d0*w(ie^S,rec_from:rec_to)+w(iem^S,rec_from:rec_to)
  
      do iw=rec_from,rec_to
         a(id^S) = 4.0d0*d(id^S,iw)-d(idp^S,iw)
         b(id^S) = 4.0d0*d(idp^S,iw)-d(id^S,iw)
         call minmod(ixI^L,id^L,a,b,tmp)
         a(id^S) = d(id^S,iw)
         b(id^S) = d(idp^S,iw)
         call minmod(ixI^L,id^L,a,b,tmp2)
         call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
         dm4(id^S,iw) = tmp3(id^S)
      end do

      flux_ul(iL^S,rec_from:rec_to) = w(iL^S,rec_from:rec_to) + mpalpha * (w(iL^S,rec_from:rec_to) - w(iLm^S,rec_from:rec_to))
      flux_md(iL^S,rec_from:rec_to) = half * (w(iL^S,rec_from:rec_to) + w(iLp^S,rec_from:rec_to) - dm4(iL^S,rec_from:rec_to))
      flux_lc(iL^S,rec_from:rec_to) = half * (3.d0 * w(iL^S,rec_from:rec_to) - w(iLm^S,rec_from:rec_to)) + mpbeta / 3.d0 * dm4(iLm^S,rec_from:rec_to)
    
      flux_min(iL^S,rec_from:rec_to) = max(min(w(iL^S,rec_from:rec_to), w(iLp^S,rec_from:rec_to), flux_md(iL^S,rec_from:rec_to)), &
                                min(w(iL^S,rec_from:rec_to), flux_ul(iL^S,rec_from:rec_to),flux_lc(iL^S,rec_from:rec_to)))
  
      flux_max(iL^S,rec_from:rec_to) = min(max(w(iL^S,rec_from:rec_to), w(iLp^S,rec_from:rec_to), flux_md(iL^S,rec_from:rec_to)), &
                                max(w(iL^S,rec_from:rec_to), flux_ul(iL^S,rec_from:rec_to),flux_lc(iL^S,rec_from:rec_to)))
      do iw=rec_from,rec_to
        a(iL^S) = flux(iL^S,iw)
        b(iL^S) = flux_min(iL^S,iw)
        c(iL^S) = flux_max(iL^S,iw)
        call median(ixI^L, iL^L, a, b, c, tmp) 
        wLC(iL^S,iw) = tmp(iL^S)
      end do
    end select

    !> right side
    !>> mmm -> pppp
    !>> mm  -> ppp
    !>> m   -> pp
    !>> 0   -> p
    !>> p   -> 0
    !>> pp  -> m
    !>> ppp -> mm
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLpppp^S,rec_from:rec_to) + u1_coeff(2) * w(iLppp^S,rec_from:rec_to) + u1_coeff(3) * w(iLpp^S,rec_from:rec_to) &
                             + u1_coeff(4) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLppp^S,rec_from:rec_to)  + u2_coeff(2) * w(iLpp^S,rec_from:rec_to)  + u2_coeff(3) * w(iLp^S,rec_from:rec_to)  &
                             + u2_coeff(4) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * w(iLpp^S,rec_from:rec_to)   + u3_coeff(2) * w(iLp^S,rec_from:rec_to)   + u3_coeff(3) * w(iL^S,rec_from:rec_to)   &
                             + u3_coeff(4) * w(iLm^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,4) = u4_coeff(1) * w(iLp^S,rec_from:rec_to)    + u4_coeff(2) * w(iL^S,rec_from:rec_to)    + u4_coeff(3) * w(iLm^S,rec_from:rec_to)  &
                             + u4_coeff(4) * w(iLmm^S,rec_from:rec_to)

    beta(iL^S,rec_from:rec_to,1) = w(iLpppp^S,rec_from:rec_to) * (547.d0 * w(iLpppp^S,rec_from:rec_to) - 3882.d0 * w(iLppp^S,rec_from:rec_to) + 4642.d0 * w(iLpp^S,rec_from:rec_to) &
                          - 1854.d0 * w(iLp^S,rec_from:rec_to)) &
                          + w(iLppp^S,rec_from:rec_to) * (7043.d0 * w(iLppp^S,rec_from:rec_to) - 17246.d0 * w(iLpp^S,rec_from:rec_to) + 7042.d0 * w(iLp^S,rec_from:rec_to)) &
                          + w(iLpp^S,rec_from:rec_to) * (11003.d0 * w(iLpp^S,rec_from:rec_to) - 9402.d0 * w(iLp^S,rec_from:rec_to)) + 2107.d0 * w(iLp^S,rec_from:rec_to)**2
    beta(iL^S,rec_from:rec_to,2) = w(iLppp^S,rec_from:rec_to) * (267.d0 * w(iLppp^S,rec_from:rec_to) - 1642.d0 * w(iLpp^S,rec_from:rec_to) + 1602.d0 * w(iLp^S,rec_from:rec_to) &
                          - 494.d0 * w(iL^S,rec_from:rec_to))  &
                          + w(iLpp^S,rec_from:rec_to) * (2843.d0 * w(iLpp^S,rec_from:rec_to) - 5966.d0 * w(iLp^S,rec_from:rec_to) + 1922.d0 * w(iL^S,rec_from:rec_to)) &
                          + w(iLp^S,rec_from:rec_to) * (3443.d0 * w(iLp^S,rec_from:rec_to) - 2522.d0 * w(iL^S,rec_from:rec_to)) + 547.d0 * w(iL^S,rec_from:rec_to) ** 2
    beta(iL^S,rec_from:rec_to,3) = w(iLpp^S,rec_from:rec_to) * (547.d0 * w(iLpp^S,rec_from:rec_to) - 2522.d0 * w(iLp^S,rec_from:rec_to) + 1922.d0 * w(iL^S,rec_from:rec_to) &
                          - 494.d0 * w(iLm^S,rec_from:rec_to))  &
                          + w(iLp^S,rec_from:rec_to) * (3443.d0 * w(iLp^S,rec_from:rec_to) - 5966.d0 * w(iL^S,rec_from:rec_to) + 1602.d0 * w(iLm^S,rec_from:rec_to)) &
                          + w(iL^S,rec_from:rec_to) * (2843.d0 * w(iL^S,rec_from:rec_to) - 1642.d0 * w(iLm^S,rec_from:rec_to)) + 267.d0 * w(iLm^S,rec_from:rec_to) ** 2
    beta(iL^S,rec_from:rec_to,4) = w(iLp^S,rec_from:rec_to) * (2107.d0 * w(iLp^S,rec_from:rec_to) - 9402.d0 * w(iL^S,rec_from:rec_to) + 7042.d0 * w(iLm^S,rec_from:rec_to) &
                          - 1854.d0 * w(iLmm^S,rec_from:rec_to))  &
                          + w(iL^S,rec_from:rec_to) * (11003.d0 * w(iL^S,rec_from:rec_to) - 17246.d0 * w(iLm^S,rec_from:rec_to) + 4642.d0 * w(iLmm^S,rec_from:rec_to)) &
                          + w(iLm^S,rec_from:rec_to) * (7043.d0 * w(iLm^S,rec_from:rec_to) - 3882.d0 * w(iLmm^S,rec_from:rec_to)) + 547.d0 * w(iLmm^S,rec_from:rec_to) ** 2

    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    do i = 1,4
       alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
       alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
    end do
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,4
       flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do

    select case(var)
    case(1)
      wRC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)
    case(2)
      idmax^D=iLmax^D+kr(idims,^D); idmin^D=iLmin^D;
      idm^L=id^L-kr(idims,^D);
      idp^L=id^L+kr(idims,^D);
  
      iemax^D=idmax^D; iemin^D=idmin^D-kr(idims,^D);
      iem^L=ie^L-kr(idims,^D);
      iep^L=ie^L+kr(idims,^D);
      iepp^L=iep^L+kr(idims,^D);
  
      d(ie^S,rec_from:rec_to) = w(ie^S,rec_from:rec_to)-2.0d0*w(iep^S,rec_from:rec_to)+w(iepp^S,rec_from:rec_to)
  
      do iw=rec_from,rec_to
        a(id^S) = 4.0d0*d(id^S,iw)-d(idm^S,iw)
        b(id^S) = 4.0d0*d(idm^S,iw)-d(id^S,iw)
        call minmod(ixI^L,id^L,a,b,tmp)
        a(id^S) = d(id^S,iw)
        b(id^S) = d(idm^S,iw)
        call minmod(ixI^L,id^L,a,b,tmp2)
        call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
        dm4(id^S,iw) = tmp3(id^S)
      end do
   
      flux_ul(iL^S,rec_from:rec_to) = w(iLp^S,rec_from:rec_to) + mpalpha * (w(iLp^S,rec_from:rec_to) - w(iLpp^S,rec_from:rec_to))
      flux_md(iL^S,rec_from:rec_to) = half * (w(iL^S,rec_from:rec_to) + w(iLp^S,rec_from:rec_to) - dm4(iL^S,rec_from:rec_to))
      flux_lc(iL^S,rec_from:rec_to) = half * (3.d0 * w(iLp^S,rec_from:rec_to) - w(iLpp^S,rec_from:rec_to)) + mpbeta / 3.d0 * dm4(iLp^S,rec_from:rec_to)
    
      flux_min(iL^S,rec_from:rec_to) = max(min(w(iLp^S,rec_from:rec_to), w(iL^S,rec_from:rec_to), flux_md(iL^S,rec_from:rec_to)), &
                                min(w(iLp^S,rec_from:rec_to), flux_ul(iL^S,rec_from:rec_to),flux_lc(iL^S,rec_from:rec_to)))
  
      flux_max(iL^S,rec_from:rec_to) = min(max(w(iLp^S,rec_from:rec_to), w(iL^S,rec_from:rec_to), flux_md(iL^S,rec_from:rec_to)), &
                                max(w(iLp^S,rec_from:rec_to), flux_ul(iL^S,rec_from:rec_to),flux_lc(iL^S,rec_from:rec_to)))
      do iw=rec_from,rec_to
        a(iL^S) = flux(iL^S,iw)
        b(iL^S) = flux_min(iL^S,iw)
        c(iL^S) = flux_max(iL^S,iw)
        call median(ixI^L, iL^L, a, b, c, tmp) 
        wRC(iL^S,iw) = tmp(iL^S)
      end do
    end select
  end subroutine WENO7limiter

  subroutine exENO7limiter(rec_from,rec_to,ixI^L,iL^L,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nprim)
    double precision, intent(inout) :: wRC(ixI^S,1:nprim),wLC(ixI^S,1:nprim) 
    !> local
    integer                         :: i, iw
    integer                         :: iLm^L, iLmm^L, iLmmm^L
    integer                         :: iLp^L, iLpp^L, iLppp^L, iLpppp^L
    integer                         :: iM^L, iMm^L, iMmm^L
    integer                         :: iMp^L, iMpp^L, iMppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer, dimension(ixI^S,1:nprim)  :: delta_sum
    integer, dimension(ixI^S,1:nprim,3):: delta
    double precision, dimension(2)            :: beta_coeff
    double precision, dimension(3)            :: d_array
    double precision, dimension(ixI^S,1:nprim)   :: gamma_sum, flux
    double precision, dimension(ixI^S,1:nprim,3) :: beta, gamma_array, kai_array
    double precision, parameter               :: exeno_ct = 1.d-1
    double precision, parameter               :: weno_eps_machine = 1.d-18

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLmmm^L=iLmm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iLpppp^L=iLppp^L+kr(idims,^D);

    iMmin^D=iLmin^D-kr(idims,^D);
    iMmax^D=iLmax^D+kr(idims,^D);
    iMm^L=iM^L-kr(idims,^D);
    iMmm^L=iMm^L-kr(idims,^D);
    iMp^L=iM^L+kr(idims,^D);
    iMpp^L=iMp^L+kr(idims,^D);
    iMppp^L=iMpp^L+kr(idims,^D);

    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)

    !>> left side
    beta(iM^S,rec_from:rec_to,1) = beta_coeff(1) * (w(iMmm^S,rec_from:rec_to) + w(iM^S,rec_from:rec_to) - 2.0d0 * w(iMm^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iMmm^S,rec_from:rec_to) - 4.0d0 * w(iMm^S,rec_from:rec_to) + 3.0d0 * w(iM^S,rec_from:rec_to))**2
    beta(iM^S,rec_from:rec_to,2) = beta_coeff(1) * (w(iMm^S,rec_from:rec_to) + w(iMp^S,rec_from:rec_to) - 2.0d0 * w(iM^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iMm^S,rec_from:rec_to) - w(iMp^S,rec_from:rec_to))**2
    beta(iM^S,rec_from:rec_to,3) = beta_coeff(1) * (w(iM^S,rec_from:rec_to) + w(iMpp^S,rec_from:rec_to) - 2.0d0 * w(iMp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * w(iM^S, rec_from:rec_to) - 4.0d0 * w(iMp^S,rec_from:rec_to) + w(iMpp^S,rec_from:rec_to))**2

    gamma_sum(iM^S,rec_from:rec_to) = 0.0d0 
    do i = 1,3
      gamma_array(iM^S,rec_from:rec_to,i) = d_array(i) / (beta(iM^S,rec_from:rec_to,i) + weno_eps_machine)**2
      gamma_sum(iM^S,rec_from:rec_to) = gamma_sum(iM^S,rec_from:rec_to) + gamma_array(iM^S,rec_from:rec_to,i)
    end do
    do i = 1,3
      kai_array(iM^S,rec_from:rec_to,i) = gamma_array(iM^S,rec_from:rec_to,i) / gamma_sum(iM^S,rec_from:rec_to)
      where(kai_array(iM^S,rec_from:rec_to,i) .lt. exeno_ct) 
        delta(iM^S,rec_from:rec_to,i) = 0
      elsewhere
        delta(iM^S,rec_from:rec_to,i) = 1
      endwhere
    end do

    delta_sum(iL^S,rec_from:rec_to) = delta(iLm^S,rec_from:rec_to,1) * 1 + delta(iLp^S,rec_from:rec_to,3) * 2 + delta(iL^S,rec_from:rec_to,1) * 4 &
                             + delta(iL^S,rec_from:rec_to,2)  * 8 + delta(iL^S,rec_from:rec_to,3) * 16

    !> f3
    where(delta_sum(iL^S,rec_from:rec_to) .eq. 31)
      flux(iL^S,rec_from:rec_to) = (- 3.d0 * w(iLmmm^S,rec_from:rec_to) + 25.d0 * w(iLmm^S,rec_from:rec_to) - 101.d0 * w(iLm^S,rec_from:rec_to) + 319.d0 * w(iL^S,rec_from:rec_to) &
                             + 214.d0 * w(iLp^S,rec_from:rec_to) - 38.d0 * w(iLpp^S,rec_from:rec_to) + 4.d0 * w(iLppp^S,rec_from:rec_to)) / 420.d0
    !> f4
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .eq. 30)
      flux(iL^S,rec_from:rec_to) = (w(iLmm^S,rec_from:rec_to) - 8.d0 * w(iLm^S,rec_from:rec_to) + 37.d0 * w(iL^S,rec_from:rec_to) + 37.d0 * w(iLp^S,rec_from:rec_to) &
                             - 8.d0 * w(iLpp^S,rec_from:rec_to) + w(iLppp^S,rec_from:rec_to)) / 60.d0
    !> f5
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .eq. 29)
      flux(iL^S,rec_from:rec_to) = (- w(iLmmm^S,rec_from:rec_to) + 7.d0 * w(iLmm^S,rec_from:rec_to) - 23.d0 * w(iLm^S,rec_from:rec_to) + 57.d0 * w(iL^S,rec_from:rec_to) &
                             + 22.d0 * w(iLp^S,rec_from:rec_to) - 2.d0 * w(iLpp^S,rec_from:rec_to)) / 60.d0
    !> f6
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .eq. 28)
      flux(iL^S,rec_from:rec_to) = (2.d0 * w(iLmm^S,rec_from:rec_to) - 13.d0 * w(iLm^S,rec_from:rec_to) + 47.d0 * w(iL^S,rec_from:rec_to) + 27.d0 * w(iLp^S,rec_from:rec_to) &
                             - 3.d0 * w(iLpp^S,rec_from:rec_to)) / 60.d0
    !> f7
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 24)
      flux(iL^S,rec_from:rec_to) = (- w(iLm^S,rec_from:rec_to) + 7.d0 * w(iL^S,rec_from:rec_to) + 7.d0 * w(iLp^S,rec_from:rec_to) - w(iLpp^S,rec_from:rec_to)) / 12.d0
    !> f9
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 16)
      flux(iL^S,rec_from:rec_to) = (2.d0 * w(iL^S,rec_from:rec_to) + 5.d0 * w(iLp^S,rec_from:rec_to) - w(iLpp^S,rec_from:rec_to)) / 6.d0
    !> f8
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 12)
      flux(iL^S,rec_from:rec_to) = (w(iLmm^S,rec_from:rec_to) - 5.d0 * w(iLm^S,rec_from:rec_to) + 13.d0 * w(iL^S,rec_from:rec_to) + 3.d0 * w(iLp^S,rec_from:rec_to)) / 12.d0
    !> f10
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 8)
      flux(iL^S,rec_from:rec_to) = (- w(iLm^S,rec_from:rec_to) + 5.d0 * w(iL^S,rec_from:rec_to) + 2.d0 * w(iLp^S,rec_from:rec_to)) / 6.d0
    !> f11
    elsewhere
     flux(iL^S,rec_from:rec_to) = (2.d0 * w(iLmm^S,rec_from:rec_to) - 7.d0 * w(iLm^S,rec_from:rec_to) + 11.d0 * w(iL^S,rec_from:rec_to)) / 6.d0
    endwhere

    wLC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

    !> right side
    beta(iM^S,rec_from:rec_to,1) = beta_coeff(1) * (w(iMppp^S,rec_from:rec_to) + w(iMp^S,rec_from:rec_to) - 2.0d0 * w(iMpp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iMppp^S,rec_from:rec_to) - 4.0d0 * w(iMpp^S,rec_from:rec_to) + 3.0d0 * w(iMp^S,rec_from:rec_to))**2
    beta(iM^S,rec_from:rec_to,2) = beta_coeff(1) * (w(iMpp^S,rec_from:rec_to) + w(iM^S,rec_from:rec_to) - 2.0d0 * w(iMp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iMpp^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to))**2
    beta(iM^S,rec_from:rec_to,3) = beta_coeff(1) * (w(iMp^S,rec_from:rec_to) + w(iMm^S,rec_from:rec_to) - 2.0d0 * w(iM^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * w(iMp^S, rec_from:rec_to) - 4.0d0 * w(iM^S,rec_from:rec_to) + w(iMm^S,rec_from:rec_to))**2

    gamma_sum(iM^S,rec_from:rec_to) = 0.0d0 
    do i = 1,3
      gamma_array(iM^S,rec_from:rec_to,i) = d_array(i) / (beta(iM^S,rec_from:rec_to,i) + weno_eps_machine)**2
      gamma_sum(iM^S,rec_from:rec_to) = gamma_sum(iM^S,rec_from:rec_to) + gamma_array(iM^S,rec_from:rec_to,i)
    end do
    do i = 1,3
      kai_array(iM^S,rec_from:rec_to,i) = gamma_array(iM^S,rec_from:rec_to,i) / gamma_sum(iM^S,rec_from:rec_to)
      where(kai_array(iM^S,rec_from:rec_to,i) .lt. exeno_ct) 
        delta(iM^S,rec_from:rec_to,i) = 0
      elsewhere
        delta(iM^S,rec_from:rec_to,i) = 1
      endwhere
    end do
 
    delta_sum(iL^S,rec_from:rec_to) = delta(iLp^S,rec_from:rec_to,1) * 1 + delta(iLm^S,rec_from:rec_to,3) * 2 + delta(iL^S,rec_from:rec_to,1) * 4 &
                             + delta(iL^S,rec_from:rec_to,2)  * 8 + delta(iL^S,rec_from:rec_to,3) * 16

    where(delta_sum(iL^S,rec_from:rec_to) .eq. 31)
      flux(iL^S,rec_from:rec_to) = (- 3.d0 * w(iLpppp^S,rec_from:rec_to) + 25.d0 * w(iLppp^S,rec_from:rec_to) - 101.d0 * w(iLpp^S,rec_from:rec_to) &
                             + 319.d0 * w(iLp^S,rec_from:rec_to) + 214.d0 * w(iL^S,rec_from:rec_to) - 38.d0 * w(iLm^S,rec_from:rec_to) &
                             + 4.d0 * w(iLmm^S,rec_from:rec_to)) / 420.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .eq. 30)
      flux(iL^S,rec_from:rec_to) = (w(iLppp^S,rec_from:rec_to) - 8.d0 * w(iLpp^S,rec_from:rec_to) + 37.d0 * w(iLp^S,rec_from:rec_to) + 37.d0 * w(iL^S,rec_from:rec_to) &
                             - 8.d0 * w(iLm^S,rec_from:rec_to) + w(iLmm^S,rec_from:rec_to)) / 60.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .eq. 29)
      flux(iL^S,rec_from:rec_to) = (- w(iLpppp^S,rec_from:rec_to) + 7.d0 * w(iLppp^S,rec_from:rec_to) - 23.d0 * w(iLpp^S,rec_from:rec_to) + 57.d0 * w(iLp^S,rec_from:rec_to) &
                             + 22.d0 * w(iL^S,rec_from:rec_to) - 2.d0 * w(iLm^S,rec_from:rec_to)) / 60.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .eq. 28)
      flux(iL^S,rec_from:rec_to) = (2.d0 * w(iLppp^S,rec_from:rec_to) - 13.d0 * w(iLpp^S,rec_from:rec_to) + 47.d0 * w(iLp^S,rec_from:rec_to) + 27.d0 * w(iL^S,rec_from:rec_to) &
                             - 3.d0 * w(iLm^S,rec_from:rec_to)) / 60.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 24)
      flux(iL^S,rec_from:rec_to) = (- w(iLpp^S,rec_from:rec_to) + 7.d0 * w(iLp^S,rec_from:rec_to) + 7.d0 * w(iL^S,rec_from:rec_to) - w(iLm^S,rec_from:rec_to)) / 12.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 16)
      flux(iL^S,rec_from:rec_to) = (2.d0 * w(iLp^S,rec_from:rec_to) + 5.d0 * w(iL^S,rec_from:rec_to) - w(iLm^S,rec_from:rec_to)) / 6.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 12)
      flux(iL^S,rec_from:rec_to) = (w(iLppp^S,rec_from:rec_to) - 5.d0 * w(iLpp^S,rec_from:rec_to) + 13.d0 * w(iLp^S,rec_from:rec_to) + 3.d0 * w(iL^S,rec_from:rec_to)) / 12.d0
    elsewhere(delta_sum(iL^S,rec_from:rec_to) .ge. 8)
      flux(iL^S,rec_from:rec_to) = (- w(iLpp^S,rec_from:rec_to) + 5.d0 * w(iLp^S,rec_from:rec_to) + 2.d0 * w(iL^S,rec_from:rec_to)) / 6.d0
    elsewhere
     flux(iL^S,rec_from:rec_to) = (2.d0 * w(iLppp^S,rec_from:rec_to) - 7.d0 * w(iLpp^S,rec_from:rec_to) + 11.d0 * w(iLp^S,rec_from:rec_to)) / 6.d0
    endwhere

    wRC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine exENO7limiter

  subroutine minmod(ixI^L,ixO^L,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S)
    double precision, intent(out):: minm(ixI^S)

    minm(ixO^S) = (sign(one,a(ixO^S))+sign(one,b(ixO^S)))/2.0d0 * &
         min(abs(a(ixO^S)),abs(b(ixO^S)))

  end subroutine minmod

  subroutine median(ixI^L,ixO^L,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S), c(ixI^S)
    double precision, intent(out):: med(ixI^S)

    double precision             :: tmp1(ixI^S),tmp2(ixI^S)

    tmp1(ixO^S) = b(ixO^S) - a(ixO^S); tmp2(ixO^S) = c(ixO^S) - a(ixO^S)

    med(ixO^S) = a(ixO^S) + (sign(one,tmp1(ixO^S))+sign(one,tmp2(ixO^S)))/2.0d0 * &
         min(abs(tmp1(ixO^S)),abs(tmp2(ixO^S)))

  end subroutine median
end module mod_weno
