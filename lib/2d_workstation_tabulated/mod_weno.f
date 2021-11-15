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

  subroutine WENO3limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
     iLmin1,iLmin2,iLmax1,iLmax2,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLpmin1,iLpmin2,iLpmax1,iLpmax2, iLppmin1,iLppmin2,iLppmax1,iLppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim,2), d_array(2)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim,2)
    double precision                :: u1_coeff(2), u2_coeff(2)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nprim,2), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-12

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    d_array(1:2) = (/ 1.0d0/4.0d0, 3.0d0/4.0d0 /)
    u1_coeff(1:2) = (/ -1.d0/2.d0, 3.d0/2.d0 /)
    u2_coeff(1:2) = (/ 1.d0/2.d0, 1.d0/2.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = u1_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + u1_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = u2_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)  + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,1) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,2) = (w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) - w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))**2
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0 
    do i = 1,2
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0
    do i = 1,2
       flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,rec_from:rec_to) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to))
    end do
  
    !> left value at right interface
    wLC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)
  
    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = u1_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + u1_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = u2_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,1) = (w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,2) = (w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) - w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))**2
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0 
    do i = 1,2
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0
    do i = 1,2
       flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,rec_from:rec_to) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to))
    end do
  
    !> right value at right interface
    wRC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)

  end subroutine WENO3limiter

  subroutine WENO5limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
     iLmin1,iLmin2,iLmax1,iLmax2,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nprim,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-18
    double precision                :: lambda
    double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
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
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = u1_coeff(1) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) + u1_coeff(2) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + u1_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = u2_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)  + u2_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       3) = u3_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) + u3_coeff(3) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = beta_coeff(1) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,rec_from:rec_to) - 4.0d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to) + 3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = beta_coeff(1) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       3) = beta_coeff(1) * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (3.0d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2, rec_from:rec_to) - 4.0d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) + w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to))**2
 
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
            i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
            i) + weno_eps_machine)**2
         alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
            rec_from:rec_to,i)
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
         abs(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
         1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,3))
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
         abs(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
         1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,3))
      do i = 1,3
        tmp(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (tau(iLmin1:iLmax1,&
           iLmin2:iLmax2,rec_from:rec_to) + weno_eps_machine) / &
           (beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) + weno_eps_machine)
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to))
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i)
      end do
    end select

    flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to))
    end do
  
    !> left value at right interface
    wLC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)
  
    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = u1_coeff(1) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) + u1_coeff(2) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + u1_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = u2_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)  + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)  + u2_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       3) = u3_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)   + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)   + u3_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = beta_coeff(1) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,rec_from:rec_to) - 4.0d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) + 3.0d0*w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = beta_coeff(1) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) - w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       3) = beta_coeff(1) * (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2, rec_from:rec_to) - 4.0d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to))**2
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) + weno_eps_machine)**2
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i)
      end do
    case(2) 
      tau(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
         abs(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
         1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,3))
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
         abs(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
         1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,3))
      do i = 1,3
        tmp(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (tau(iLmin1:iLmax1,&
           iLmin2:iLmax2,rec_from:rec_to) + weno_eps_machine) / &
           (beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) + weno_eps_machine)
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
           i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to))
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
           rec_from:rec_to,i)
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to))
    end do
  
    !> right value at right interface
    wRC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)

  end subroutine WENO5limiter

  subroutine WENO7limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
     iLmin1,iLmin2,iLmax1,iLmax2,idims,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims, var
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLmmmmin1,iLmmmmin2,iLmmmmax1,&
       iLmmmmax2
    integer                         :: iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2, iLppppmin1,iLppppmin2,iLppppmax1,iLppppmax2
    integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
       idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
       idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
       iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
       ieppmax1,ieppmax2

    double precision, dimension(4)  :: d_array, u1_coeff, u2_coeff, u3_coeff,&
        u4_coeff
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim,&
       4)  :: f_array, beta, alpha_array
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)         :: a,&
        b, c, tmp, tmp2, tmp3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)    :: alpha_sum, d, dm4
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)    :: flux, flux_min, flux_max, flux_ul, flux_md, flux_lc
    integer                         :: i, iw
    double precision, parameter     :: mpalpha = 2.d0, mpbeta = 4.d0
    double precision, parameter     :: weno_eps_machine = 1.0d-18

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLmmmmin1=iLmmmin1-kr(idims,1);iLmmmmin2=iLmmmin2-kr(idims,2)
    iLmmmmax1=iLmmmax1-kr(idims,1);iLmmmmax2=iLmmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    iLppppmin1=iLpppmin1+kr(idims,1);iLppppmin2=iLpppmin2+kr(idims,2)
    iLppppmax1=iLpppmax1+kr(idims,1);iLppppmax2=iLpppmax2+kr(idims,2);

    d_array(1:4) = (/ 1.d0/35.d0, 12.d0/35.d0, 18.d0/35.d0, 4.d0/35.d0 /)
    u1_coeff(1:4) = (/ -1.d0/4.d0, 13.d0/12.d0, -23.d0/12.d0, 25.d0/12.d0 /)
    u2_coeff(1:4) = (/ 1.d0/12.d0, -5.d0/12.d0, 13.d0/12.d0, 1.d0/4.d0 /)
    u3_coeff(1:4) = (/ -1.d0/12.d0, 7.d0/12.d0, 7.d0/12.d0, -1.d0/12.d0 /)
    u4_coeff(1:4) = (/ 1.d0/4.d0, 13.d0/12.d0, -5.d0/12.d0, 1.d0/12.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = u1_coeff(1) * w(iLmmmmin1:iLmmmmax1,iLmmmmin2:iLmmmmax2,&
       rec_from:rec_to) + u1_coeff(2) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) + u1_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + u1_coeff(4) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = u2_coeff(1) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to)  + u2_coeff(2) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)  + u2_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)  + u2_coeff(4) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       3) = u3_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)   + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)   + u3_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)   + u3_coeff(4) * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       4) = u4_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)    + u4_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)    + u4_coeff(3) * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to)  + u4_coeff(4) * &
       w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,rec_from:rec_to)
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = w(iLmmmmin1:iLmmmmax1,iLmmmmin2:iLmmmmax2,&
       rec_from:rec_to) * (547.d0 * w(iLmmmmin1:iLmmmmax1,iLmmmmin2:iLmmmmax2,&
       rec_from:rec_to) - 3882.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) + 4642.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) - 1854.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)) + w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) * (7043.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) - 17246.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + 7042.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) * (11003.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) - 9402.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)) + 2107.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,2) = w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,rec_from:rec_to) * (267.d0 * w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,rec_from:rec_to) - 1642.d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to) + 1602.d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) - 494.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to))  + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) * (2843.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) - 5966.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) + 1922.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) * (3443.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 2522.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)) + 547.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,3) = w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to) * (547.d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to) - 2522.d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) + 1922.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) - 494.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to))  + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) * (3443.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 5966.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) + 1602.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) * (2843.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 1642.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)) + 267.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,4) = w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) * (2107.d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) - 9402.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) + 7042.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) - 1854.d0 * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,rec_from:rec_to))  + w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) * (11003.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) - 17246.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) + 4642.d0 * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,rec_from:rec_to)) + w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) * (7043.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) - 3882.d0 * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,rec_from:rec_to)) + 547.d0 * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,rec_from:rec_to) ** 2

    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0 
    do i = 1,4
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0
    do i = 1,4
       flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,rec_from:rec_to) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to))
    end do
    
    select case(var)
    ! case1 for wenojs, case2 for mpweno
    case(1) 
      wLC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to)
    case(2)
      idmax1=iLmax1;idmax2=iLmax2; idmin1=iLmin1-kr(idims,1)
      idmin2=iLmin2-kr(idims,2);
      idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
      idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
      idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
      idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);
  
      iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2); iemin1=idmin1
      iemin2=idmin2;
      iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
      iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
      iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
      iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);
  
      d(iemin1:iemax1,iemin2:iemax2,rec_from:rec_to) = w(iepmin1:iepmax1,&
         iepmin2:iepmax2,rec_from:rec_to)-2.0d0*w(iemin1:iemax1,iemin2:iemax2,&
         rec_from:rec_to)+w(iemmin1:iemmax1,iemmin2:iemmax2,rec_from:rec_to)
  
      do iw=rec_from,rec_to
         a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,&
            iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
         b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idpmin1:idpmax1,&
            idpmin2:idpmax2,iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
         call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
            idmax2,a,b,tmp)
         a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
         b(idmin1:idmax1,idmin2:idmax2) = d(idpmin1:idpmax1,idpmin2:idpmax2,&
            iw)
         call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
            idmax2,a,b,tmp2)
         call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
            idmax2,tmp,tmp2,tmp3)
         dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,&
            idmin2:idmax2)
      end do

      flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = w(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + mpalpha * (w(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to))
      flux_md(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = half * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - dm4(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to))
      flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = half * (3.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to)) + mpbeta / 3.d0 * dm4(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,rec_from:rec_to)
    
      flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = max(min(w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to), w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),&
          flux_md(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)),&
          min(w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
          flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
         flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))
  
      flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = min(max(w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to), w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),&
          flux_md(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)),&
          max(w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
          flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
         flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))
      do iw=rec_from,rec_to
        a(iLmin1:iLmax1,iLmin2:iLmax2) = flux(iLmin1:iLmax1,iLmin2:iLmax2,iw)
        b(iLmin1:iLmax1,iLmin2:iLmax2) = flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        c(iLmin1:iLmax1,iLmin2:iLmax2) = flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        call median(ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,iLmin2,iLmax1,&
           iLmax2, a, b, c, tmp) 
        wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
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
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = u1_coeff(1) * w(iLppppmin1:iLppppmax1,iLppppmin2:iLppppmax2,&
       rec_from:rec_to) + u1_coeff(2) * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,rec_from:rec_to) + u1_coeff(3) * &
       w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + u1_coeff(4) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = u2_coeff(1) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to)  + u2_coeff(2) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)  + u2_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)  + u2_coeff(4) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       3) = u3_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)   + u3_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)   + u3_coeff(4) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       4) = u4_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)    + u4_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)    + u4_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)  + u4_coeff(4) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to)

    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       1) = w(iLppppmin1:iLppppmax1,iLppppmin2:iLppppmax2,&
       rec_from:rec_to) * (547.d0 * w(iLppppmin1:iLppppmax1,&
       iLppppmin2:iLppppmax2,rec_from:rec_to) - 3882.d0 * &
       w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) + 4642.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) - 1854.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)) + w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) * (7043.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) - 17246.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + 7042.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) * (11003.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) - 9402.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)) + 2107.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
       2) = w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) * (267.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) - 1642.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + 1602.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 494.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))  + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) * (2843.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) - 5966.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) + 1922.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) * (3443.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 2522.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)) + 547.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,3) = w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) * (547.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,rec_from:rec_to) - 2522.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) + 1922.d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) - 494.d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to))  + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) * (3443.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 5966.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) + 1602.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) * (2843.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 1642.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)) + 267.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,4) = w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) * (2107.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to) - 9402.d0 * w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to) + 7042.d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,rec_from:rec_to) - 1854.d0 * w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,rec_from:rec_to))  + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) * (11003.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 17246.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + 4642.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to)) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) * (7043.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) - 3882.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to)) + 547.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) ** 2

    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0 
    do i = 1,4
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 0.0d0
    do i = 1,4
       flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,rec_from:rec_to) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to,i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to))
    end do

    select case(var)
    case(1)
      wRC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to)
    case(2)
      idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2); idmin1=iLmin1
      idmin2=iLmin2;
      idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
      idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
      idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
      idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);
  
      iemax1=idmax1;iemax2=idmax2; iemin1=idmin1-kr(idims,1)
      iemin2=idmin2-kr(idims,2);
      iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
      iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
      iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
      iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);
      ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
      ieppmax1=iepmax1+kr(idims,1);ieppmax2=iepmax2+kr(idims,2);
  
      d(iemin1:iemax1,iemin2:iemax2,rec_from:rec_to) = w(iemin1:iemax1,&
         iemin2:iemax2,rec_from:rec_to)-2.0d0*w(iepmin1:iepmax1,&
         iepmin2:iepmax2,rec_from:rec_to)+w(ieppmin1:ieppmax1,&
         ieppmin2:ieppmax2,rec_from:rec_to)
  
      do iw=rec_from,rec_to
        a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,&
           iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
        b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmmin1:idmmax1,&
           idmmin2:idmmax2,iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
        call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
           idmax2,a,b,tmp)
        a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
        b(idmin1:idmax1,idmin2:idmax2) = d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
        call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
           idmax2,a,b,tmp2)
        call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
           idmax2,tmp,tmp2,tmp3)
        dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,&
           idmin2:idmax2)
      end do
   
      flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to) + mpalpha * (w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to) - w(iLppmin1:iLppmax1,&
         iLppmin2:iLppmax2,rec_from:rec_to))
      flux_md(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = half * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - dm4(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to))
      flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = half * (3.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to)) + mpbeta / 3.d0 * dm4(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to)
    
      flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to), w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
          flux_md(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)),&
          min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),&
          flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
         flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))
  
      flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) = min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to), w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
          flux_md(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)),&
          max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),&
          flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
         flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))
      do iw=rec_from,rec_to
        a(iLmin1:iLmax1,iLmin2:iLmax2) = flux(iLmin1:iLmax1,iLmin2:iLmax2,iw)
        b(iLmin1:iLmax1,iLmin2:iLmax2) = flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        c(iLmin1:iLmax1,iLmin2:iLmax2) = flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        call median(ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,iLmin2,iLmax1,&
           iLmax2, a, b, c, tmp) 
        wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
      end do
    end select
  end subroutine WENO7limiter

  subroutine exENO7limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
     iLmin1,iLmin2,iLmax1,iLmax2,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim) 
    !> local
    integer                         :: i, iw
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLmmmmin1,iLmmmmin2,iLmmmmax1,&
       iLmmmmax2
    integer                         :: iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2, iLppppmin1,iLppppmin2,iLppppmax1,iLppppmax2
    integer                         :: iMmin1,iMmin2,iMmax1,iMmax2, iMmmin1,&
       iMmmin2,iMmmax1,iMmmax2, iMmmmin1,iMmmmin2,iMmmmax1,iMmmmax2
    integer                         :: iMpmin1,iMpmin2,iMpmax1,iMpmax2,&
        iMppmin1,iMppmin2,iMppmax1,iMppmax2, iMpppmin1,iMpppmin2,iMpppmax1,&
       iMpppmax2
    integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
       idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
       idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
       iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
       ieppmax1,ieppmax2
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)  :: delta_sum
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim,3):: delta
    double precision, dimension(2)            :: beta_coeff
    double precision, dimension(3)            :: d_array
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)   :: gamma_sum, flux
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim,&
       3) :: beta, gamma_array, kai_array
    double precision, parameter               :: exeno_ct = 1.d-1
    double precision, parameter               :: weno_eps_machine = 1.d-18

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLmmmmin1=iLmmmin1-kr(idims,1);iLmmmmin2=iLmmmin2-kr(idims,2)
    iLmmmmax1=iLmmmax1-kr(idims,1);iLmmmmax2=iLmmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    iLppppmin1=iLpppmin1+kr(idims,1);iLppppmin2=iLpppmin2+kr(idims,2)
    iLppppmax1=iLpppmax1+kr(idims,1);iLppppmax2=iLpppmax2+kr(idims,2);

    iMmin1=iLmin1-kr(idims,1);iMmin2=iLmin2-kr(idims,2);
    iMmax1=iLmax1+kr(idims,1);iMmax2=iLmax2+kr(idims,2);
    iMmmin1=iMmin1-kr(idims,1);iMmmin2=iMmin2-kr(idims,2)
    iMmmax1=iMmax1-kr(idims,1);iMmmax2=iMmax2-kr(idims,2);
    iMmmmin1=iMmmin1-kr(idims,1);iMmmmin2=iMmmin2-kr(idims,2)
    iMmmmax1=iMmmax1-kr(idims,1);iMmmmax2=iMmmax2-kr(idims,2);
    iMpmin1=iMmin1+kr(idims,1);iMpmin2=iMmin2+kr(idims,2)
    iMpmax1=iMmax1+kr(idims,1);iMpmax2=iMmax2+kr(idims,2);
    iMppmin1=iMpmin1+kr(idims,1);iMppmin2=iMpmin2+kr(idims,2)
    iMppmax1=iMpmax1+kr(idims,1);iMppmax2=iMpmax2+kr(idims,2);
    iMpppmin1=iMppmin1+kr(idims,1);iMpppmin2=iMppmin2+kr(idims,2)
    iMpppmax1=iMppmax1+kr(idims,1);iMpppmax2=iMppmax2+kr(idims,2);

    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)

    !>> left side
    beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
       1) = beta_coeff(1) * (w(iMmmmin1:iMmmmax1,iMmmmin2:iMmmmax2,&
       rec_from:rec_to) + w(iMmin1:iMmax1,iMmin2:iMmax2,&
       rec_from:rec_to) - 2.0d0 * w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iMmmmin1:iMmmmax1,&
       iMmmmin2:iMmmmax2,rec_from:rec_to) - 4.0d0 * w(iMmmin1:iMmmax1,&
       iMmmin2:iMmmax2,rec_from:rec_to) + 3.0d0 * w(iMmin1:iMmax1,&
       iMmin2:iMmax2,rec_from:rec_to))**2
    beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
       2) = beta_coeff(1) * (w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,&
       rec_from:rec_to) + w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
       rec_from:rec_to) - 2.0d0 * w(iMmin1:iMmax1,iMmin2:iMmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iMmmin1:iMmmax1,&
       iMmmin2:iMmmax2,rec_from:rec_to) - w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
       rec_from:rec_to))**2
    beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
       3) = beta_coeff(1) * (w(iMmin1:iMmax1,iMmin2:iMmax2,&
       rec_from:rec_to) + w(iMppmin1:iMppmax1,iMppmin2:iMppmax2,&
       rec_from:rec_to) - 2.0d0 * w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (3.0d0 * w(iMmin1:iMmax1,&
       iMmin2:iMmax2, rec_from:rec_to) - 4.0d0 * w(iMpmin1:iMpmax1,&
       iMpmin2:iMpmax2,rec_from:rec_to) + w(iMppmin1:iMppmax1,&
       iMppmin2:iMppmax2,rec_from:rec_to))**2

    gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to) = 0.0d0 
    do i = 1,3
      gamma_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) = d_array(i) / (beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) + weno_eps_machine)**2
      gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,&
         rec_from:rec_to) = gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,&
         rec_from:rec_to) + gamma_array(iMmin1:iMmax1,iMmin2:iMmax2,&
         rec_from:rec_to,i)
    end do
    do i = 1,3
      kai_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) = gamma_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) / gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to)
      where(kai_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) .lt. exeno_ct) 
        delta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,i) = 0
      elsewhere
        delta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,i) = 1
      endwhere
    end do

    delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) = delta(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to,1) * 1 + delta(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to,3) * 2 + delta(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to,1) * 4 + delta(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to,2)  * 8 + delta(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to,3) * 16

    !> f3
    where(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 31)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- 3.d0 * &
         w(iLmmmmin1:iLmmmmax1,iLmmmmin2:iLmmmmax2,&
         rec_from:rec_to) + 25.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         rec_from:rec_to) - 101.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to) + 319.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) + 214.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - 38.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 4.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         rec_from:rec_to)) / 420.d0
    !> f4
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 30)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (w(iLmmmin1:iLmmmax1,&
         iLmmmin2:iLmmmax2,rec_from:rec_to) - 8.d0 * w(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,rec_from:rec_to) + 37.d0 * w(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + 37.d0 * w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to) - 8.d0 * w(iLppmin1:iLppmax1,&
         iLppmin2:iLppmax2,rec_from:rec_to) + w(iLpppmin1:iLpppmax1,&
         iLpppmin2:iLpppmax2,rec_from:rec_to)) / 60.d0
    !> f5
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 29)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- &
         w(iLmmmmin1:iLmmmmax1,iLmmmmin2:iLmmmmax2,&
         rec_from:rec_to) + 7.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         rec_from:rec_to) - 23.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to) + 57.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) + 22.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - 2.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to)) / 60.d0
    !> f6
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 28)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (2.d0 * &
         w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         rec_from:rec_to) - 13.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to) + 47.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) + 27.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - 3.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to)) / 60.d0
    !> f7
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 24)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- w(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,rec_from:rec_to) + 7.d0 * w(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + 7.d0 * w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to) - w(iLppmin1:iLppmax1,&
         iLppmin2:iLppmax2,rec_from:rec_to)) / 12.d0
    !> f9
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 16)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (2.d0 * &
         w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) + 5.d0 * &
         w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to)) / 6.d0
    !> f8
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 12)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (w(iLmmmin1:iLmmmax1,&
         iLmmmin2:iLmmmax2,rec_from:rec_to) - 5.d0 * w(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,rec_from:rec_to) + 13.d0 * w(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + 3.d0 * w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to)) / 12.d0
    !> f10
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 8)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- w(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,rec_from:rec_to) + 5.d0 * w(iLmin1:iLmax1,&
         iLmin2:iLmax2,rec_from:rec_to) + 2.d0 * w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,rec_from:rec_to)) / 6.d0
    !> f11
    elsewhere
     flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (2.d0 * &
        w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
        rec_from:rec_to) - 7.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
        rec_from:rec_to) + 11.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
        rec_from:rec_to)) / 6.d0
    endwhere

    wLC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)

    !> right side
    beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
       1) = beta_coeff(1) * (w(iMpppmin1:iMpppmax1,iMpppmin2:iMpppmax2,&
       rec_from:rec_to) + w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
       rec_from:rec_to) - 2.0d0 * w(iMppmin1:iMppmax1,iMppmin2:iMppmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iMpppmin1:iMpppmax1,&
       iMpppmin2:iMpppmax2,rec_from:rec_to) - 4.0d0 * w(iMppmin1:iMppmax1,&
       iMppmin2:iMppmax2,rec_from:rec_to) + 3.0d0 * w(iMpmin1:iMpmax1,&
       iMpmin2:iMpmax2,rec_from:rec_to))**2
    beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
       2) = beta_coeff(1) * (w(iMppmin1:iMppmax1,iMppmin2:iMppmax2,&
       rec_from:rec_to) + w(iMmin1:iMmax1,iMmin2:iMmax2,&
       rec_from:rec_to) - 2.0d0 * w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (w(iMppmin1:iMppmax1,&
       iMppmin2:iMppmax2,rec_from:rec_to) - w(iMmin1:iMmax1,iMmin2:iMmax2,&
       rec_from:rec_to))**2
    beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
       3) = beta_coeff(1) * (w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
       rec_from:rec_to) + w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,&
       rec_from:rec_to) - 2.0d0 * w(iMmin1:iMmax1,iMmin2:iMmax2,&
       rec_from:rec_to))**2 + beta_coeff(2) * (3.0d0 * w(iMpmin1:iMpmax1,&
       iMpmin2:iMpmax2, rec_from:rec_to) - 4.0d0 * w(iMmin1:iMmax1,&
       iMmin2:iMmax2,rec_from:rec_to) + w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,&
       rec_from:rec_to))**2

    gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to) = 0.0d0 
    do i = 1,3
      gamma_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) = d_array(i) / (beta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) + weno_eps_machine)**2
      gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,&
         rec_from:rec_to) = gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,&
         rec_from:rec_to) + gamma_array(iMmin1:iMmax1,iMmin2:iMmax2,&
         rec_from:rec_to,i)
    end do
    do i = 1,3
      kai_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) = gamma_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) / gamma_sum(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to)
      where(kai_array(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,&
         i) .lt. exeno_ct) 
        delta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,i) = 0
      elsewhere
        delta(iMmin1:iMmax1,iMmin2:iMmax2,rec_from:rec_to,i) = 1
      endwhere
    end do
 
    delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) = delta(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to,1) * 1 + delta(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to,3) * 2 + delta(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to,1) * 4 + delta(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to,2)  * 8 + delta(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to,3) * 16

    where(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 31)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- 3.d0 * &
         w(iLppppmin1:iLppppmax1,iLppppmin2:iLppppmax2,&
         rec_from:rec_to) + 25.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         rec_from:rec_to) - 101.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 319.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 214.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - 38.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to) + 4.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         rec_from:rec_to)) / 420.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 30)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
         (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         rec_from:rec_to) - 8.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 37.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 37.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - 8.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to) + w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         rec_from:rec_to)) / 60.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 29)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- &
         w(iLppppmin1:iLppppmax1,iLppppmin2:iLppppmax2,&
         rec_from:rec_to) + 7.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         rec_from:rec_to) - 23.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 57.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 22.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - 2.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to)) / 60.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .eq. 28)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (2.d0 * &
         w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         rec_from:rec_to) - 13.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 47.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 27.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - 3.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to)) / 60.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 24)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- &
         w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 7.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 7.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to)) / 12.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 16)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (2.d0 * &
         w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 5.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         rec_from:rec_to)) / 6.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 12)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
         (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         rec_from:rec_to) - 5.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 13.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 3.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to)) / 12.d0
    elsewhere(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) .ge. 8)
      flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (- &
         w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         rec_from:rec_to) + 5.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         rec_from:rec_to) + 2.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         rec_from:rec_to)) / 6.d0
    elsewhere
     flux(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (2.d0 * &
        w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
        rec_from:rec_to) - 7.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
        rec_from:rec_to) + 11.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
        rec_from:rec_to)) / 6.d0
    endwhere

    wRC(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)

  end subroutine exENO7limiter

  subroutine minmod(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2),&
        b(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: minm(ixImin1:ixImax1,ixImin2:ixImax2)

    minm(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (sign(one,a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))+sign(one,b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)))/2.0d0 * min(abs(a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
       abs(b(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

  end subroutine minmod

  subroutine median(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2),&
        b(ixImin1:ixImax1,ixImin2:ixImax2), c(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: med(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision             :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = c(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    med(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + (sign(one,tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))+sign(one,tmp2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)))/2.0d0 * min(abs(tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)),abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

  end subroutine median
end module mod_weno
