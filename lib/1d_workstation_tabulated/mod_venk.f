module mod_venk
  ! Venkatakrishnan limiter
  !
  ! 2019.10.11 coded up by nanami;
  !
  ! see Venkatakrishnan 1993 for this limiter;


  implicit none
  private

  public :: venklimiter

contains

  subroutine venklimiter(rec_from,rec_to,ixImin1,ixImax1,iLmin1,iLmax1,idims,&
     dxdim,w,wLp,wRp)
    use mod_global_parameters

    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nprim)
    double precision, intent(inout) :: wRp(ixImin1:ixImax1,1:nprim),&
       wLp(ixImin1:ixImax1,1:nprim) 
    !> local
    integer                         :: iMmin1,iMmax1, iMmmin1,iMmmax1, iMpmin1,&
       iMpmax1
    integer                         :: iLmmin1,iLmmax1, iLpmin1,iLpmax1,&
        iLppmin1,iLppmax1
    double precision                :: wmax(ixImin1:ixImax1,1:nprim),&
       wmin(ixImin1:ixImax1,1:nprim)
    double precision                :: westp(ixImin1:ixImax1,1:nprim),&
       westm(ixImin1:ixImax1,1:nprim)
    double precision                :: phi1(ixImin1:ixImax1,1:nprim),&
       phi2(ixImin1:ixImax1,1:nprim)
    double precision                :: phi3(ixImin1:ixImax1,1:nprim),&
       phi4(ixImin1:ixImax1,1:nprim)
    double precision                :: phim(ixImin1:ixImax1,1:nprim),&
       phip(ixImin1:ixImax1,1:nprim),phi(ixImin1:ixImax1,1:nprim)
    double precision                :: deltap(ixImin1:ixImax1,1:nprim),&
       deltam(ixImin1:ixImax1,1:nprim)
    double precision                :: eps2
    double precision, parameter     :: venk_omega = 1.d-12
    double precision, parameter     :: venk_k = 0.3d0

    iMmin1=iLmin1;iMmax1=iLmax1;
    iMmax1=iLmax1+kr(idims,1);
    iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
    iMmmin1=iMmin1-kr(idims,1);iMmmax1=iMmax1-kr(idims,1);
    iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
    iMpmin1=iMmin1+kr(idims,1);iMpmax1=iMmax1+kr(idims,1);
    iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

    eps2 = (venk_k * dxdim) ** 3

    wmax(iMmin1:iMmax1,rec_from:rec_to) = max(w(iMmmin1:iMmmax1,&
       rec_from:rec_to), w(iMmin1:iMmax1,rec_from:rec_to), w(iMpmin1:iMpmax1,&
       rec_from:rec_to))
    wmin(iMmin1:iMmax1,rec_from:rec_to) = min(w(iMmmin1:iMmmax1,&
       rec_from:rec_to), w(iMmin1:iMmax1,rec_from:rec_to), w(iMpmin1:iMpmax1,&
       rec_from:rec_to))
    !> use central difference approximation as (eq.5) and take phi = 1
    westp(iMmin1:iMmax1,rec_from:rec_to) = w(iMmin1:iMmax1,&
       rec_from:rec_to) + (w(iMpmin1:iMpmax1,&
       rec_from:rec_to) - w(iMmmin1:iMmmax1,rec_from:rec_to)) * 0.25d0
    westm(iMmin1:iMmax1,rec_from:rec_to) = w(iMmin1:iMmax1,&
       rec_from:rec_to) - (w(iMpmin1:iMpmax1,&
       rec_from:rec_to) - w(iMmmin1:iMmmax1,rec_from:rec_to)) * 0.25d0

    !> (eq.30) & (eq.31)
    deltap = 0
    deltam = 0
    phi1 = 0
    phi2 = 0
    phip = 1
    where(westp(iMmin1:iMmax1,rec_from:rec_to) .gt. w(iMmin1:iMmax1,&
       rec_from:rec_to))
      deltap(iMmin1:iMmax1,rec_from:rec_to) = wmax(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to)
      deltam(iMmin1:iMmax1,rec_from:rec_to) = westp(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to)
      deltam(iMmin1:iMmax1,rec_from:rec_to) = sign(dabs(deltam(iMmin1:iMmax1,&
         rec_from:rec_to)) + venk_omega, deltam(iMmin1:iMmax1,&
         rec_from:rec_to))
      phi1(iMmin1:iMmax1,rec_from:rec_to) = (deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to)
      phi2(iMmin1:iMmax1,rec_from:rec_to) = deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + 2.d0 * deltam(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to) + eps2
      phip(iMmin1:iMmax1,rec_from:rec_to) = phi1(iMmin1:iMmax1,&
         rec_from:rec_to) / phi2(iMmin1:iMmax1, rec_from:rec_to)
    elsewhere(westp(iMmin1:iMmax1,rec_from:rec_to) .lt. w(iMmin1:iMmax1,&
       rec_from:rec_to))
      deltap(iMmin1:iMmax1,rec_from:rec_to) = wmin(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to)
      deltam(iMmin1:iMmax1,rec_from:rec_to) = westp(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to)
      deltam(iMmin1:iMmax1,rec_from:rec_to) = sign(dabs(deltam(iMmin1:iMmax1,&
         rec_from:rec_to)) + venk_omega, deltam(iMmin1:iMmax1,&
         rec_from:rec_to))
      phi1(iMmin1:iMmax1,rec_from:rec_to) = (deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to)
      phi2(iMmin1:iMmax1,rec_from:rec_to) = deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + 2.d0 * deltam(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to) + eps2
      phip(iMmin1:iMmax1,rec_from:rec_to) = phi1(iMmin1:iMmax1,&
         rec_from:rec_to) / phi2(iMmin1:iMmax1, rec_from:rec_to)
    elsewhere
      phip(iMmin1:iMmax1,rec_from:rec_to) = one
   endwhere

    deltap = 0
    deltam = 0
    phi3 = 0
    phi4 = 0
    phim = 0
    where(westm(iMmin1:iMmax1,rec_from:rec_to) .lt. w(iMmin1:iMmax1,&
       rec_from:rec_to))
      deltap(iMmin1:iMmax1,rec_from:rec_to) = - (wmax(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to))
      deltam(iMmin1:iMmax1,rec_from:rec_to) = westm(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to)
      deltam(iMmin1:iMmax1,rec_from:rec_to) = sign(dabs(deltam(iMmin1:iMmax1,&
         rec_from:rec_to)) + venk_omega, deltam(iMmin1:iMmax1,&
         rec_from:rec_to))
      phi3(iMmin1:iMmax1,rec_from:rec_to) = (deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to)
      phi4(iMmin1:iMmax1,rec_from:rec_to) = deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + 2.d0 * deltam(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to) + eps2
      phim(iMmin1:iMmax1,rec_from:rec_to) = phi3(iMmin1:iMmax1,&
         rec_from:rec_to) / phi4(iMmin1:iMmax1, rec_from:rec_to)
    elsewhere(westm(iMmin1:iMmax1,rec_from:rec_to) .gt. w(iMmin1:iMmax1,&
       rec_from:rec_to))
      deltap(iMmin1:iMmax1,rec_from:rec_to) = - (wmin(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to))
      deltam(iMmin1:iMmax1,rec_from:rec_to) = westm(iMmin1:iMmax1,&
         rec_from:rec_to) - w(iMmin1:iMmax1,rec_from:rec_to)
      deltam(iMmin1:iMmax1,rec_from:rec_to) = sign(dabs(deltam(iMmin1:iMmax1,&
         rec_from:rec_to)) + venk_omega, deltam(iMmin1:iMmax1,&
         rec_from:rec_to))
      phi3(iMmin1:iMmax1,rec_from:rec_to) = (deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to)
      phi4(iMmin1:iMmax1,rec_from:rec_to) = deltap(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + 2.d0 * deltam(iMmin1:iMmax1,&
         rec_from:rec_to)**2 + deltap(iMmin1:iMmax1,&
         rec_from:rec_to) * deltam(iMmin1:iMmax1,rec_from:rec_to) + eps2
      phim(iMmin1:iMmax1,rec_from:rec_to) = phi3(iMmin1:iMmax1,&
         rec_from:rec_to) / phi4(iMmin1:iMmax1, rec_from:rec_to)
    elsewhere
      phim(iMmin1:iMmax1,rec_from:rec_to) = one
    endwhere
    !> (eq.3)
    phi(iMmin1:iMmax1,rec_from:rec_to) = min(phim(iMmin1:iMmax1,&
       rec_from:rec_to),phip(iMmin1:iMmax1,rec_from:rec_to))
    !> (eq.5)
    wLp(iLmin1:iLmax1,rec_from:rec_to) = w(iLmin1:iLmax1,&
       rec_from:rec_to) + 0.25d0 * (w(iLpmin1:iLpmax1,&
       rec_from:rec_to)-w(iLmmin1:iLmmax1,&
       rec_from:rec_to)) * phi(iLmin1:iLmax1,rec_from:rec_to)
    wRp(iLmin1:iLmax1,rec_from:rec_to) = w(iLpmin1:iLpmax1,&
       rec_from:rec_to) - 0.25d0 * (w(iLppmin1:iLppmax1,&
       rec_from:rec_to)-w(iLmin1:iLmax1,rec_from:rec_to)) * &
       phi(iLpmin1:iLpmax1,rec_from:rec_to)

  end subroutine venklimiter

end module mod_venk
