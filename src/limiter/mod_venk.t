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

  subroutine venklimiter(rec_from,rec_to,ixI^L,iL^L,idims,dxdim,w,wLp,wRp)
    use mod_global_parameters

    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nprim)
    double precision, intent(inout) :: wRp(ixI^S,1:nprim),wLp(ixI^S,1:nprim) 
    !> local
    integer                         :: iM^L, iMm^L, iMp^L
    integer                         :: iLm^L, iLp^L, iLpp^L
    double precision                :: wmax(ixI^S,1:nprim),wmin(ixI^S,1:nprim)
    double precision                :: westp(ixI^S,1:nprim),westm(ixI^S,1:nprim)
    double precision                :: phi1(ixI^S,1:nprim),phi2(ixI^S,1:nprim)
    double precision                :: phi3(ixI^S,1:nprim),phi4(ixI^S,1:nprim)
    double precision                :: phim(ixI^S,1:nprim),phip(ixI^S,1:nprim),phi(ixI^S,1:nprim)
    double precision                :: deltap(ixI^S,1:nprim),deltam(ixI^S,1:nprim)
    double precision                :: eps2
    double precision, parameter     :: venk_omega = 1.d-12
    double precision, parameter     :: venk_k = 0.3d0

    iM^L=iL^L;
    iMmax^D=iLmax^D+kr(idims,^D);
    iLm^L=iL^L-kr(idims,^D);
    iMm^L=iM^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iMp^L=iM^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);

    eps2 = (venk_k * dxdim) ** 3

    wmax(iM^S,rec_from:rec_to) = max(w(iMm^S,rec_from:rec_to), w(iM^S,rec_from:rec_to), w(iMp^S,rec_from:rec_to))
    wmin(iM^S,rec_from:rec_to) = min(w(iMm^S,rec_from:rec_to), w(iM^S,rec_from:rec_to), w(iMp^S,rec_from:rec_to))
    !> use central difference approximation as (eq.5) and take phi = 1
    westp(iM^S,rec_from:rec_to) = w(iM^S,rec_from:rec_to) + (w(iMp^S,rec_from:rec_to) - w(iMm^S,rec_from:rec_to)) * 0.25d0
    westm(iM^S,rec_from:rec_to) = w(iM^S,rec_from:rec_to) - (w(iMp^S,rec_from:rec_to) - w(iMm^S,rec_from:rec_to)) * 0.25d0

    !> (eq.30) & (eq.31)
    deltap = 0
    deltam = 0
    phi1 = 0
    phi2 = 0
    phip = 1
    where(westp(iM^S,rec_from:rec_to) .gt. w(iM^S,rec_from:rec_to))
      deltap(iM^S,rec_from:rec_to) = wmax(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to)
      deltam(iM^S,rec_from:rec_to) = westp(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to)
      deltam(iM^S,rec_from:rec_to) = sign(dabs(deltam(iM^S,rec_from:rec_to)) + venk_omega, deltam(iM^S,rec_from:rec_to))
      phi1(iM^S,rec_from:rec_to) = (deltap(iM^S,rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to)
      phi2(iM^S,rec_from:rec_to) = deltap(iM^S,rec_from:rec_to)**2 + 2.d0 * deltam(iM^S,rec_from:rec_to)**2 + deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to) + eps2
      phip(iM^S,rec_from:rec_to) = phi1(iM^S,rec_from:rec_to) / phi2(iM^S, rec_from:rec_to)
    elsewhere(westp(iM^S,rec_from:rec_to) .lt. w(iM^S,rec_from:rec_to))
      deltap(iM^S,rec_from:rec_to) = wmin(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to)
      deltam(iM^S,rec_from:rec_to) = westp(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to)
      deltam(iM^S,rec_from:rec_to) = sign(dabs(deltam(iM^S,rec_from:rec_to)) + venk_omega, deltam(iM^S,rec_from:rec_to))
      phi1(iM^S,rec_from:rec_to) = (deltap(iM^S,rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to)
      phi2(iM^S,rec_from:rec_to) = deltap(iM^S,rec_from:rec_to)**2 + 2.d0 * deltam(iM^S,rec_from:rec_to)**2 + deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to) + eps2
      phip(iM^S,rec_from:rec_to) = phi1(iM^S,rec_from:rec_to) / phi2(iM^S, rec_from:rec_to)
    elsewhere
      phip(iM^S,rec_from:rec_to) = one
   endwhere

    deltap = 0
    deltam = 0
    phi3 = 0
    phi4 = 0
    phim = 0
    where(westm(iM^S,rec_from:rec_to) .lt. w(iM^S,rec_from:rec_to))
      deltap(iM^S,rec_from:rec_to) = - (wmax(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to))
      deltam(iM^S,rec_from:rec_to) = westm(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to)
      deltam(iM^S,rec_from:rec_to) = sign(dabs(deltam(iM^S,rec_from:rec_to)) + venk_omega, deltam(iM^S,rec_from:rec_to))
      phi3(iM^S,rec_from:rec_to) = (deltap(iM^S,rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to)
      phi4(iM^S,rec_from:rec_to) = deltap(iM^S,rec_from:rec_to)**2 + 2.d0 * deltam(iM^S,rec_from:rec_to)**2 + deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to) + eps2
      phim(iM^S,rec_from:rec_to) = phi3(iM^S,rec_from:rec_to) / phi4(iM^S, rec_from:rec_to)
    elsewhere(westm(iM^S,rec_from:rec_to) .gt. w(iM^S,rec_from:rec_to))
      deltap(iM^S,rec_from:rec_to) = - (wmin(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to))
      deltam(iM^S,rec_from:rec_to) = westm(iM^S,rec_from:rec_to) - w(iM^S,rec_from:rec_to)
      deltam(iM^S,rec_from:rec_to) = sign(dabs(deltam(iM^S,rec_from:rec_to)) + venk_omega, deltam(iM^S,rec_from:rec_to))
      phi3(iM^S,rec_from:rec_to) = (deltap(iM^S,rec_from:rec_to)**2 + eps2) + 2.d0 * deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to)
      phi4(iM^S,rec_from:rec_to) = deltap(iM^S,rec_from:rec_to)**2 + 2.d0 * deltam(iM^S,rec_from:rec_to)**2 + deltap(iM^S,rec_from:rec_to) * deltam(iM^S,rec_from:rec_to) + eps2
      phim(iM^S,rec_from:rec_to) = phi3(iM^S,rec_from:rec_to) / phi4(iM^S, rec_from:rec_to)
    elsewhere
      phim(iM^S,rec_from:rec_to) = one
    endwhere
    !> (eq.3)
    phi(iM^S,rec_from:rec_to) = min(phim(iM^S,rec_from:rec_to),phip(iM^S,rec_from:rec_to))
    !> (eq.5)
    wLp(iL^S,rec_from:rec_to) = w(iL^S,rec_from:rec_to) + 0.25d0 * (w(iLp^S,rec_from:rec_to)-w(iLm^S,rec_from:rec_to)) * phi(iL^S,rec_from:rec_to)
    wRp(iL^S,rec_from:rec_to) = w(iLp^S,rec_from:rec_to) - 0.25d0 * (w(iLpp^S,rec_from:rec_to)-w(iL^S,rec_from:rec_to)) * phi(iLp^S,rec_from:rec_to)

  end subroutine venklimiter

end module mod_venk
