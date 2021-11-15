!> Module containing the MP5 (fifth order) flux scheme
module mod_mp5

  implicit none
  private

  public :: MP5limiter

contains

  !> MP5 limiter from Suresh & Huynh 1997 Following the convention of Mignone et
  !> al. 2010. Needs at least three ghost cells
  subroutine MP5limiter(rec_from,rec_to,ixImin1,ixImax1,iLmin1,iLmax1,idims,w,&
     wLC,wRC)

    use mod_global_parameters
    use mod_physics, only: phys_check_prim

    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nprim)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,1:nprim),&
       wLC(ixImin1:ixImax1,1:nprim) 
    ! .. local ..
    integer                         :: iLmmin1,iLmmax1, iLmmmin1,iLmmmax1,&
        iLpmin1,iLpmax1, iLppmin1,iLppmax1, iLpppmin1,iLpppmax1
    integer                         :: idmin1,idmax1, idpmin1,idpmax1,&
        idppmin1,idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1,&
        iepmin1,iepmax1, ieppmin1,ieppmax1
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,1:nprim)  :: f, fmp, fmin,&
        fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1,1:nprim)  :: wRCtmp, wLCtmp
    double precision, dimension(ixImin1:ixImax1) :: tmp, tmp2, tmp3, a, b, c
    integer                         :: flagL(ixImin1:ixImax1),&
        flagR(ixImin1:ixImax1)
    double precision, parameter     :: eps=1.0d-20, alpha=4.0d0
    !double precision                :: alpha
    !----------------------------------------------------------------------------

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    ! range to process:
    !iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

    !{#IFDEF HALL
    ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
    ! also add one ghost zone!
    !   {iL^L=iL^L^LADD1;}
    !}

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

    iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmax1=iLmmax1-kr(idims,1);
    iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
    iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

    f(iLmin1:iLmax1,rec_from:rec_to) = 1.0d0/60.0d0 * (2.0d0* &
       w(iLmmmin1:iLmmmax1,rec_from:rec_to) - 13.0d0* w(iLmmin1:iLmmax1,&
       rec_from:rec_to) + 47.0d0* w(iLmin1:iLmax1,&
       rec_from:rec_to) + 27.0d0* w(iLpmin1:iLpmax1,&
       rec_from:rec_to) - 3.0d0*  w(iLppmin1:iLppmax1,rec_from:rec_to))

    ! get fmp and ful:
    do iw=rec_from, rec_to
       a(iLmin1:iLmax1) = w(iLpmin1:iLpmax1,iw)-w(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = alpha*(w(iLmin1:iLmax1,iw)-w(iLmmin1:iLmmax1,iw))
       call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
       fmp(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + tmp(iLmin1:iLmax1)
       ful(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + b(iLmin1:iLmax1)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1; idmin1=iLmin1-kr(idims,1);
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1+kr(idims,1); iemin1=idmin1;
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);

    d(iemin1:iemax1,rec_from:rec_to) = w(iepmin1:iepmax1,&
       rec_from:rec_to)-2.0d0*w(iemin1:iemax1,&
       rec_from:rec_to)+w(iemmin1:iemmax1,rec_from:rec_to)

    do iw=rec_from, rec_to
       a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idpmin1:idpmax1,iw)
       b(idmin1:idmax1) = 4.0d0*d(idpmin1:idpmax1,iw)-d(idmin1:idmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
       a(idmin1:idmax1) = d(idmin1:idmax1,iw)
       b(idmin1:idmax1) = d(idpmin1:idpmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,rec_from:rec_to) = (w(iLmin1:iLmax1,&
       rec_from:rec_to)+w(iLpmin1:iLpmax1,&
       rec_from:rec_to))/2.0d0-dm4(iLmin1:iLmax1,rec_from:rec_to)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,rec_from:rec_to) = half*(3.0d0*w(iLmin1:iLmax1,&
       rec_from:rec_to) - w(iLmmin1:iLmmax1,&
       rec_from:rec_to)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,rec_from:rec_to)

    fmin(iLmin1:iLmax1,rec_from:rec_to) = max(min(w(iLmin1:iLmax1,&
       rec_from:rec_to),w(iLpmin1:iLpmax1,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       rec_from:rec_to)),min(w(iLmin1:iLmax1,rec_from:rec_to),&
       ful(iLmin1:iLmax1,rec_from:rec_to),flc(iLmin1:iLmax1,rec_from:rec_to)))

    fmax(iLmin1:iLmax1,rec_from:rec_to) = min(max(w(iLmin1:iLmax1,&
       rec_from:rec_to),w(iLpmin1:iLpmax1,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       rec_from:rec_to)),max(w(iLmin1:iLmax1,rec_from:rec_to),&
       ful(iLmin1:iLmax1,rec_from:rec_to),flc(iLmin1:iLmax1,rec_from:rec_to)))

    do iw=rec_from, rec_to
       a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
       c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
       call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
       flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,rec_from:rec_to)-w(iLmin1:iLmax1,&
       rec_from:rec_to))*(f(iLmin1:iLmax1,rec_from:rec_to)-fmp(iLmin1:iLmax1,&
       rec_from:rec_to)) .le. eps)
       wLCtmp(iLmin1:iLmax1,rec_from:rec_to) = f(iLmin1:iLmax1,&
          rec_from:rec_to)
    elsewhere
       wLCtmp(iLmin1:iLmax1,rec_from:rec_to) = flim(iLmin1:iLmax1,&
          rec_from:rec_to)
    end where

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand fromula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLpppmin1=iLppmin1+kr(idims,1);iLpppmax1=iLppmax1+kr(idims,1);

    f(iLmin1:iLmax1,rec_from:rec_to) = 1.0d0/60.0d0 * (2.0d0* &
       w(iLpppmin1:iLpppmax1,rec_from:rec_to) - 13.0d0* w(iLppmin1:iLppmax1,&
       rec_from:rec_to) + 47.0d0* w(iLpmin1:iLpmax1,&
       rec_from:rec_to) + 27.0d0* w(iLmin1:iLmax1,&
       rec_from:rec_to) - 3.0d0*  w(iLmmin1:iLmmax1,rec_from:rec_to))

    ! get fmp and ful:
    do iw=rec_from, rec_to
       a(iLmin1:iLmax1) = w(iLmin1:iLmax1,iw)-w(iLpmin1:iLpmax1,iw)
       b(iLmin1:iLmax1) = alpha*(w(iLpmin1:iLpmax1,iw)-w(iLppmin1:iLppmax1,&
          iw))
       call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
       fmp(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + tmp(iLmin1:iLmax1)
       ful(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + b(iLmin1:iLmax1)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1+kr(idims,1); idmin1=iLmin1;
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1; iemin1=idmin1-kr(idims,1);
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);
    ieppmin1=iepmin1+kr(idims,1);ieppmax1=iepmax1+kr(idims,1);

    d(iemin1:iemax1,rec_from:rec_to) = w(iemin1:iemax1,&
       rec_from:rec_to)-2.0d0*w(iepmin1:iepmax1,&
       rec_from:rec_to)+w(ieppmin1:ieppmax1,rec_from:rec_to)

    do iw=rec_from, rec_to
       a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idmmin1:idmmax1,iw)
       b(idmin1:idmax1) = 4.0d0*d(idmmin1:idmmax1,iw)-d(idmin1:idmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
       a(idmin1:idmax1) = d(idmin1:idmax1,iw)
       b(idmin1:idmax1) = d(idmmin1:idmmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,rec_from:rec_to) = (w(iLmin1:iLmax1,&
       rec_from:rec_to)+w(iLpmin1:iLpmax1,&
       rec_from:rec_to))/2.0d0-dm4(iLmin1:iLmax1,rec_from:rec_to)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,rec_from:rec_to) = half*(3.0d0*w(iLpmin1:iLpmax1,&
       rec_from:rec_to) - w(iLppmin1:iLppmax1,&
       rec_from:rec_to)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,rec_from:rec_to)

    fmin(iLmin1:iLmax1,rec_from:rec_to) = max(min(w(iLpmin1:iLpmax1,&
       rec_from:rec_to),w(iLmin1:iLmax1,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       rec_from:rec_to)),min(w(iLpmin1:iLpmax1,rec_from:rec_to),&
       ful(iLmin1:iLmax1,rec_from:rec_to),flc(iLmin1:iLmax1,rec_from:rec_to)))

    fmax(iLmin1:iLmax1,rec_from:rec_to) = min(max(w(iLpmin1:iLpmax1,&
       rec_from:rec_to),w(iLmin1:iLmax1,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       rec_from:rec_to)),max(w(iLpmin1:iLpmax1,rec_from:rec_to),&
       ful(iLmin1:iLmax1,rec_from:rec_to),flc(iLmin1:iLmax1,rec_from:rec_to)))

    do iw=rec_from, rec_to
       a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
       c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
       call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
       flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,rec_from:rec_to)-w(iLpmin1:iLpmax1,&
       rec_from:rec_to))*(f(iLmin1:iLmax1,rec_from:rec_to)-fmp(iLmin1:iLmax1,&
       rec_from:rec_to))  .le. eps)
       wRCtmp(iLmin1:iLmax1,rec_from:rec_to) = f(iLmin1:iLmax1,&
          rec_from:rec_to)
    elsewhere
       wRCtmp(iLmin1:iLmax1,rec_from:rec_to) = flim(iLmin1:iLmax1,&
          rec_from:rec_to)
    end where

    ! Since limiter not TVD, negative pressures or densities could result.  
    ! Fall back to flat interpolation (minmod would also work). 
    call phys_check_prim(ixGlo1,ixGhi1,iLmin1,iLmax1,wLCtmp,flagL)
    call phys_check_prim(ixGlo1,ixGhi1,iLmin1,iLmax1,wRCtmp,flagR)

    do iw=rec_from, rec_to
       where (flagL(iLmin1:iLmax1) == 0 .and. flagR(iLmin1:iLmax1) == 0)
          wLC(iLmin1:iLmax1,iw)=wLCtmp(iLmin1:iLmax1,iw)
          wRC(iLmin1:iLmax1,iw)=wRCtmp(iLmin1:iLmax1,iw)
       end where
    end do

  end subroutine MP5limiter

  subroutine minmod(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: a(ixImin1:ixImax1), b(ixImin1:ixImax1)
    double precision, intent(out):: minm(ixImin1:ixImax1)

    minm(ixOmin1:ixOmax1) = (sign(one,a(ixOmin1:ixOmax1))+sign(one,&
       b(ixOmin1:ixOmax1)))/2.0d0 * min(abs(a(ixOmin1:ixOmax1)),&
       abs(b(ixOmin1:ixOmax1)))

  end subroutine minmod

  subroutine median(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: a(ixImin1:ixImax1), b(ixImin1:ixImax1),&
        c(ixImin1:ixImax1)
    double precision, intent(out):: med(ixImin1:ixImax1)

    double precision             :: tmp1(ixImin1:ixImax1),&
       tmp2(ixImin1:ixImax1)

    tmp1(ixOmin1:ixOmax1) = b(ixOmin1:ixOmax1) - a(ixOmin1:ixOmax1)
    tmp2(ixOmin1:ixOmax1) = c(ixOmin1:ixOmax1) - a(ixOmin1:ixOmax1)

    med(ixOmin1:ixOmax1) = a(ixOmin1:ixOmax1) + (sign(one,&
       tmp1(ixOmin1:ixOmax1))+sign(one,tmp2(ixOmin1:ixOmax1)))/2.0d0 * &
       min(abs(tmp1(ixOmin1:ixOmax1)),abs(tmp2(ixOmin1:ixOmax1)))

  end subroutine median

end module mod_mp5
