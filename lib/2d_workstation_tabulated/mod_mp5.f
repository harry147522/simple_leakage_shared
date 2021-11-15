!> Module containing the MP5 (fifth order) flux scheme
module mod_mp5

  implicit none
  private

  public :: MP5limiter

contains

  !> MP5 limiter from Suresh & Huynh 1997 Following the convention of Mignone et
  !> al. 2010. Needs at least three ghost cells
  subroutine MP5limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,&
     iLmin2,iLmax1,iLmax2,idims,w,wLC,wRC)

    use mod_global_parameters
    use mod_physics, only: phys_check_prim

    integer, intent(in)             :: rec_from, rec_to
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim) 
    ! .. local ..
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
       idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
       idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
       iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
       ieppmax1,ieppmax2
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim)  :: f,&
        fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim)  :: wRCtmp, wLCtmp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: tmp, tmp2,&
        tmp3, a, b, c
    integer                         :: flagL(ixImin1:ixImax1,ixImin2:ixImax2),&
        flagR(ixImin1:ixImax1,ixImin2:ixImax2)
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

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);

    f(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 1.0d0/60.0d0 * (2.0d0* &
       w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       rec_from:rec_to) - 13.0d0* w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to) + 47.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to))

    ! get fmp and ful:
    do iw=rec_from, rec_to
       a(iLmin1:iLmax1,iLmin2:iLmax2) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
          iw)-w(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2) = alpha*(w(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iw))
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,&
          a,b,tmp)
       fmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2)
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw) + b(iLmin1:iLmax1,iLmin2:iLmax2)
    end do ! iw loop

    ! get dm4:
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

    do iw=rec_from, rec_to
       a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,&
          iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
       b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idpmin1:idpmax1,&
          idpmin2:idpmax2,iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
          a,b,tmp)
       a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
       b(idmin1:idmax1,idmin2:idmax2) = d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
          a,b,tmp2)
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
          tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,idmin2:idmax2)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
       half*(3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
       max(min(w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
       w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)),min(w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to),ful(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
       min(max(w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
       w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)),max(w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to),ful(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))

    do iw=rec_from, rec_to
       a(iLmin1:iLmax1,iLmin2:iLmax2) = fmin(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2) = f(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       c(iLmin1:iLmax1,iLmin2:iLmax2) = fmax(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       call median(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,&
          a,b,c,tmp)
       flim(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)-w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)) .le. eps)
       wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,rec_from:rec_to)
    elsewhere
       wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) = flim(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)
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

    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);

    f(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = 1.0d0/60.0d0 * (2.0d0* &
       w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       rec_from:rec_to) - 13.0d0* w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to) + 47.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) + 27.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to) - 3.0d0*  w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       rec_from:rec_to))

    ! get fmp and ful:
    do iw=rec_from, rec_to
       a(iLmin1:iLmax1,iLmin2:iLmax2) = w(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2) = alpha*(w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iw)-w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iw))
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,&
          a,b,tmp)
       fmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
          iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2)
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
          iw) + b(iLmin1:iLmax1,iLmin2:iLmax2)
    end do ! iw loop

    ! get dm4:
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
       iemin2:iemax2,rec_from:rec_to)-2.0d0*w(iepmin1:iepmax1,iepmin2:iepmax2,&
       rec_from:rec_to)+w(ieppmin1:ieppmax1,ieppmin2:ieppmax2,rec_from:rec_to)

    do iw=rec_from, rec_to
       a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,&
          iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
       b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmmin1:idmmax1,&
          idmmin2:idmmax2,iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
          a,b,tmp)
       a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
       b(idmin1:idmax1,idmin2:idmax2) = d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
          a,b,tmp2)
       call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
          tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,idmin2:idmax2)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
       half*(3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       rec_from:rec_to)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
       max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),&
       w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)),min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to),ful(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = &
       min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,rec_from:rec_to),&
       w(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),fmd(iLmin1:iLmax1,&
       iLmin2:iLmax2,rec_from:rec_to)),max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       rec_from:rec_to),ful(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)))

    do iw=rec_from, rec_to
       a(iLmin1:iLmax1,iLmin2:iLmax2) = fmin(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2) = f(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       c(iLmin1:iLmax1,iLmin2:iLmax2) = fmax(iLmin1:iLmax1,iLmin2:iLmax2,iw)
       call median(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,&
          a,b,c,tmp)
       flim(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)-w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,rec_from:rec_to))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,&
       rec_from:rec_to))  .le. eps)
       wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,rec_from:rec_to)
    elsewhere
       wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,&
          rec_from:rec_to) = flim(iLmin1:iLmax1,iLmin2:iLmax2,rec_from:rec_to)
    end where

    ! Since limiter not TVD, negative pressures or densities could result.  
    ! Fall back to flat interpolation (minmod would also work). 
    call phys_check_prim(ixGlo1,ixGlo2,ixGhi1,ixGhi2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,flagL)
    call phys_check_prim(ixGlo1,ixGlo2,ixGhi1,ixGhi2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wRCtmp,flagR)

    do iw=rec_from, rec_to
       where (flagL(iLmin1:iLmax1,iLmin2:iLmax2) == 0 .and. &
          flagR(iLmin1:iLmax1,iLmin2:iLmax2) == 0)
          wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLCtmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
          wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRCtmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
       end where
    end do

  end subroutine MP5limiter

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

end module mod_mp5
