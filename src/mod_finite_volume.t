!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume

contains

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,sold,fC,fE,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_cfc

    character(len=*), intent(in)                          :: method
    double precision, intent(in)                          :: qdt, qtC, qt, dx^D
    integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixI^S,1:ncons,1:ndim)     :: fC
    double precision, dimension(ixI^S,7-2*ndim:3)         :: fE

    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:ncons) :: consL, consR
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nprim) :: primL, primR
    ! High order flux
    double precision, dimension(ixI^S,1:ncons) :: fLC, fRC
    ! Low order flux
    double precision, dimension(ixI^S)      :: cmaxC
    double precision, dimension(ixI^S)      :: cminC
    double precision, dimension(ixO^S)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    ! cell-face location coordinates
    double precision, dimension(ixI^S,1:ndim) :: xi
    integer, dimension(ixI^S)               :: patchf
    integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L

    ! variables that related to positivity preserving limiter
    double precision, dimension(ixI^S,1:ncons,1:ndim) :: fClow
    double precision, parameter :: epsD = 1.0d-16, epstau = 1.0d-16

    associate(primCT=>sCT%prim, &
              consCT=>sCT%cons, &
              cons_new=>snew%cons,&
              cons_old=>sold%cons)

    fC=0.d0
    fClow=0.d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idims= idims^LIM
       ! use interface value of w0 at idims
       block%iw0=idims

       hxO^L=ixO^L-kr(idims,^D);

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       if(stagger_grid) then
          ! ct needs all transverse cells
          ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if

       ! Determine stencil size
       {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
       {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

       ! get cell-face coordinates
       xi=x
       xi(ixI^S,idims)=xi(ixI^S,idims)+0.5d0*sCT%dx(ixI^S,idims)

       ! primR and primL are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'higher' direction.
       primR(kxC^S,1:nprim)=primCT(kxR^S,1:nprim)
       primL(kxC^S,1:nprim)=primCT(kxC^S,1:nprim)

       ! low order TVD-Lax-Friedrich flux first.
       if (positivity_preserving) then
         ! fixme: maybe this part can be faster
         call phys_to_conserved(ixI^L,ixCR^L,consR,primR,xi)
         call phys_to_conserved(ixI^L,ixCR^L,consL,primL,xi)
         call phys_get_flux(consL,primL,xi,ixI^L,ixC^L,idims,fLC)
         call phys_get_flux(consR,primR,xi,ixI^L,ixC^L,idims,fRC)
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,cmaxC)
         call get_Riemann_flux_tvdlf()
         if ( slab_uniform ) then
           do iw=nc_hydro_lo,nc_hydro_hi
             fClow(ixI^S,iw,idims)=fC(ixI^S,iw,idims)
           end do
         else 
           do iw=nc_hydro_lo,nc_hydro_hi
             fClow(ixI^S,iw,idims)=fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
           end do
         end if
       end if

       if ( nmetric > 0 ) then
          ! fixme: make it faster, fewer index

          ! get cell-face metric
          if ( reconstruct_cfc ) then
             primR(kxC^S,nmetric_lo:nmetric_hi)=primCT(kxR^S,nmetric_lo:nmetric_hi)
             primL(kxC^S,nmetric_lo:nmetric_hi)=primCT(kxC^S,nmetric_lo:nmetric_hi)
             call reconstruct_LR(limiter_wenozp5, nmetric_lo, nmetric_hi, &
                 ixI^L,ixCR^L,ixCR^L,idims,primCT,consL,consR,primL,primR,xi,dxdim(idims))
          else
             primL(kxC^S,nmetric_lo:nmetric_hi)=primCT(kxC^S,nmetric_lo:nmetric_hi)
             call cfc_metric_interpolation(ixI^L,ixCR^L,idims,primCT,x,primL,xi)
             primR(ixCR^S,nmetric_lo:nmetric_hi)=primL(ixCR^S,nmetric_lo:nmetric_hi)
          end if
       end if
       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(typelimiter, nhydro_lo, nhydro_lo+nreconstruct-1, &
         ixI^L,ixCR^L,ixCR^L,idims,primCT,consL,consR,primL,primR,xi,dxdim(idims))

       ! special modification of left and right status before flux evaluation
       !call phys_modify_wLR(ixI^L,ixCR^L,consL,consR,primL,primR,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(consL,primL,xi,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(consR,primR,xi,ixI^L,ixC^L,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,cmaxC)
       else
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,cmaxC,cminC)
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case('tvdmu')
         call get_Riemann_flux_tvdmu()
       case('tvdlf')
         call get_Riemann_flux_tvdlf()
       case('hll')
         call get_Riemann_flux_hll()
       case('hllc','hllcd')
         call get_Riemann_flux_hllc()
       case default
         call mpistop('unkown Riemann flux')
       end select

       ! special modification of flux
       ! fixme: for M1, dont use it at the moment
       !call phys_modify_flux(ixI^L,ixCR^L,consL,consR,primL,primR,xi,sCT,idims,fC)

       if (.not.slab_uniform) then
          do iw=nc_hydro_lo,nc_hydro_hi
            fC(ixI^S,iw,idims)=fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
          end do
       end if

    end do ! Next idims
    block%iw0=0

    ! If use positivity preserving limiter, with fC and fClow, work out the modify the flux
    if (positivity_preserving) then
       call positivity_preserving_limiter()
    end if

    if(associated(usr_set_flux)) then
       do idims= idims^LIM
          call usr_set_flux(ixI^L,ixC^L,idims,fC,xi)
       end do ! Next idims
    end if

    ! fixme: this might not right as I modified the sturcture
    if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qdt,primCT,fC,fE,sCT,snew)

    do idims= idims^LIM
       hxO^L=ixO^L-kr(idims,^D);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          do iw = nc_hydro_lo,nc_hydro_hi
            if ( flux_type(idims, iw) /= flux_nul ) then
               if (associated(phys_iw_methods(iw)%inv_capacity)) then
                 call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, cons_new, inv_volume)
                 cons_new(ixO^S,iw)=cons_new(ixO^S,iw) + inv_volume * &
                      dxinv(idims)*(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
               else
                 cons_new(ixO^S,iw)=cons_new(ixO^S,iw) &
                      + dxinv(idims)*(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
               end if
            end if
          end do
        else
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, cons_new, inv_volume)
          else
            inv_volume(ixO^S) = 1.0d0
          end if
          inv_volume(ixO^S) = inv_volume(ixO^S)/block%dvolume(ixO^S)

          do iw=nc_hydro_lo,nc_hydro_hi
            if ( flux_type(idims, iw) /= flux_nul ) then
               cons_new(ixO^S,iw)=cons_new(ixO^S,iw) &
                  - qdt * inv_volume(ixO^S) * (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) 
            end if
          enddo
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,consL,consR,cons_new,x,fC,dx^D)

    end do ! Next idims

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,consCT,primCT,cons_new,x)

    if(stagger_grid) call phys_face_to_center(ixO^L,snew)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nprim,qtC,primCT,qt,cons_new,x,.false.)

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()
      do iw=nc_hydro_lo,nc_hydro_hi
         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixC^S)

      fac = -0.5d0*tvdlfeps*cmaxC(ixC^S)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=nc_hydro_lo,nc_hydro_hi

         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + fac*(consR(ixC^S,iw)-consL(ixC^S,iw))
         end if

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixC^S), div(ixC^S)

      where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
      elsewhere
        patchf(ixC^S) =  1
        fac(ixC^S) = tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)
        div(ixC^S) = 1/(cmaxC(ixC^S)-cminC(ixC^S))
      endwhere

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=nc_hydro_lo,nc_hydro_hi
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S, iw) = 0.5d0*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                 (consR(ixC^S,iw)-consL(ixC^S,iw)))
         else
            where(patchf(ixC^S)==1)
               ! Add hll dissipation to the flux
               fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                    +fac(ixC^S)*(consR(ixC^S,iw)-consL(ixC^S,iw))) * div(ixC^S)
            elsewhere(patchf(ixC^S)== 2)
               fLC(ixC^S, iw)=fRC(ixC^S, iw)
            elsewhere(patchf(ixC^S)==-2)
               fLC(ixC^S, iw)=fLC(ixC^S, iw)
            endwhere
         endif

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      implicit none
      double precision, dimension(ixI^S,1:ncons)     :: whll, Fhll, fCD
      double precision, dimension(ixI^S)              :: lambdaCD

      patchf(ixC^S) =  1
      where(cminC(ixC^S) >= zero)
         patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
         patchf(ixC^S) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') &
           call phys_diffuse_hllcd(ixI^L,ixC^L,idims,consL,consR,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixC^S)==1)) &
           call phys_get_lCD(consL,consR,fLC,fRC,cminC,cmaxC,idims,ixI^L,ixC^L, &
           whll,Fhll,lambdaCD,patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixC^S))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(consL,consR,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
              cminC,cmaxC,ixI^L,ixC^L,idims,fCD)
      endif ! Calculate the CD flux

      do iw=nc_hydro_lo,nc_hydro_hi
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S,iw) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
                 max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
                 (consR(ixC^S,iw) - consL(ixC^S,iw)))
         else
            where(patchf(ixC^S)==-2)
               fLC(ixC^S,iw)=fLC(ixC^S,iw)
            elsewhere(abs(patchf(ixC^S))==1)
               fLC(ixC^S,iw)=fCD(ixC^S,iw)
            elsewhere(patchf(ixC^S)==2)
               fLC(ixC^S,iw)=fRC(ixC^S,iw)
            elsewhere(patchf(ixC^S)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixC^S,iw)=Fhll(ixC^S,iw)
            elsewhere(patchf(ixC^S)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixC^S,iw) = 0.5d0*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                    -tvdlfeps * max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                    (consR(ixC^S,iw)-consL(ixC^S,iw)))
            endwhere
         end if

         fC(ixC^S,iw,idims)=fLC(ixC^S,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    subroutine positivity_preserving_limiter()
      !use mod_eos, only: small_D, small_tau
      implicit none
      double precision, dimension(ixI^S)   :: inv_lambda
      double precision, dimension(ixI^S)   :: theta, theta_tau

       do idims= idims^LIM
         if (slab_uniform) then
            !inv_lambda(ixI^S) = -dxinv(idims)
            inv_lambda(ixI^S) = dxdim(idims)/(qdt)
         else
            inv_lambda(ixI^S) = block%dvolume(ixI^S)/(qdt)
         end if

         call get_theta(ixI^L,ixC^L,idims,epsD,inv_lambda(ixI^S),sCT%cons(ixI^S,D_),fClow(ixI^S,D_,idims),fC(ixI^S,D_,idims),theta)
         call get_theta(ixI^L,ixC^L,idims,epstau,inv_lambda(ixI^S),sCT%cons(ixI^S,tau_),fClow(ixI^S,tau_,idims),fC(ixI^S,tau_,idims),theta_tau)

         theta(ixC^S) = min(theta(ixC^S),theta_tau(ixC^S))

         do iw = nc_hydro_lo,nc_hydro_hi
            fC(ixC^S,iw,idims) = theta(ixC^S)*fC(ixC^S,iw,idims) &
                 + (1.0d0-theta(ixC^S))*fClow(ixC^S,iw,idims)
         end do
       end do ! Next idims
    end subroutine positivity_preserving_limiter

  end subroutine finite_volume

  !> Determine the upwinded consL(ixL) and consR(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(typelimiter_in,rec_from,rec_to,ixI^L,ixL^L,ixR^L,idims,w,consL,consR,primL,primR,x,dxdim)
    use mod_physics
    use mod_eos
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: typelimiter_in
    integer, intent(in) :: rec_from, rec_to
    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
    double precision, intent(in) :: dxdim
    double precision, dimension(ixI^S,1:nprim) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:ncons) :: consL, consR
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S,1:nprim) :: primL, primR
    double precision, dimension(ixI^S,1:ndim) :: x

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    logical            :: fattening = .False.

    select case (typelimiter_in)
    case (limiter_venk)
       call venklimiter(rec_from,rec_to,ixI^L,ixL^L,idims,dxdim,w,primL,primR) 
    case (limiter_mp5)
       call MP5limiter(rec_from,rec_to,ixI^L,ixL^L,idims,w,primL,primR)
    case (limiter_weno3)
       call WENO3limiter(rec_from,rec_to,ixI^L,ixL^L,idims,w,primL,primR)
    case (limiter_weno5)
       call WENO5limiter(rec_from,rec_to,ixI^L,ixL^L,idims,dxdim,w,primL,primR,1)
    case (limiter_wenoz5)
       call WENO5limiter(rec_from,rec_to,ixI^L,ixL^L,idims,dxdim,w,primL,primR,2)
    case (limiter_wenozp5)
       call WENO5limiter(rec_from,rec_to,ixI^L,ixL^L,idims,dxdim,w,primL,primR,3)
    case (limiter_weno7)
       call WENO7limiter(rec_from,rec_to,ixI^L,ixL^L,idims,w,primL,primR,1)
    case (limiter_mpweno7)
       call WENO7limiter(rec_from,rec_to,ixI^L,ixL^L,idims,w,primL,primR,2)
    case (limiter_exeno7)
       call exENO7limiter(rec_from,rec_to,ixI^L,ixL^L,idims,w,primL,primR)
    case (limiter_ppm)
       if ( (rec_from >= nhydro_lo) .and. (rec_to <= nhydro_hi) ) fattening = .True.
       ! our fattening is only available for hydro
       call PPMlimiter(rec_from,rec_to,ixI^L,ixM^LL,idims,w,w,primL,primR,fattening)
    case default
       jxR^L=ixR^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=rec_from, rec_to
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             primL(ixL^S,iw)=dlog10(primL(ixL^S,iw))
             primR(ixR^S,iw)=dlog10(primR(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,ldw,rdw)
          primL(ixL^S,iw)=primL(ixL^S,iw)+0.5d0*ldw(ixL^S)
          primR(ixR^S,iw)=primR(ixR^S,iw)-0.5d0*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             primL(ixL^S,iw)=10.0d0**primL(ixL^S,iw)
             primR(ixR^S,iw)=10.0d0**primR(ixR^S,iw)
          end if
       end do
    end select

    if ( (rec_from >= nhydro_lo) .and. (rec_to <= nhydro_hi) ) then
       ! reconstructing hydro variables, 
       ! first check if all prim are OK
       call phys_handle_small_values(primL,x,ixI^L,ixL^L, .False., 'reconstruction')
       call phys_handle_small_values(primR,x,ixI^L,ixR^L, .False., 'reconstruction')

       ! we need to update rest of the prim variables
       call phys_update_eos(ixI^L,ixL^L,primL)
       call phys_update_eos(ixI^L,ixR^L,primR)
   
       call phys_to_conserved(ixI^L,ixL^L,consL,primL,x)
       call phys_to_conserved(ixI^L,ixR^L,consR,primR,x)
    end if
  end subroutine reconstruct_LR

  subroutine get_theta(ixI^L,ixC^L,idims,eps,inv_lambda,u,flow,fhigh,theta)
    use mod_global_parameters, only: kr
    integer, intent(in)                             :: ixI^L, ixC^L, idims
    double precision, intent(in)                    :: eps
    double precision, dimension(ixI^S), intent(in)  :: inv_lambda, u
    double precision, dimension(ixI^S), intent(in)  :: flow, fhigh
    double precision, dimension(ixI^S), intent(out) :: theta

    integer                                         :: ixCp^L, ixOp^L
    double precision, dimension(ixI^S)              :: tmp, thp, thm
    double precision, dimension(ixI^S)              :: diff_fdA

    ! Note: here we assume that u( i=0 ) is given
    ixCp^L=ixC^L+kr(idims,^D);
    ixOpmin^D=ixCmin^D; ixOpmax^D=ixCpmax^D;
    
    thm(ixC^S) = 1.0d0
    thp(ixC^S) = 1.0d0
    
    tmp(ixOp^S) = 0.5d0*inv_lambda(ixOp^S)*(u(ixOp^S)-eps)

    diff_fdA(ixC^S) = -flow(ixC^S) + fhigh(ixC^S)
    where (diff_fdA(ixC^S) == 0.0d0)
       diff_fdA(ixC^S) = epsilon(0.0d0) ! avoid flow = fhight case
    end where
    
    where (tmp(ixC^S) < fhigh(ixC^S))
       thm(ixC^S) = tmp(ixC^S) - flow(ixC^S)
       thm(ixC^S) = thm(ixC^S) / (diff_fdA(ixC^S))
    end where
    
    where (tmp(ixCp^S) < -fhigh(ixC^S))
       !thp(ixCp^S) = - tmp(ixCp^S) - flow(ixC^S)
       !thp(ixCp^S) = thp(ixCp^S) / (diff_fdA(ixC^S))
       thp(ixC^S) = - tmp(ixCp^S) - flow(ixC^S)
       thp(ixC^S) = thp(ixC^S) / (diff_fdA(ixC^S))
    end where

    theta(ixC^S) = min(thm(ixC^S),thp(ixC^S))
    theta(ixC^S) = min(max(theta(ixC^S),0.0d0),1.0d0)
    
  end subroutine get_theta

end module mod_finite_volume
