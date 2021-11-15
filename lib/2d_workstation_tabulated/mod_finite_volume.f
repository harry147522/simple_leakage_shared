!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume

contains

  !> finite volume method
  subroutine finite_volume(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idimsmin,idimsmax,qtC,sCT,qt,snew,sold,fC,fE,dx1,&
     dx2,x)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_cfc

    character(len=*), intent(in)                          :: method
    double precision, intent(in)                          :: qdt, qtC, qt, dx1,&
       dx2
    integer, intent(in)                                   :: ixImin1,ixImin2,&
       ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimsmin,idimsmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in) ::  x
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ncons,&
       1:ndim)     :: fC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)         :: fE

    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ncons) :: consL, consR
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim) :: primL, primR
    ! High order flux
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ncons) :: fLC, fRC
    ! Low order flux
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cminC
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    ! cell-face location coordinates
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: xi
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)               :: &
       patchf
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
       ixCRmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,&
       kxRmax2

    ! variables that related to positivity preserving limiter
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ncons,&
       1:ndim) :: fClow
    double precision, parameter :: epsD = 1.0d-16, epstau = 1.0d-16

    associate(primCT=>sCT%prim, consCT=>sCT%cons, cons_new=>snew%cons,&
       cons_old=>sold%cons)

    fC=0.d0
    fClow=0.d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in fv : Nonconforming input limits")

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
    dxdim(1)=dx1;dxdim(2)=dx2;
    do idims= idimsmin,idimsmax
       ! use interface value of w0 at idims
       block%iw0=idims

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

       if(stagger_grid) then
          ! ct needs all transverse cells
          ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
          ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
          ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
          ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;
       end if

       ! Determine stencil size
       ixCRmin1 = max(ixCmin1 - phys_wider_stencil,ixGlo1)
       ixCRmin2 = max(ixCmin2 - phys_wider_stencil,ixGlo2)
       ixCRmax1 = min(ixCmax1 + phys_wider_stencil,ixGhi1)
       ixCRmax2 = min(ixCmax2 + phys_wider_stencil,ixGhi2)

       ! get cell-face coordinates
       xi=x
       xi(ixImin1:ixImax1,ixImin2:ixImax2,idims)=xi(ixImin1:ixImax1,&
          ixImin2:ixImax2,idims)+0.5d0*sCT%dx(ixImin1:ixImax1,ixImin2:ixImax2,&
          idims)

       ! primR and primL are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'higher' direction.
       primR(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nprim)=primCT(kxRmin1:kxRmax1,&
          kxRmin2:kxRmax2,1:nprim)
       primL(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nprim)=primCT(kxCmin1:kxCmax1,&
          kxCmin2:kxCmax2,1:nprim)

       ! low order TVD-Lax-Friedrich flux first.
       if (positivity_preserving) then
         ! fixme: maybe this part can be faster
         call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,&
            ixCRmin2,ixCRmax1,ixCRmax2,consR,primR,xi)
         call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,&
            ixCRmin2,ixCRmax1,ixCRmax2,consL,primL,xi)
         call phys_get_flux(consL,primL,xi,ixImin1,ixImin2,ixImax1,ixImax2,&
            ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,fLC)
         call phys_get_flux(consR,primR,xi,ixImin1,ixImin2,ixImax1,ixImax2,&
            ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,fRC)
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixImin1,ixImin2,&
            ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC)
         call get_Riemann_flux_tvdlf()
         if ( slab_uniform ) then
           do iw=nc_hydro_lo,nc_hydro_hi
             fClow(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                idims)=fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)
           end do
         else 
           do iw=nc_hydro_lo,nc_hydro_hi
             fClow(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                idims)=fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                idims)*block%surfaceC(ixImin1:ixImax1,ixImin2:ixImax2,idims)
           end do
         end if
       end if

       if ( nmetric > 0 ) then
          ! fixme: make it faster, fewer index

          ! get cell-face metric
          if ( reconstruct_cfc ) then
             primR(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
                nmetric_lo:nmetric_hi)=primCT(kxRmin1:kxRmax1,kxRmin2:kxRmax2,&
                nmetric_lo:nmetric_hi)
             primL(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
                nmetric_lo:nmetric_hi)=primCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
                nmetric_lo:nmetric_hi)
             call reconstruct_LR(limiter_wenozp5, nmetric_lo, nmetric_hi,&
                 ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,ixCRmax1,&
                ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,primCT,&
                consL,consR,primL,primR,xi,dxdim(idims))
          else
             primL(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
                nmetric_lo:nmetric_hi)=primCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
                nmetric_lo:nmetric_hi)
             call cfc_metric_interpolation(ixImin1,ixImin2,ixImax1,ixImax2,&
                ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,primCT,x,primL,xi)
             primR(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
                nmetric_lo:nmetric_hi)=primL(ixCRmin1:ixCRmax1,&
                ixCRmin2:ixCRmax2,nmetric_lo:nmetric_hi)
          end if
       end if
       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(typelimiter, nhydro_lo, nhydro_lo+nreconstruct-1,&
           ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,&
          ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,primCT,consL,consR,primL,&
          primR,xi,dxdim(idims))

       ! special modification of left and right status before flux evaluation
       !call phys_modify_wLR(ixI^L,ixCR^L,consL,consR,primL,primR,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(consL,primL,xi,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,fLC)
       call phys_get_flux(consR,primR,xi,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixImin1,ixImin2,&
            ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC)
       else
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixImin1,ixImin2,&
            ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC)
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
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)=fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)*block%surfaceC(ixImin1:ixImax1,&
               ixImin2:ixImax2,idims)
          end do
       end if

    end do ! Next idims
    block%iw0=0

    ! If use positivity preserving limiter, with fC and fClow, work out the modify the flux
    if (positivity_preserving) then
       call positivity_preserving_limiter()
    end if

    if(associated(usr_set_flux)) then
       do idims= idimsmin,idimsmax
          call usr_set_flux(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,fC,xi)
       end do ! Next idims
    end if

    ! fixme: this might not right as I modified the sturcture
    if(stagger_grid) call phys_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,qdt,primCT,fC,fE,sCT,snew)

    do idims= idimsmin,idimsmax
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          do iw = nc_hydro_lo,nc_hydro_hi
            if ( flux_type(idims, iw) /= flux_nul ) then
               if (associated(phys_iw_methods(iw)%inv_capacity)) then
                 call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,&
                    ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, cons_new,&
                     inv_volume)
                 cons_new(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    iw)=cons_new(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    iw) + inv_volume * dxinv(idims)*(fC(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,iw,idims)-fC(hxOmin1:hxOmax1,&
                    hxOmin2:hxOmax2,iw,idims))
               else
                 cons_new(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    iw)=cons_new(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    iw) + dxinv(idims)*(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                    idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
               end if
            end if
          end do
        else
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,&
               ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, cons_new, inv_volume)
          else
            inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.0d0
          end if
          inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
             inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/block%dvolume(&
             ixOmin1:ixOmax1,ixOmin2:ixOmax2)

          do iw=nc_hydro_lo,nc_hydro_hi
            if ( flux_type(idims, iw) /= flux_nul ) then
               cons_new(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  iw)=cons_new(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  iw) - qdt * inv_volume(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2) * (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                  idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims)) 
            end if
          enddo
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') call tvdlimit2(method,qdt,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,ixOmax1,&
          ixOmax2,idims,consL,consR,cons_new,x,fC,dx1,dx2)

    end do ! Next idims

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,consCT,primCT,&
       cons_new,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       snew)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nprim,qtC,primCT,qt,&
       cons_new,x,.false.)

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()
      do iw=nc_hydro_lo,nc_hydro_hi
         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=0.5d0*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      fac = -0.5d0*tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=nc_hydro_lo,nc_hydro_hi

         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=0.5d0*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2, iw) + fac*(consR(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw)-consL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
         end if

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
          div(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      elsewhere
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
        fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = tvdlfeps*cminC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        div(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = 1/(cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      endwhere

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=nc_hydro_lo,nc_hydro_hi
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) = 0.5d0*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                dabs(cminC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))) * (consR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)-consL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)
               ! Add hll dissipation to the flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw) = (cmaxC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw)-cminC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2) * fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw) +fac(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2)*(consR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-consL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw))) * div(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)== 2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2, iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2, iw)
            endwhere
         endif

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ncons)     :: whll, Fhll, fCD
      double precision, dimension(ixImin1:ixImax1,&
         ixImin2:ixImax2)              :: lambdaCD

      patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') call phys_diffuse_hllcd(ixImin1,ixImin2,ixImax1,&
         ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,consL,consR,fLC,fRC,&
         patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)) call &
         phys_get_lCD(consL,consR,fLC,fRC,cminC,cmaxC,idims,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2, whll,Fhll,lambdaCD,&
         patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(consL,consR,whll,fRC,fLC,Fhll,patchf,lambdaCD,cminC,&
            cmaxC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
            ixCmax2,idims,fCD)
      endif ! Calculate the CD flux

      do iw=nc_hydro_lo,nc_hydro_hi
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) = 0.5d0 * (fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                abs(cminC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))) * (consR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - consL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))==1)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fCD(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=Fhll(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw) = 0.5d0*((fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)) -tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                   dabs(cminC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2))) * (consR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-consL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
            endwhere
         end if

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    subroutine positivity_preserving_limiter()
      !use mod_eos, only: small_D, small_tau
      implicit none
      double precision, dimension(ixImin1:ixImax1,&
         ixImin2:ixImax2)   :: inv_lambda
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)   :: theta,&
          theta_tau

       do idims= idimsmin,idimsmax
         if (slab_uniform) then
            !inv_lambda(ixI^S) = -dxinv(idims)
            inv_lambda(ixImin1:ixImax1,ixImin2:ixImax2) = dxdim(idims)/(qdt)
         else
            inv_lambda(ixImin1:ixImax1,ixImin2:ixImax2) = &
               block%dvolume(ixImin1:ixImax1,ixImin2:ixImax2)/(qdt)
         end if

         call get_theta(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
            ixCmax1,ixCmax2,idims,epsD,inv_lambda(ixImin1:ixImax1,&
            ixImin2:ixImax2),sCT%cons(ixImin1:ixImax1,ixImin2:ixImax2,D_),&
            fClow(ixImin1:ixImax1,ixImin2:ixImax2,D_,idims),fC(ixImin1:ixImax1,&
            ixImin2:ixImax2,D_,idims),theta)
         call get_theta(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
            ixCmax1,ixCmax2,idims,epstau,inv_lambda(ixImin1:ixImax1,&
            ixImin2:ixImax2),sCT%cons(ixImin1:ixImax1,ixImin2:ixImax2,tau_),&
            fClow(ixImin1:ixImax1,ixImin2:ixImax2,tau_,idims),&
            fC(ixImin1:ixImax1,ixImin2:ixImax2,tau_,idims),theta_tau)

         theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = min(theta(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2),theta_tau(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

         do iw = nc_hydro_lo,nc_hydro_hi
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
               idims) = theta(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
               idims) + (1.0d0-theta(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))*fClow(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
               idims)
         end do
       end do ! Next idims
    end subroutine positivity_preserving_limiter

  end subroutine finite_volume

  !> Determine the upwinded consL(ixL) and consR(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(typelimiter_in,rec_from,rec_to,ixImin1,ixImin2,&
     ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,&
     ixRmax2,idims,w,consL,consR,primL,primR,x,dxdim)
    use mod_physics
    use mod_eos
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: typelimiter_in
    integer, intent(in) :: rec_from, rec_to
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2, idims
    double precision, intent(in) :: dxdim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nprim) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ncons) :: consL, consR
    ! left and right constructed status in primitive form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nprim) :: primL, primR
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: x

    integer            :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, iw
    double precision   :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
        rdw(ixImin1:ixImax1,ixImin2:ixImax2), dwC(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    logical            :: fattening = .False.

    select case (typelimiter_in)
    case (limiter_venk)
       call venklimiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,dxdim,w,primL,primR) 
    case (limiter_mp5)
       call MP5limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,&
          ixLmin2,ixLmax1,ixLmax2,idims,w,primL,primR)
    case (limiter_weno3)
       call WENO3limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,w,primL,primR)
    case (limiter_weno5)
       call WENO5limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,dxdim,w,primL,primR,1)
    case (limiter_wenoz5)
       call WENO5limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,dxdim,w,primL,primR,2)
    case (limiter_wenozp5)
       call WENO5limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,dxdim,w,primL,primR,3)
    case (limiter_weno7)
       call WENO7limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,w,primL,primR,1)
    case (limiter_mpweno7)
       call WENO7limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,w,primL,primR,2)
    case (limiter_exeno7)
       call exENO7limiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2,idims,w,primL,primR)
    case (limiter_ppm)
       if ( (rec_from >= nhydro_lo) .and. (rec_to <= nhydro_hi) ) fattening = &
          .True.
       ! our fattening is only available for hydro
       call PPMlimiter(rec_from,rec_to,ixImin1,ixImin2,ixImax1,ixImax2,ixMlo1,&
          ixMlo2,ixMhi1,ixMhi2,idims,w,w,primL,primR,fattening)
    case default
       jxRmin1=ixRmin1+kr(idims,1);jxRmin2=ixRmin2+kr(idims,2)
       jxRmax1=ixRmax1+kr(idims,1);jxRmax2=ixRmax2+kr(idims,2);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=ixLmin1-kr(idims,1)
       ixCmin2=ixLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

       do iw=rec_from, rec_to
          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=dlog10(w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw))
             primL(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=dlog10(primL(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw))
             primR(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=dlog10(primR(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw))
          end if

          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,typelimiter_in,ldw,rdw)
          primL(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=primL(ixLmin1:ixLmax1,&
             ixLmin2:ixLmax2,iw)+0.5d0*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2)
          primR(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=primR(ixRmin1:ixRmax1,&
             ixRmin2:ixRmax2,iw)-0.5d0*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)

          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=10.0d0**w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw)
             primL(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=10.0d0**primL(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)
             primR(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=10.0d0**primR(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)
          end if
       end do
    end select

    if ( (rec_from >= nhydro_lo) .and. (rec_to <= nhydro_hi) ) then
       ! reconstructing hydro variables, 
       ! first check if all prim are OK
       call phys_handle_small_values(primL,x,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixLmin1,ixLmin2,ixLmax1,ixLmax2, .False., 'reconstruction')
       call phys_handle_small_values(primR,x,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixRmin1,ixRmin2,ixRmax1,ixRmax2, .False., 'reconstruction')

       ! we need to update rest of the prim variables
       call phys_update_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,primL)
       call phys_update_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,&
          ixRmax1,ixRmax2,primR)
   
       call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,consL,primL,x)
       call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,&
          ixRmax1,ixRmax2,consR,primR,x)
    end if
  end subroutine reconstruct_LR

  subroutine get_theta(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
     ixCmax2,idims,eps,inv_lambda,u,flow,fhigh,theta)
    use mod_global_parameters, only: kr
    integer, intent(in)                             :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, idims
    double precision, intent(in)                    :: eps
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)  :: inv_lambda, u
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)  :: flow, fhigh
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(out) :: theta

    integer                                         :: ixCpmin1,ixCpmin2,&
       ixCpmax1,ixCpmax2, ixOpmin1,ixOpmin2,ixOpmax1,ixOpmax2
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)              :: tmp, thp, thm
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)              :: diff_fdA

    ! Note: here we assume that u( i=0 ) is given
    ixCpmin1=ixCmin1+kr(idims,1);ixCpmin2=ixCmin2+kr(idims,2)
    ixCpmax1=ixCmax1+kr(idims,1);ixCpmax2=ixCmax2+kr(idims,2);
    ixOpmin1=ixCmin1;ixOpmin2=ixCmin2; ixOpmax1=ixCpmax1;ixOpmax2=ixCpmax2;
    
    thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = 1.0d0
    thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = 1.0d0
    
    tmp(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2) = &
       0.5d0*inv_lambda(ixOpmin1:ixOpmax1,&
       ixOpmin2:ixOpmax2)*(u(ixOpmin1:ixOpmax1,ixOpmin2:ixOpmax2)-eps)

    diff_fdA(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -flow(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2) + fhigh(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
    where (diff_fdA(ixCmin1:ixCmax1,ixCmin2:ixCmax2) == 0.0d0)
       diff_fdA(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = epsilon(0.0d0) !avoid flow = fhight case
    end where
    
    where (tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2) < fhigh(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2))
       thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = tmp(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2) - flow(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       thm(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = thm(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2) / (diff_fdA(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
    end where
    
    where (tmp(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2) < -fhigh(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2))
       !thp(ixCp^S) = - tmp(ixCp^S) - flow(ixC^S)
       !thp(ixCp^S) = thp(ixCp^S) / (diff_fdA(ixC^S))
       thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = - tmp(ixCpmin1:ixCpmax1,&
          ixCpmin2:ixCpmax2) - flow(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = thp(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2) / (diff_fdA(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
    end where

    theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = min(thm(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2),thp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
    theta(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = min(max(theta(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2),0.0d0),1.0d0)
    
  end subroutine get_theta

end module mod_finite_volume
