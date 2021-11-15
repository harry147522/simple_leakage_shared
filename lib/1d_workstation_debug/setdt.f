!>setdt  - set dt for all levels between levmin and levmax. 
!>         dtpar>0  --> use fixed dtpar for all level
!>         dtpar<=0 --> determine CFL limited timestep 
subroutine setdt()
  use mod_global_parameters
  use mod_physics
  use mod_usr_methods, only: usr_get_dt

  integer :: iigrid, igrid, ncycle, ncycle2, ifile, idim
  double precision :: dtnew, qdtnew, dtmin_mype, factor, dx1, dxmin1

  double precision :: dtmax, dxmin, cmax_mype, v(ixGlo1:ixGhi1)

  if (dtpar<=zero) then
     dtmin_mype=bigdouble
     cmax_mype = zero
     !$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx1)
     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        dtnew=bigdouble
        dx1=rnode(rpdx1_,igrid);
        dxlevel(1)=rnode(rpdx1_,igrid);
        saveigrid = igrid
        block=>ps(igrid)
        block%iw0=0

        call getdt_courant(ps(igrid)%prim,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qdtnew,&
           ps(igrid)%x)
        dtnew=min(dtnew,qdtnew)

        call phys_get_dt(ps(igrid)%prim,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qdtnew,dx1,&
           ps(igrid)%x)
        dtnew=min(dtnew,qdtnew)

        if (associated(usr_get_dt)) then
           call usr_get_dt(ps(igrid)%prim,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qdtnew,&
              dx1,ps(igrid)%x)
        end if

        dtnew          = min(dtnew,qdtnew)
        dtmin_mype     = min(dtmin_mype,dtnew)
        dt_grid(igrid) = dtnew
     end do
     !$OMP END PARALLEL DO
  else
     dtmin_mype=dtpar
  end if

  if (dtmin_mype<dtmin) then
     write(unitterm,*)"Error: Time step too small!", dtmin_mype
     write(unitterm,*)"on processor:", mype, "at time:", global_time," step:",&
         it
     write(unitterm,*)"Lower limit of time step:", dtmin
     crash=.true.
  end if

  if (slowsteps>it-it_init+1) then
     factor=one-(one-dble(it-it_init+1)/dble(slowsteps))**2
     dtmin_mype=dtmin_mype*factor
  end if


  dtmin_mype=min(dtmin_mype,time_max-global_time)

  if (dtpar<=zero) then
     call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, icomm,&
        ierrmpi)
  else
     dt=dtmin_mype
  end if

  if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),&
     1:nfile)<bigdouble))then
     dtmax = minval(ceiling(global_time/dtsave(1:nfile))*dtsave(1:nfile))-&
        global_time
     do ifile=1,nfile
        dtmax = min(tsave(isavet(ifile),ifile)-global_time,dtmax)
     end do
     if(dtmax > smalldouble)then 
       dt=min(dt,dtmax)
     else
       ! dtmax=0 means dtsave is divisible by global_time
       dt=min(dt,minval(dtsave(1:nfile)))
     end if      
  end if

  if(mype==0) then
    if(any(dtsave(1:nfile)<dt)) then
      write(unitterm,*) 'Warning: timesteps: ',dt,&
         ' exceeding output intervals ', dtsave(1:nfile)
    endif
  endif   

  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     dt_grid(igrid)=dt
  end do
  !$OMP END PARALLEL DO
       

  ! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
  ! so does GLM: 
  if(need_global_cmax) call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
     MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)

  contains

    !> compute CFL limited dt (for variable time stepping)
    subroutine getdt_courant(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,x)
      use mod_global_parameters
      use mod_physics, only: phys_get_cmax
      
      integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim), dtnew
      
      integer :: idims
      double precision :: courantmax, dxinv(1:ndim), courantmaxtot,&
          courantmaxtots
      double precision :: cmax(ixImin1:ixImax1), cmaxtot(ixImin1:ixImax1),&
          tmp(ixImin1:ixImax1)

      dtnew=bigdouble
      
      courantmax=zero
      courantmaxtot=zero
      courantmaxtots=zero
      
      dxinv(1)=one/dx1;
      
      cmaxtot(ixOmin1:ixOmax1)=zero
      
      do idims=1,ndim
        call phys_get_cmax(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,cmax)
        if(need_global_cmax) cmax_mype = max(cmax_mype,&
           maxval(cmax(ixOmin1:ixOmax1)))
        if(slab_uniform) then
          cmaxtot(ixOmin1:ixOmax1)=cmaxtot(ixOmin1:ixOmax1)+&
             cmax(ixOmin1:ixOmax1)*dxinv(idims)
          courantmax=max(courantmax,maxval(cmax(ixOmin1:ixOmax1)*dxinv(&
             idims)))
        else
          tmp(ixOmin1:ixOmax1)=cmax(ixOmin1:ixOmax1)/block%ds(ixOmin1:ixOmax1,&
             idims)
          cmaxtot(ixOmin1:ixOmax1)=cmaxtot(ixOmin1:ixOmax1)+&
             tmp(ixOmin1:ixOmax1)
          courantmax=max(courantmax,maxval(tmp(ixOmin1:ixOmax1)))
        end if
        courantmaxtot=courantmaxtot+courantmax
      end do
      
      select case (typecourant)
      case ('minimum')
         ! courantmax='max(c/dx)'
         if (courantmax>smalldouble)     dtnew=min(dtnew,&
            courantpar/courantmax)
      case ('summax')
         ! courantmaxtot='summed max(c/dx)'
         if (courantmaxtot>smalldouble)  dtnew=min(dtnew,&
            courantpar/courantmaxtot)
      case ('maxsum')
         ! courantmaxtots='max(summed c/dx)'
         courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ixOmin1:ixOmax1)))
         if (courantmaxtots>smalldouble) dtnew=min(dtnew,&
            courantpar/courantmaxtots)
      case default
         write(unitterm,*)'Unknown typecourant=',typecourant
         call mpistop("Error from getdt_courant: no such typecourant!")
      end select
      
    end subroutine getdt_courant

end subroutine setdt
