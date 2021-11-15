module mod_usr
  use mod_physics
!  use mod_grhd
  use mod_grhd_ccsn_leakage
  use mod_grhd_ccsn_leakage_simple_leakage
!  use mod_eos_hybrid
!  use mod_eos_polytrope
  use mod_eos_tabulated
  use mod_map_profile
  use mod_multigrid_coupling
  use mod_cfc
  use mod_yeofrho  
 
  implicit none

  logical, parameter  ::  initialize_metric = .True.
  logical, parameter  ::  evolve_metric = .True.
  !logical, parameter  ::  initialize_metric = .False.
  !logical, parameter  ::  evolve_metric = .False.
  logical :: do_ye_of_rho = .true.   ! true
  logical :: do_yeofrhofit = .true.   ! true
  logical :: do_highcorrection = .false.

!  logical :: bounce = .False.
!  double precision :: t_bounce
contains

  subroutine usr_init()
    !use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc
    use mod_physics, only: phys_step_adv_global

!    character*128  ::  profile_path = './profile/s15s7b2.short'
    character*128  ::  profile_path = './profile/s40WW95.dat'    !evan's40ww95
!    character*128  ::  profile_path = './profile/s40presn'
!    character*128  ::  profile_path = './profile/s20presn_evan'

    integer  ::  Nr   ! number of grid of the profile
    double precision :: rho_index
    

    usr_init_one_grid => ccsn_init_one_grid
    usr_improve_initial_condition => ccsn_improve_initial_condition
    usr_print_log => printlog
    usr_refine_grid => my_refine
    usr_process_adv_global => my_poststep_analysis
    phys_step_adv_global => my_physstep_analysis


    call set_coordinate_system("spherical")
!    call grhd_activate()
    call grhd_ccsn_leakage_activate()
    call grhd_simple_leakage_activate()

    call eos_tabulated_activate()
!    call eos_hybrid_activate()
!    call eos_polytrope_activate
    ! use atmosphere
    call eos_atmo_activate()

    !call initialize_gmunu()

    call mod_read_profile(profile_path)

    ! setup atmosphere density from the profile
!    call eos_initialize_atmo(minval(prof%rho))
!   2.42985d-15  = 1.5d3  cgs
    call eos_initialize_atmo(2.42985d-15)
if (mype == 0) then
    write(*,*) "atmosphere quantitiies"
    write(*,*) "small_rho_thr    small_rho     small_eps  (code unit)"
    write(*,*) small_rho_thr, small_rho, small_eps
endif

    ! to use metric solver, we need to activate multigrid
    use_multigrid = .true.
    if (evolve_metric) then
       ! to use metric solver, we need to activate multigrid
       call cfc_solver_activate()
    end if

    call set_eos_epsmin_max


    if (do_ye_of_rho) then
    call init_yeofrho
    endif



  end subroutine usr_init

  ! Initialize one grid
  subroutine ccsn_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    use mod_global_parameters
    use mod_eos

    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nprim)

    integer  ::  ix1

    w = 0.0d0

    ! initally flat spacetime
    w(ix^S, alp_) = 1.0d0
    w(ix^S, psi_) = 1.0d0


    {do ix^D = ixM^LLIM^D \}
       call mod_map(w(ix1,rho_),x(ix1,r_),prof%rho,prof%radius,prof%Nr)
       call mod_map(w(ix1,W_vel(1)),x(ix1,r_),prof%vel,prof%radius,prof%Nr)
       call mod_map(w(ix1,temp_),x(ix1,r_),prof%temp,prof%radius,prof%Nr)
       call mod_map(w(ix1,ye_),x(ix1,r_),prof%ye,prof%radius,prof%Nr)
       call mod_map(w(ix1,press_),x(ix1,r_),prof%press,prof%radius,prof%Nr)
       call mod_map(w(ix1,eps_),x(ix1,r_),prof%eps,prof%radius,prof%Nr)
       
    {end do^D&\}


    {do ix^D = ixM^LLIM^D \}
       ! get_eps
       call eos_get_eps_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_),&
                temp=w(ix1,temp_), ye=w(ix1, ye_))
!       call eos_get_pressure_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
!       call eos_get_cs2_one_grid(w(ix1, cs2_),w(ix1, rho_),w(ix1, eps_))
       call eos_eps_get_all_one_grid(rho = w(ix1,rho_),temp = w(ix1,temp_), &
                                 ye = w(ix1,ye_),eps = w(ix1,eps_),&
                                 prs = w(ix1,press_),ent = w(ix1,ent_), &
                                 cs2 = w(ix1,cs2_))

    {end do^D&\}

    ! fill atmosphere
    ! here we assume that the smallest rho in the profile is the atmosphere
    where (w(ix^S,rho_) <= small_rho_thr )
       w(ix^S,rho_)   = small_rho
       w(ix^S,W_vel(1)) = 0.0d0
       w(ix^S,press_)   = small_press
       w(ix^S,eps_)     = small_eps
    end where

    !perturbation
    !w(ix^S, veloc(1)) = 1.0d-8
    !where (x(ix^S,rho_) <= 1.0d1 )
       !w(ix^S,veloc(1)) = -5.0d-5 * dsin(dpi * x(ix^S,1)/1.0d1)
    !end where

  end subroutine ccsn_init_one_grid

  !> before the main loop, we improve the initial data here
  subroutine ccsn_improve_initial_condition()
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_cfc

    integer              :: ix^D
    integer              :: igrid, iigrid

    ! deallocate profile varibles to save memories
    deallocate(prof%radius)
    deallocate(prof%mass)
    deallocate(prof%rho)
    deallocate(prof%press)
    deallocate(prof%temp)
    deallocate(prof%eps)
    deallocate(prof%vel)
    deallocate(prof%ye)
    deallocate(prof%omega)

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       {do ix^D = ixM^LLIM^D \}
          ! reset eps
          call eos_get_eps_one_grid(ps(igrid)%prim(ix^D, press_),ps(igrid)%prim(ix^D, rho_),ps(igrid)%prim(ix^D, eps_),&
          temp=ps(igrid)%prim(ix^D, temp_), ye=ps(igrid)%prim(ix^D, ye_))
          ps(igrid)%prim(ix^D, eps_) = max( ps(igrid)%prim(ix^D, eps_), eos_epsmin )

          ! fill atmosphere
          if (ps(igrid)%prim(ix^D,rho_) <= small_rho_thr ) then
             ps(igrid)%prim(ix^D,rho_)   = small_rho
             ps(igrid)%prim(ix^D,eps_)   = small_eps
             ps(igrid)%prim(ix^D,press_)   = small_press
             ps(igrid)%prim(ix^D,W_vel(1)) = 0.0d0
             !ps(igrid)%prim(ix^D,W_vel(2)) = 0.0d0
             !ps(igrid)%prim(ix^D,W_vel(3)) = 0.0d0
          end if
       {end do^D&\}
   
       ! update rest of the prim
       call phys_update_eos(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,1:nprim))
    end do

    ! initialize metric
    !if (evolve_metric) then
       if ( initialize_metric ) call cfc_metric_init()
    !end if
  end subroutine ccsn_improve_initial_condition

  subroutine my_physstep_analysis(iit, qt)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    integer              :: igrid, iigrid
!    double precision     :: rho_max, rho_max_local
!    double precision     :: eps_max, eps_max_local
!code test
  if (bounce) then
        ! 0.005 second after bounce  (the shock still propagating)
! code test
     if (do_ye_of_rho .and. (global_time .lt. t_bounce*time_gf+0.005d0*time_gf)) then
            call adjust_ye      !after bounce only adjust the uns  (the shock still propagating)hocked region
     endif
  else
   !not bounce
    if (do_ye_of_rho) then
            call adjust_ye
    endif
  endif


  end subroutine 

  subroutine my_poststep_analysis(iit, qt)
    use mod_global_parameters
    use mod_eos
    use mod_eos_tabulated
!    use mod_grhd_ccsn_leakage_phys_parameters
    use mod_grhd_ccsn_leakage_phys_parameters
!    use mod_input_output, only : grid_rmax
!    use mod_eos_hybrid, only: rho_nuc


    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer              :: ix^D, ix^L, idir, pns_radius_local_i, veloc1_min_local_i
    integer              :: igrid, iigrid
    double precision     :: rho_max, ent_max, rho_max_local, ent_max_local
    double precision     :: shock_veloc, prev_shock_veloc, ye_min_local, ye_min
    double precision     :: pns_radius_local, pns_radius_local_tmp
    double precision     :: shock_radius_local
    double precision     :: veloc1_min_local, veloc1_min
    double precision     :: out_ye, delta_ye, delta_S, out_munu, r_max

    double precision     :: veloc_local(ixG^T,1:ndir), diff(ixG^T)
    double precision     :: lfac(ixG^T)
    double precision     :: gamma(ixG^T,1:3,1:3)

      character*1024 filename
     filename = trim(base_filename)//".outinfo"

!bounce = .true.
  ent_max = 0.0d0
  ent_max_local = 0.0d0
!for code test
  rho_max = 0.0d0
  rho_max_local = 0.0d0
!Bounce checker
  if (.not. bounce) then
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

      {do ix^D = ixM^LL \}
        if (ps(igrid)%x(ix^D,1) .le. 3.0d6*length_gf) then
           if (ps(igrid)%prim(ix^D,ent_) .ge. ent_max_local) then
           ent_max_local = ps(igrid)%prim(ix^D,ent_)

           endif
        endif
      {enddo^D&\}
    end do
    call MPI_ALLREDUCE(ent_max_local, ent_max, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)

      if (ent_max >= 3.0d0) then

!for code test
!    do iigrid = 1, igridstail
!       igrid = igrids(iigrid)
!      {do ix^D = ixM^LL \}
!           if (ps(igrid)%prim(ix^D,rho_) .ge. rho_max_local) then
!           rho_max_local = ps(igrid)%prim(ix^D,rho_)
!           endif
!      {enddo^D&\}
!    end do
!    call MPI_ALLREDUCE(rho_max_local, rho_max, 1, mpi_double_precision, &
!          MPI_MAX, icomm, ierrmpi)
!
!     if (mype == 0) then
!        write(*,*) rho_max , "rhoMax"
!     endif
!  if (.not. bounce) then
!       if (rho_max >= 5.0d-2) then
!        stop  "bounce then stop excess rho=5.0d-2"
!!!!
       bounce = .True.
    write(*,*) "bounce! entropy > 3.0d0", "global_time = ", global_time
        t_bounce = global_time/time_gf   ! in second
        
!        t_dump = global_time/time_gf + 1d-5
       if (mype == 0) then
          write(*,'(a,i7,a,es12.4)') ' core bounce at it=',it,' global_time=',global_time
          !write(*,'(a,i7,a,i7,a,es12.4)') ' save a snapshot No.',&
          !     snapshotnext,' at it=',it,' global_time=',global_time

          open(unit=666,file='t_bounce.dat',status="unknown",&
               form='formatted',position="append")
          write(666,*) 'it=',it,' global_time=',global_time
          close(666)

          open(unit=666,file='savenow',status="unknown",&
               form='formatted',position="append")
          write(666,*) ''
          close(666)

!          open(unit=666,file='hydro_at_tbounce.dat',status="unknown",&
!               form='formatted',position="append")
!          write(666,*) 'it=',it,' global_time=',global_time
!          write(666,*) "    lx1(km)  lye  leps  lrho   cs2(code unit)  ltemp v1"
!          do iigrid = 1, igridstail
!           igrid = igrids(iigrid)
!             write(666,*) ps(igrid)%x(igrid,1)*1.0d-5/length_gf,ps(igrid)%prim(igrid,ye_), ps(igrid)%prim(igrid,eps_),&
!                          ps(igrid)%prim(igrid,rho_),ps(igrid)%prim(igrid,cs2_),ps(igrid)%prim(igrid,temp_),ps(igrid)%prim(igrid,veloc(1))
!          enddo
!          close(666)

!          open(unit=666,file='savenow',status="unknown",&
!               form='formatted',position="append")
!          write(666,*) ''
!          close(666)
       end if
      endif
  endif


!! shock radius
  if (.not.bounce) then
!  if (bounce) then  !just for try not real
   !  ishock = nghostcells+1

     shock_radius_local = 0.0d0
     shock_radius = 0.0d0
  else
!!NOT TESTED YET
!
!!1D
!! find !find
!!  whole domain, where is the lcoation of minimum value of v1


        shock_radius_local = 0.0d0
        shock_radius = 0.0d0
        lfac = 0.0d0
        veloc1_min_local = 9.9d10
        veloc1_min = 9.9d10

    do iigrid = 1, igridstail
        igrid = igrids(iigrid)
     !need change!!!!!!   need use  W2v2  method to find  W
!         call grhd_ccsn_leakage_get_intermediate_variables(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T, 1:nprim), ps(igrid)%x(ixG^T, 1:ndim), &
         call grhd_ccsn_leakage_get_intermediate_variables(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T, 1:nprim), ps(igrid)%x(ixG^T, 1:ndim), &
                   gamma=gamma(ixG^T,1:3,1:3), &
                   lfac=lfac(ixG^T) )
!   咁寫 出事
!        veloc1_min_local = minval(ps(igrid)%prim(ixM^T, W_vel(1)) / lfac(ixM^T))

        veloc1_min_local = min(minval(ps(igrid)%prim(ixM^T, W_vel(1)) / lfac(ixM^T)), veloc1_min_local)
!        write(*,*) veloc1_min_local, "veloc1_min_local"
    enddo

    call MPI_ALLREDUCE(veloc1_min_local, veloc1_min, 1, mpi_double_precision, &
          MPI_MIN, icomm, ierrmpi)
!        write(*,*) veloc1_min, "veloc1_min_"


    do iigrid = 1, igridstail
        igrid = igrids(iigrid)
!         call grhd_ccsn_leakage_get_intermediate_variables(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T, 1:nprim), ps(igrid)%x(ixG^T, 1:ndim), &
         call grhd_ccsn_leakage_get_intermediate_variables(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T, 1:nprim), ps(igrid)%x(ixG^T, 1:ndim), &
                   gamma=gamma(ixG^T,1:3,1:3), &
                   lfac=lfac(ixG^T) )

!        veloc_local(ixM^T, 1) = ps(igrid)%prim(ixM^T,W_vel(1))/lfac(ixM^T)
       {do ix^D = ixM^LL \}

!          if (veloc_local(ix^D,1) .eq. veloc1_min ) then
          if ((ps(igrid)%prim(ix^D, W_vel(1))/lfac(ix^D)) .eq. veloc1_min ) then
             shock_radius_local = max(shock_radius_local, ps(igrid)%x(ix^D, r_))
!             write(*,*) shock_radius, ix1,"shock-0radius inside loop"
          endif
       {enddo^D&\}
    enddo

    call MPI_ALLREDUCE(shock_radius_local, shock_radius, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)

     shock_radius = shock_radius/length_gf  ! in cm , need nghostcells or not?

    if (.not. bounce) then
!            write(*,*) shock_radius/1.d5, "shock_radius"
    else
!           write(*,*) shock_radius/1.d5, "shock_radius  bounce"
    endif
  endif



!! shock radius
!!pass
  if (.not.bounce) then
   !  ishock = nghostcells+1
     pns_radius_local = 0.0d0
  else

!NOT TESTED YET

!1D
! find !find
!  whole domain, where is the lcoation of the atmosphere of PNS (rho_cgs <  1.0d11)
        pns_radius = 0.0d0
        pns_radius_local = 0.0d0

    do iigrid = 1, igridstail
        igrid = igrids(iigrid)


!        ix1 = ixMlo1
!!         write(*,*) minval(ps(igrid)%prim(ixM^T, rho_)/rho_gf), minval(ps(igrid)%x(ixM^T,r_))/1.d5/length_gf
!        do while ( ps(igrid)%prim(ix1, rho_)/rho_gf .ge. 1.0d8 .and. ix1 .lt. ixMhi1 )
!          pns_radius_local = ps(igrid)%x(ix1, r_)
!          write(*,*) ps(igrid)%prim(ix1, rho_)/rho_gf/1.0d8, pns_radius_local/1.d5/length_gf, ix1, "rho pns_local"
!          ix1 = ix1+1
!        end do

       {do ix^D = ixM^LL \}

!          if (ps(igrid)%prim(ix1, rho_)/rho_gf .gt. 1.0d11) then
!             write(*,*) ps(igrid)%x(ix1, r_)/1.d5/length_gf ,"greater than 1.d11"
!                stop
!          endif

!          if (ps(igrid)%prim(ix1, rho_)/rho_gf .lt. 1.0d4 .and. .gt. 9.95d10) then
          if (ps(igrid)%prim(ix^D, rho_)/rho_gf .gt. 1.0d11) then
               pns_radius_local = max(pns_radius_local, ps(igrid)%x(ix^D, r_))
!                pns_radius_local = ps(igrid)%x(ix^D, r_)
              !  write(*,*) pns_radius_local/1.d5/length_gf, "pns_local"
!                write(*,*) ps(igrid)%prim(ix^D, rho_)/rho_gf/1.d8, pns_radius_local/1.d5/length_gf, "rho pns_local"
!             endif
          endif
       {enddo^D&\}

!       {do ix^D = ixM^LL \}
!          diff(ix^D) = abs(ps(igrid)%prim(ix^D, rho_)/1.0d8 - 1.0d0)
!       {enddo^D&\}
!        pns_radius_local_i = minloc(diff,1)
!        pns_radius_local = ps(igrid)%x(pns_radius_local_i, r_)
!
!  !     write(*,*) pns_radius_local/1.d5/length_gf, "pns_radius_local"
!                write(*,*) ps(igrid)%prim(pns_radius_local_i, rho_)/rho_gf/1.d8, pns_radius_local/1.d5/length_gf, "rho pns_local"
    end do

    call MPI_ALLREDUCE(pns_radius_local, pns_radius, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)

     pns_radius = pns_radius/length_gf

    if (.not. bounce) then
    !        write(*,*) pns_radius/1.d5, "pns_radius"
    else
!           write(*,*) pns_radius/1.d5, "pns_radius  bounce"
    endif
  endif
!
!  rho_max = 0.0d0
!  rho_max_local = 0.0d0
!
!    do iigrid = 1, igridstail
!       igrid = igrids(iigrid)
!      {do ix^D = ixM^LL \}
!           if (ps(igrid)%prim(ix^D,rho_) .ge. rho_max_local) then
!           rho_max_local = ps(igrid)%prim(ix^D,rho_)
!           endif
!      {enddo^D&\}
!    end do
!    call MPI_ALLREDUCE(rho_max_local, rho_max, 1, mpi_double_precision, &
!          MPI_MAX, icomm, ierrmpi)
!
!
!
!        if (mype==0) then
!              open(unit=666,file=filename,status="unknown",form='formatted',position="append")
!            if (.not. bounce) then
!               write(666,*) global_time/time_gf , dt/time_gf, rho_max,&
!                           shock_radius/1.d5, pns_radius/1.d5
!            else
!               write(666,*) global_time/time_gf , dt/time_gf, rho_max,&
!                           shock_radius/1.d5, pns_radius/1.d5, &
!                           t_bounce,  "Bounced"
!
!            endif
!              close(666)
!        endif
!!!!!!!!!!!!!!!!!!!!!!!!!code test no need

  end subroutine my_poststep_analysis


  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nprim), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    double precision             :: ratio(ixI^S)
    integer                      :: i_level, to_level
    double precision             :: phi_grav ! relativistic gravitational potential :=  1 - alp
    double precision             :: phi_grav_cut
    double precision, parameter  :: phi_grav_max = 0.2d0

    refine = 0
    coarsen = 0

    phi_grav = minval( w(ixO^S, alp_) )
    phi_grav = 1.0d0 - phi_grav
 
    to_level = refine_max_level
    phi_grav_cut = phi_grav_max 
    do i_level = refine_max_level-1, 1, -1
       phi_grav_cut = phi_grav_cut / 2.0d0
       if ( phi_grav < phi_grav_cut ) then
          to_level = i_level
       end if
    end do

    if ( level > to_level ) then
       refine = -1
       coarsen = 1
    else if ( level < to_level ) then
       refine = 1
       coarsen = -1
    end if

    ratio(ixO^S) = dx(1, level) / x(ixO^S,1) 

    if ( any( ratio(ixO^S) > 1.0d-2 ) ) then
       ! the resolution is too low
       refine = 1
       coarsen = -1
    else if ( all( ratio(ixO^S) < 1.0d-2 ) ) then
       ! the resolution is too high
       refine = -1
       coarsen = 1
    end if
  end subroutine my_refine

  subroutine adjust_ye
    use mod_eos
   ! use mod_grhd_ccsn_phys_parameters

    integer              :: igrid, iigrid
    integer              :: ix^D
    integer              :: ix^L

    double precision     :: out_ye, delta_ye, delta_S, out_munu, ltemp
!  yeofrho part
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

         {do ix^D = ixM^LL \}
!           ps(igrid)%prim(ix^D, rho_) = 1.0d12&rho_gf

            out_munu = 0.0d0
          ! nuchem(i,j,k) = 0.0d0
           if (ps(igrid)%prim(ix^D,rho_).gt.1.0d6*rho_gf) then
 !             if((.not.bounce).or.(i.gt.ishock(1))) then
               if((.not.bounce).or.(ps(igrid)%x(ix^D,r_).gt.shock_radius*length_gf)) then
        !code test
         !     if(.not.bounce) then
                 if (do_yeofrhofit) then
                    call fit_ye(ps(igrid)%prim(ix^D,rho_)/rho_gf,out_ye)
                 else
                    call interpolate_ye(ps(igrid)%prim(ix^D,rho_)/rho_gf,out_ye)
                 endif
                 delta_ye = min(0.0d0, out_ye-ps(igrid)%prim(ix^D,ye_))


                !    call eos(out_munu,rho(i,j,k),temp(i,j,k),ye(i,j,k),eps(i,j,k), &!no need temp ??
                !              keytemp,keyerr,eosflag,eos_rf_prec)
                   ltemp = ps(igrid)%prim(ix^D, temp_)   !  do not want to update prim(temp_)
                   call eos_eps_get_all_one_grid(rho = ps(igrid)%prim(ix^D,rho_), temp=ltemp, &
                                     ye = ps(igrid)%prim(ix^D,ye_),eps = ps(igrid)%prim(ix^D,eps_),&
                                     munu = out_munu)


!                    nuchem(i,j,k) = out_munu
!                  if (mpye==0) then
!                        write(*,*) out_munu,"out_munu"
!                  endif

                    if (out_munu.lt.10.0d0) then
                       delta_S = 0.0d0
                    elseif (ps(igrid)%prim(ix^D,rho_).lt.2.0d12*rho_gf) then
                       delta_S = -delta_ye*(out_munu-10.0d0)/ps(igrid)%prim(ix^D,temp_)
                    else
                       delta_S = 0.0d0
                    endif
                    
!                  if (mpye==0) then
!                        write(*,*)  delta_S , delta_ye, "delta_S   delta_ye"
!                  endif
                    !this updates temp using constant entropy, then finds new energy
                    !associated with that new temp, the one call does both.


                    ps(igrid)%prim(ix^D,ent_) =  ps(igrid)%prim(ix^D,ent_) + delta_S
                    ps(igrid)%prim(ix^D,ye_) = ps(igrid)%prim(ix^D,ye_) + delta_ye

                     !using new entropy to find all new things
!                    write(*,*)  delta_S, delta_ye, "delta"
    !                write(*,*) "pass good way ye)ofrho.f90 !!!!!!!!!!!!!"

                   call eos_temp_get_all_one_grid(rho = ps(igrid)%prim(ix^D,rho_),&
                                     ye = ps(igrid)%prim(ix^D,ye_), temp = ps(igrid)%prim(ix^D,temp_), &
                                     eps = ps(igrid)%prim(ix^D,eps_), prs =ps(igrid)%prim(ix^D,press_),&
                                     cs2 = ps(igrid)%prim(ix^D,cs2_), ent = ps(igrid)%prim(ix^D, ent_))

!    broken find_from_entropy!!!  need fix
!                   call eos_get_all_from_ent_one_grid(rho = ps(igrid)%prim(ix^D,rho_), ent = ps(igrid)%prim(ix^D,ent_),&
!                                     ye = ps(igrid)%prim(ix^D,ye_), temp = ps(igrid)%prim(ix^D,temp_), &
!                                     eps = ps(igrid)%prim(ix^D,eps_), prs =ps(igrid)%prim(ix^D,press_),&
!                                     cs2 = ps(igrid)%prim(ix^D,cs2_))

                endif
             endif
         {enddo^D&\}
      enddo
!end part of yeofrho

  end subroutine adjust_ye


!   old printlog
  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters

    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(1:nprim), total_volume
    double precision     :: volume_coverage(refine_max_level)
    integer              :: nx^D, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    character(len=40)    :: fmt_string
    character(len=80)    :: filename
    character(len=2048)  :: line
    logical, save        :: opened  = .false.
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20

    integer              :: igrid, iigrid
    double precision     :: tmp
    integer              :: igrid_locate_r, ilocate_r
    double precision     :: locate_r


    tmp = bigdouble ! some non-zero number
    locate_r = 0.0d0 ! some non-zero number
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       locate_r = minval(ps(igrid)%x(ixM^T,1))
       if (locate_r <= tmp) then
          tmp = locate_r
          igrid_locate_r = igrid
          ilocate_r = minloc(ps(igrid)%x(ixM^T,1), dim=1,mask=(ps(igrid)%x(ixM^T,1) >0))
       end if
    end do

    call MPI_ALLREDUCE(tmp, locate_r, 1, mpi_double_precision, &
          MPI_MIN, icomm, ierrmpi)

    if (tmp /= locate_r) return ! Sorry, this processor is not the chosen one XD

    wmean(1:nprim) = ps(igrid_locate_r)%prim(ilocate_r,1:nprim)

    !if (mype == 0) then

       ! On first entry, open the file and generate the header
       if (.not. opened) then

          filename = trim(base_filename) // ".log"

          ! Delete the log when not doing a restart run
          if (restart_from_file == undefined) then
             open(unit=my_unit,file=trim(filename),status='replace')
             close(my_unit, status='delete')
          end if

          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)

          opened = .true.

          ! Start of file headern
          line = "it global_time dt"
          do level=1,nprim
             i = len_trim(line) + 2
             write(line(i:),"(a,a)") trim(prim_names(level)), " "
          end do

          ! Only write header if not restarting
          if (restart_from_file == undefined) then
            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
                 len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       ! Construct the line to be added to the log

       fmt_string = '(' // fmt_i // ',2' // fmt_r // ')'
       write(line, fmt_string) it, global_time, dt
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', nprim, fmt_r // ')'
       write(line(i:), fmt_string) wmean(1:nprim)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    !end if

  end subroutine printlog




!!!! new printlog
!  subroutine printlog()
!    use mod_input_output
!    use mod_timing
!    use mod_forest, only: nleafs, nleafs_active, nleafs_level
!    use mod_global_parameters
!    logical              :: fileopen
!    integer              :: i, iw, level
!    double precision     :: wmean(0:nprim), wbuffer(0:nprim), total_volume
!    double precision     :: volume_coverage(refine_max_level)
!    integer              :: nx^D, nc, ncells, dit
!    double precision     :: dtTimeLast, now, cellupdatesPerSecond
!    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
!    character(len=40)    :: fmt_string
!    character(len=80)    :: filename
!    character(len=2048)  :: line
!    logical, save        :: opened  = .false.
!    integer              :: amode, istatus(MPI_STATUS_SIZE)
!    integer, parameter   :: my_unit = 20
!
!    integer              :: igrid, iigrid
!    double precision     :: tmp(1:2), mass_buffer, mass
!    double precision     :: lfac(ixG^T)
!    double precision     :: rho_max, rho_max_local
!    double precision     :: alp_min, alp_min_local
!
!    mass_buffer = 0.0d0
!
!    ! calculate total rest mass
!    do iigrid = 1, igridstail
!       igrid = igrids(iigrid)
!       call phys_get_lfac2(ps(igrid), ixG^LL, ixM^LL, lfac(ixG^T))
!       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
!       mass_buffer = mass_buffer + &
!          sum(lfac(ixM^T) * ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%prim(ixM^T,rho_)*ps(igrid)%dvolume(ixM^T))
!    end do
!
!    ! gather all the results!
!    call MPI_ALLREDUCE(mass_buffer, mass, 1, mpi_double_precision, &
!          MPI_SUM, icomm, ierrmpi)
!
!    ! find max rho
!    rho_max_local = 0.0d0
!    do iigrid = 1, igridstail
!       igrid = igrids(iigrid)
!       rho_max_local = max(rho_max_local, maxval(ps(igrid)%prim(ixM^T,rho_)) )
!    end do
!    call MPI_ALLREDUCE(rho_max_local, rho_max, 1, mpi_double_precision, &
!          MPI_MAX, icomm, ierrmpi)
!
!    ! find min alpha
!    alp_min_local = 1.0d0
!    do iigrid = 1, igridstail
!       igrid = igrids(iigrid)
!       alp_min_local = min(alp_min_local, minval(ps(igrid)%prim(ixM^T,alp_)) )
!    end do
!    call MPI_ALLREDUCE(alp_min_local, alp_min, 1, mpi_double_precision, &
!          MPI_MIN, icomm, ierrmpi)
!
!    if (mype == 0) then
!
!       ! On first entry, open the file and generate the header
!       if (.not. opened) then
!
!          filename = trim(base_filename) // ".log"
!
!          ! Delete the log when not doing a restart run
!          if (restart_from_file == undefined) then
!             open(unit=my_unit,file=trim(filename),status='replace')
!             close(my_unit, status='delete')
!          end if
!
!          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
!          amode    = ior(amode,MPI_MODE_APPEND)
!
!          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
!               MPI_INFO_NULL, log_fh, ierrmpi)
!
!          opened = .true.
!
!          ! Start of file headern
!          line = "it global_time dt"
!          i = len_trim(line) + 2
!          write(line(i:),"(a,a)") trim("Mass"), " "
!          write(line(i:),"(a,a)") trim("rho_max"), " "
!          write(line(i:),"(a,a)") trim("alp_min"), " "
!
!          ! Only write header if not restarting
!          if (restart_from_file == undefined) then
!            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
!                 len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
!          end if
!       end if
!
!       ! Construct the line to be added to the log
!
!       fmt_string = '(' // fmt_i // ',2' // fmt_r // ')'
!       write(line, fmt_string) it, global_time, dt
!       i = len_trim(line) + 2
!
!       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
!       write(line(i:), fmt_string) mass
!       i = len_trim(line) + 2
!
!       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
!       write(line(i:), fmt_string) rho_max
!       i = len_trim(line) + 2
!
!       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
!       write(line(i:), fmt_string) alp_min
!       i = len_trim(line) + 2
!
!       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
!            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
!    end if
!
!  end subroutine printlog

end module mod_usr
