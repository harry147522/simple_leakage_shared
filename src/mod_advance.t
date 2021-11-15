!> Module containing all the time stepping schemes
module mod_advance

  implicit none
  private

  logical :: firstsweep, lastsweep

  !> Whether to conserve fluxes at the current sub-step
  logical :: fix_conserve_at_step = .true.

  public :: advance
  public :: process
  public :: process_advanced

contains

  !> Advance all the grids over one time step, including all sources
  subroutine advance(iit)
    use mod_global_parameters
    use mod_source, only: add_split_source
    use mod_physics, only: phys_req_diagonal &
                         , phys_to_primitive, phys_to_conserved &
                         , phys_step_adv_global 
    use mod_ghostcells_update, only: getbc
    use mod_cfc
    use mod_usr_methods
    integer, intent(in) :: iit
    integer :: iigrid, igrid, idimsplit
    integer :: ix^L

    if (positivity_preserving) then
       ! Note that the cons in one extra ghost cell are needed if we want to use PP limiter
       ix^L=ixM^LL^LADD1;
    else
       ix^L=ixM^LL;
    end if

    ! old solution values at t_n-1 no longer needed: make copy of w(t_n)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call phys_to_conserved(ixG^LL,ix^L,ps(igrid)%cons(ixG^T,1:ncons),ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))
       pso(igrid)%prim(ixG^T,1:nprim)=ps(igrid)%prim(ixG^T,1:nprim)
       pso(igrid)%cons(ix^S,1:ncons)=ps(igrid)%cons(ix^S,1:ncons)
       if(stagger_grid) pso(igrid)%prims=ps(igrid)%prims
    end do
    !$OMP END PARALLEL DO


! single core do
    if (associated(usr_before_addsource)) then
    ! leakage add source
!    write(*,*) "read usr_before_addsource"
    call usr_before_addsource(iit, global_time)
!     write(*,*) "passed usr_b4"
    end if

   
!    ! split source addition
    call add_split_source(prior=.true.)

    firstsweep=.true.
    if (dimsplit) then
       if ((iit/2)*2==iit .or. typedimsplit=='xy') then
          ! do the sweeps in order of increasing idim,
          do idimsplit=1,ndim
             lastsweep= idimsplit==ndim
             call advect(idimsplit,idimsplit)
          end do
       else
          ! If the parity of "iit" is odd and typedimsplit=xyyx,
          ! do sweeps backwards
          do idimsplit=ndim,1,-1
             lastsweep= idimsplit==1
             call advect(idimsplit,idimsplit)
          end do
       end if
    else
       ! Add fluxes from all directions at once
       lastsweep= .true.
       call advect(1,ndim)
    end if
    ! split source addition
    call add_split_source(prior=.false.)

   
    ! check if it's time to update metric
    if (cfc_update_flag()) then
       ! primitive will be updated when solving cfc
       call cfc_solve_metric()
    else
       ! transform the variables to primitive
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          call phys_to_primitive(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))
       end do
       !$OMP END PARALLEL DO
    end if

    if (associated(phys_step_adv_global)) then
!        write(*,*)"passed phys_step_adv_global"
       ! update the ghosts
       call getbc(global_time,0.0d0,ps,1,nprim,phys_req_diagonal)
       call phys_step_adv_global(iit,global_time)
    end if

    ! after updating primitive variables, update the ghosts
    call getbc(global_time,0.0d0,ps,1,nprim,phys_req_diagonal)

  end subroutine advance

  !> Advance all grids over one time step, but without taking dimensional
  !> splitting or split source terms into account
  subroutine advect(idim^LIM)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal

    integer, intent(in) :: idim^LIM
    integer             :: iigrid, igrid

    call init_comm_fix_conserve(idim^LIM,ncons)
    fix_conserve_at_step = time_advance .and. levmax>levmin

    ! copy w instead of wold because of potential use of dimsplit or sourcesplit
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ps1(igrid)%prim=ps(igrid)%prim
       ps1(igrid)%cons=ps(igrid)%cons
       ! Note: ps1,2,3,4 are needed because the metric needed to be passed.
       ! also, for GRRMHD, W_vel of preivous step is needed as an initial guess 
       if ( allocated(ps2(igrid)%prim) ) then
         ps2(igrid)%prim=ps(igrid)%prim
         ps2(igrid)%cons=ps(igrid)%cons
       end if
       if ( allocated(ps3(igrid)%prim) ) then
         ps3(igrid)%prim=ps(igrid)%prim
         ps3(igrid)%cons=ps(igrid)%cons
       end if
       if ( allocated(ps4(igrid)%prim) ) then
         ps4(igrid)%prim=ps(igrid)%prim
         ps4(igrid)%cons=ps(igrid)%cons
       end if
       if(stagger_grid) ps1(igrid)%prims=ps(igrid)%prims
    end do
    !$OMP END PARALLEL DO

    istep = 0

    select case (time_stepper)
    case ("onestep")
       select case (time_integrator)
       case ("Forward_Euler")
          call advect1(flux_scheme,1.0d0,idim^LIM,global_time,ps1,global_time,ps)

       case ("IMEX_Euler")
          call advect1(flux_scheme,1.0d0,idim^LIM,global_time,ps,global_time,ps1)
          call global_implicit_update(1.0d0,dt,global_time+dt,ps,ps1)

       case ("IMEX_SP")
          call global_implicit_update(1.0d0,dt,global_time,ps,ps1)
          call advect1(flux_scheme,1.0d0,idim^LIM,global_time,ps1,global_time,ps)

       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          write(unitterm,*) "Error in advect: Unknown time integration method"
          call mpistop("Correct time_integrator")
       end select

    case ("twostep")
       select case (time_integrator)
       case ("Predictor_Corrector")
          ! PC or explicit midpoint 
          ! predictor step
          fix_conserve_at_step = .false.
          call advect1(typepred1,0.5d0,idim^LIM,global_time,ps,global_time,ps1)
          ! corrector step
          fix_conserve_at_step = time_advance .and. levmax>levmin
          call advect1(flux_scheme,1.0d0,idim^LIM,global_time+0.5d0*dt,ps1,global_time,ps)

       case ("RK2_alfa")
          ! RK2 with alfa parameter, where rk_a21=alfa
          call advect1(flux_scheme,rk_a21, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons = ps(igrid)%cons+rk_b1*(ps1(igrid)%cons-ps(igrid)%cons)/rk_a21
             if(stagger_grid) ps(igrid)%prims = ps(igrid)%prims+(one-rk_b2)*(ps1(igrid)%prims-ps(igrid)%prims)/rk_a21
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_b2,idim^LIM,global_time+rk_a21*dt,ps1,global_time+rk_b1*dt,ps)

       case ("ssprk2")
          ! ssprk2 or Heun's method
          call advect1(flux_scheme,1.0d0, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons = 0.5d0*ps(igrid)%cons+0.5d0*ps1(igrid)%cons
             if(stagger_grid) ps(igrid)%prims = 0.5d0*ps(igrid)%prims+0.5d0*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,0.5d0,idim^LIM,global_time+dt,ps1,global_time+0.5d0*dt,ps)

       case ("IMEX_Midpoint")
          ! IMEX Midpoint (1,2,2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons = ps(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims = ps(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,0.5d0, idim^LIM,global_time,ps,global_time,ps1)
          call global_implicit_update(0.5d0,dt,global_time+0.5d0*dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons = ps(igrid)%cons+2.0d0*(ps2(igrid)%cons-ps1(igrid)%cons)
             if(stagger_grid) ps(igrid)%prims = ps(igrid)%prims+2.0d0*(ps2(igrid)%prims-ps1(igrid)%prims)
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,1.0d0, idim^LIM,global_time+0.5d0*dt,ps2,global_time,ps)

       case ("IMEX_Trapezoidal")
          call advect1(flux_scheme,1.0d0, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons = 0.5d0*(ps(igrid)%cons+ps1(igrid)%cons)
             if(stagger_grid) ps2(igrid)%prims = 0.5d0*(ps(igrid)%prims+ps1(igrid)%prims)
          end do
          !$OMP END PARALLEL DO
          call evaluate_implicit(global_time,ps)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons = ps1(igrid)%cons+0.5d0*dt*ps(igrid)%cons
             if(stagger_grid) ps1(igrid)%prims = ps1(igrid)%prims+0.5d0*dt*ps(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons = ps2(igrid)%cons+0.5d0*dt*ps(igrid)%cons
             if(stagger_grid) ps(igrid)%prims = ps2(igrid)%prims+0.5d0*dt*ps(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          ! fixme: maybe we dont need getbc here? at least for GRRMHD?
          !call getbc(global_time+dt,dt,ps1,1,nprim,phys_req_diagonal)
          call global_implicit_update(0.5d0,dt,global_time+dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons = ps(igrid)%cons+ps2(igrid)%cons-ps1(igrid)%cons
             if(stagger_grid) ps(igrid)%prims = ps(igrid)%prims+ps2(igrid)%prims-ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,0.5d0, idim^LIM,global_time+dt,ps2,global_time+0.5d0*dt,ps)

       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          write(unitterm,*) "Error in advect: Unknown time integration method"
          call mpistop("Correct time_integrator")
       end select

    case ("threestep")
       select case (time_integrator)
       case ("ssprk3")
          ! this is SSPRK(3,3) Gottlieb-Shu 1998 or SSP(3,2) depending on ssprk_order (3 vs 2)
          call advect1(flux_scheme,rk_beta11, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=rk_alfa21*ps(igrid)%cons+rk_alfa22*ps1(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims=rk_alfa21*ps(igrid)%prims+rk_alfa22*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,global_time+rk_alfa22*rk_c2*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=rk_alfa31*ps(igrid)%cons+rk_alfa33*ps2(igrid)%cons
             if(stagger_grid) ps(igrid)%prims=rk_alfa31*ps(igrid)%prims+rk_alfa33*ps2(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta33, idim^LIM,global_time+rk_c3*dt,ps2,global_time+(1.0d0-rk_beta33)*dt,ps)

       case ("RK3_BT")
          ! this is a general threestep RK according to its Butcher Table
          call advect1(flux_scheme,rk3_a21, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%cons=(ps1(igrid)%cons-ps(igrid)%cons)/rk3_a21
             if(stagger_grid) ps3(igrid)%prims=(ps1(igrid)%prims-ps(igrid)%prims)/rk3_a21
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=ps(igrid)%cons+rk3_a31*ps3(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims=ps(igrid)%prims+rk3_a31*ps3(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk3_a32, idim^LIM,global_time+rk3_c2*dt,ps1,global_time+rk3_a31*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=ps(igrid)%cons+rk3_b1*ps3(igrid)%cons &
                  +rk3_b2*(ps2(igrid)%cons-(ps(igrid)%cons+rk3_a31*ps3(igrid)%cons))/rk3_a32
             if(stagger_grid)then
                 ps(igrid)%prims=ps(igrid)%prims+rk3_b1*ps3(igrid)%prims &
                   +rk3_b2*(ps2(igrid)%prims-(ps(igrid)%prims+rk3_a31*ps3(igrid)%prims))/rk3_a32
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk3_b3, idim^LIM,global_time+rk3_c3*dt,ps2,global_time+(1.0d0-rk3_b3)*dt,ps)

       case ("IMEX_ARS3")
          ! this is IMEX scheme ARS3
          call advect1(flux_scheme,ars_gamma, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps4(igrid)%cons=(ps1(igrid)%cons-ps(igrid)%cons)/ars_gamma
             if(stagger_grid) ps4(igrid)%prims=(ps1(igrid)%prims-ps(igrid)%prims)/ars_gamma
          end do
          !$OMP END PARALLEL DO
          call global_implicit_update(ars_gamma,dt,global_time+ars_gamma*dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons=(ps2(igrid)%cons-ps1(igrid)%cons)/ars_gamma
             if(stagger_grid) ps1(igrid)%prims=(ps2(igrid)%prims-ps1(igrid)%prims)/ars_gamma
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%cons=ps(igrid)%cons+(ars_gamma-1.0d0)*ps4(igrid)%cons+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%cons
             if(stagger_grid) then
                ps3(igrid)%prims=ps(igrid)%prims+(ars_gamma-1.0d0)*ps4(igrid)%prims+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%prims
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,2.0d0*(1.0d0-ars_gamma), idim^LIM,global_time+ars_gamma*dt,ps2,global_time+(ars_gamma-1.0d0)*dt,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=ps1(igrid)%cons+(ps3(igrid)%cons-(ps(igrid)%cons+ &
               (ars_gamma-1.0d0)*ps4(igrid)%cons+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%cons))/(2.0d0*(1.0d0-ars_gamma))
             if(stagger_grid) then
             ps2(igrid)%prims=ps1(igrid)%prims+(ps3(igrid)%prims-(ps(igrid)%prims+ &
               (ars_gamma-1.0d0)*ps4(igrid)%prims+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%prims))/(2.0d0*(1.0d0-ars_gamma))
             endif
          end do
          !$OMP END PARALLEL DO
          call global_implicit_update(ars_gamma,dt,global_time+(1.0d0-ars_gamma)*dt,ps4,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=ps(igrid)%cons+0.5d0*ps2(igrid)%cons &
                +0.5d0*(ps4(igrid)%cons-ps3(igrid)%cons)/ars_gamma
             if(stagger_grid) then
                ps(igrid)%prims=ps(igrid)%prims+0.5d0*ps2(igrid)%prims &
                    +0.5d0*(ps4(igrid)%prims-ps3(igrid)%prims)/ars_gamma
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,0.5d0, idim^LIM,global_time+(1.0d0-ars_gamma)*dt,ps4,global_time+0.5d0*dt,ps)

       case ("IMEX_232")
          ! this is IMEX_ARK(2,3,2) or IMEX_SSP(2,3,2)
          call advect1(flux_scheme,imex_a21, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps4(igrid)%cons=(ps1(igrid)%cons-ps(igrid)%cons)/imex_a21
             ps3(igrid)%cons=ps(igrid)%cons
             if(stagger_grid) then
               ps4(igrid)%prims=(ps1(igrid)%prims-ps(igrid)%prims)/imex_a21
               ps3(igrid)%prims=ps(igrid)%prims
             endif
          end do
          !$OMP END PARALLEL DO
          call evaluate_implicit(global_time,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons=ps1(igrid)%cons+imex_ha21*dt*ps3(igrid)%cons
             if(stagger_grid) ps1(igrid)%prims=ps1(igrid)%prims+imex_ha21*dt*ps3(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          !call getbc(global_time+imex_a21*dt,dt,ps1,1,nprim,phys_req_diagonal)
          call global_implicit_update(imex_ha22,dt,global_time+imex_c2*dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=ps(igrid)%cons+imex_a31*ps4(igrid)%cons &
                +imex_b1*dt*ps3(igrid)%cons+imex_b2*(ps2(igrid)%cons-ps1(igrid)%cons)/imex_ha22
             if(stagger_grid) then
             ps(igrid)%prims=ps(igrid)%prims+imex_a31*ps4(igrid)%prims &
                +imex_b1*dt*ps3(igrid)%prims+imex_b2*(ps2(igrid)%prims-ps1(igrid)%prims)/imex_ha22
             endif
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%cons=ps1(igrid)%cons-imex_a21*ps4(igrid)%cons &
                -imex_ha21*dt*ps3(igrid)%cons+imex_b1*dt*ps3(igrid)%cons
             if(stagger_grid) then
             ps3(igrid)%prims=ps1(igrid)%prims-imex_a21*ps4(igrid)%prims &
                -imex_ha21*dt*ps3(igrid)%prims+imex_b1*dt*ps3(igrid)%prims
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,imex_a32, idim^LIM,global_time+imex_c2*dt,ps2,global_time+imex_a31*dt,ps)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=(ps(igrid)%cons-ps3(igrid)%cons-imex_a31*ps4(igrid)%cons)/imex_a32 &
                +(1.0d0-imex_b2/imex_a32)*(ps2(igrid)%cons-ps1(igrid)%cons)/imex_ha22
             if(stagger_grid) then
             ps2(igrid)%prims=(ps(igrid)%prims-ps3(igrid)%prims-imex_a31*ps4(igrid)%prims)/imex_a32 &
                +(1.0d0-imex_b2/imex_a32)*(ps2(igrid)%prims-ps1(igrid)%prims)/imex_ha22
             endif
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons=ps3(igrid)%cons+imex_b1*ps4(igrid)%cons+imex_b2*ps2(igrid)%cons
             if(stagger_grid) then
             ps1(igrid)%prims=ps3(igrid)%prims+imex_b1*ps4(igrid)%prims+imex_b2*ps2(igrid)%prims
             endif
          end do
          !$OMP END PARALLEL DO
          call global_implicit_update(imex_b3,dt,global_time+imex_c3*dt,ps2,ps)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=ps1(igrid)%cons+ps2(igrid)%cons-ps(igrid)%cons
             if(stagger_grid) then
             ps(igrid)%prims=ps1(igrid)%prims+ps2(igrid)%prims-ps(igrid)%prims
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,imex_b3, idim^LIM,global_time+imex_c3*dt,ps2,global_time+(1.0d0-imex_b3)*dt,ps)

       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          write(unitterm,*) "Error in advect: Unknown time integration method"
          call mpistop("Correct time_integrator")
       end select

    case ("fourstep")
       select case (time_integrator)
       case ("ssprk4")
          ! SSPRK(4,3) or SSP(4,2) depending on ssprk_order (3 vs 2)
          ! ssprk43: Strong stability preserving 4 stage RK 3rd order by Ruuth and Spiteri
          !    Ruuth & Spiteri J. S C, 17 (2002) p. 211 - 220
          !    supposed to be stable up to CFL=2.
          ! ssp42: stable up to CFL=3
          call advect1(flux_scheme,rk_beta11, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=rk_alfa21*ps(igrid)%cons+rk_alfa22*ps1(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims=rk_alfa21*ps(igrid)%prims+rk_alfa22*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,global_time+rk_alfa22*rk_c2*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons=rk_alfa31*ps(igrid)%cons+rk_alfa33*ps2(igrid)%cons
             if(stagger_grid) ps1(igrid)%prims=rk_alfa31*ps(igrid)%prims+rk_alfa33*ps2(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta33, idim^LIM,global_time+rk_c3*dt,ps2,global_time+rk_alfa33*rk_c3*dt,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=rk_alfa41*ps(igrid)%cons+rk_alfa44*ps1(igrid)%cons
             if(stagger_grid) ps(igrid)%prims=rk_alfa41*ps(igrid)%prims+rk_alfa44*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta44, idim^LIM,global_time+rk_c4*dt,ps1,global_time+(1.0d0-rk_beta44)*dt,ps)

       case ("rk4")
          ! the standard RK(4,4) method
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=ps(igrid)%cons
             ps3(igrid)%cons=ps(igrid)%cons
             if(stagger_grid) then
                ps2(igrid)%prims=ps(igrid)%prims
                ps3(igrid)%prims=ps(igrid)%prims
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,0.5d0, idim^LIM,global_time,ps,global_time,ps1)
          call advect1(flux_scheme,0.5d0, idim^LIM,global_time+0.5d0*dt,ps1,global_time,ps2)
          call advect1(flux_scheme,1.0d0, idim^LIM,global_time+0.5d0*dt,ps2,global_time,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=(1.0d0/3.0d0)*(-ps(igrid)%cons+ps1(igrid)%cons+2.0d0*ps2(igrid)%cons+ps3(igrid)%cons)
             if(stagger_grid) ps(igrid)%prims=(1.0d0/3.0d0) &
                 *(-ps(igrid)%prims+ps1(igrid)%prims+2.0d0*ps2(igrid)%prims+ps3(igrid)%prims)
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,1.0d0/6.0d0, idim^LIM,global_time+dt,ps3,global_time+dt*5.0d0/6.0d0,ps)

       case ("jameson")
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=ps(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims=ps(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,1.0d0/4.0d0, idim^LIM,global_time,ps,global_time,ps1)
          call advect1(flux_scheme,1.0d0/3.0d0, idim^LIM,global_time+dt/4.0d0,ps1,global_time,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons=ps(igrid)%cons
             if(stagger_grid) ps1(igrid)%prims=ps(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,0.5d0, idim^LIM,global_time+dt/3.0d0,ps2,global_time,ps1)
          call advect1(flux_scheme,1.0d0, idim^LIM,global_time+0.5d0*dt,ps1,global_time,ps)

       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          write(unitterm,*) "Error in advect: Unknown time integration method"
          call mpistop("Correct time_integrator")
       end select

    case ("fivestep")
       select case (time_integrator)
       case ("ssprk5")
          ! SSPRK(5,4) by Ruuth and Spiteri
          call advect1(flux_scheme,rk_beta11, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=rk_alfa21*ps(igrid)%cons+rk_alfa22*ps1(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims=rk_alfa21*ps(igrid)%prims+rk_alfa22*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,global_time+rk_alfa22*rk_c2*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%cons=rk_alfa31*ps(igrid)%cons+rk_alfa33*ps2(igrid)%cons
             if(stagger_grid) ps1(igrid)%prims=rk_alfa31*ps(igrid)%prims+rk_alfa33*ps2(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta33, idim^LIM,global_time+rk_c3*dt,ps2,global_time+rk_alfa33*rk_c3*dt,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%cons=rk_alfa53*ps2(igrid)%cons+rk_alfa54*ps1(igrid)%cons
             if(stagger_grid) ps3(igrid)%prims=rk_alfa53*ps2(igrid)%prims+rk_alfa54*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%cons=rk_alfa41*ps(igrid)%cons+rk_alfa44*ps1(igrid)%cons
             if(stagger_grid) ps2(igrid)%prims=rk_alfa41*ps(igrid)%prims+rk_alfa44*ps1(igrid)%prims
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta44, idim^LIM,global_time+rk_c4*dt,ps1,global_time+rk_alfa44*rk_c4*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%cons=ps3(igrid)%cons+rk_alfa55*ps2(igrid)%cons &
                +(rk_beta54/rk_beta44)*(ps2(igrid)%cons-(rk_alfa41*ps(igrid)%cons+rk_alfa44*ps1(igrid)%cons))
             if(stagger_grid) then
             ps(igrid)%prims=ps3(igrid)%prims+rk_alfa55*ps2(igrid)%prims &
                +(rk_beta54/rk_beta44)*(ps2(igrid)%prims-(rk_alfa41*ps(igrid)%prims+rk_alfa44*ps1(igrid)%prims))
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_scheme,rk_beta55, idim^LIM,global_time+rk_c5*dt,ps2,global_time+(1.0d0-rk_beta55)*dt,ps)

       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          write(unitterm,*) "Error in advect: Unknown time integration method"
          call mpistop("Correct time_integrator")
       end select

    case default
       write(unitterm,*) "time_stepper=",time_stepper
       write(unitterm,*) "Error in advect: Unknown time stepping method"
       call mpistop("Correct time_stepper")
    end select

    firstsweep=.false.
  end subroutine advect

  !> Implicit global update step within IMEX schemes, advance psa=psb+dtfactor*qdt*F_im(psa)
  subroutine global_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_implicit_update, phys_req_diagonal

    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks)   !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC      !< Both states psa and psb at this time level
    double precision, intent(in) :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)
    integer :: iigrid, igrid
   
    ! note that here is treat implicit step as a step
    istep = istep+1

    if (associated(phys_implicit_update)) then
       call phys_implicit_update(dtfactor,qdt,qtC,psa,psb)
    else
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          psa(igrid)%cons = psb(igrid)%cons
       end do
       !$OMP END PARALLEL DO
    end if
    
    ! enforce boundary conditions for psa
    !call getbc(qtC,0.d0,psa,nhydro_lo,nhydro,phys_req_diagonal)

  end subroutine global_implicit_update

  !> Evaluate Implicit part in place, i.e. psa==>F_im(psa)
  subroutine evaluate_implicit(qtC,psa)
    use mod_global_parameters
    use mod_physics, only: phys_evaluate_implicit

    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    double precision, intent(in) :: qtC      !< psa at this time level
    integer :: iigrid, igrid

    if (associated(phys_evaluate_implicit)) then
       call phys_evaluate_implicit(qtC,psa)
    else
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          psa(igrid)%cons = 0.0d0
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine evaluate_implicit

  !> Integrate all grids by one partial step
  subroutine advect1(method,dtfactor,idim^LIM,qtC,psa,qt,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics

    integer, intent(in) :: idim^LIM
    type(state), target :: psa(max_blocks) !< Compute fluxes based on this state
    type(state), target :: psb(max_blocks) !< Update solution on this state
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    character(len=*), intent(in) :: method(nlevelshi)

    double precision :: qdt
    integer :: iigrid, igrid, level

    logical :: setigrid

    ! For the first step, it is not necessarry to do the con2prim and getbc
    if (istep > 0) then
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          call phys_to_primitive(ixG^LL,ixM^LL,psa(igrid)%cons(ixG^T,1:ncons),&
                        psa(igrid)%prim(ixG^T,1:nprim),psa(igrid)%x(ixG^T,1:ndim))
       end do
       !$OMP END PARALLEL DO
       ! fixme
       !call getbc(qt,qdt,psa,nhydro_lo,nhydro,phys_req_diagonal)
       call getbc(qt,qdt,psa,1,nprim,phys_req_diagonal)
    end if

    istep = istep+1

    ! fixme: what is this?
    if(associated(phys_special_advance)) then
      call phys_special_advance(qdt,qt,psa)
    end if

    ! opedit: Just advance the active grids:
    !$OMP PARALLEL DO PRIVATE(igrid,level,qdt)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       level=node(plevel_,igrid)
       qdt=dtfactor*dt_grid(igrid)

       call process1_grid(method(level),igrid,qdt,ixG^LL,idim^LIM,qtC,&
            psa(igrid),qt,psb(igrid),pso(igrid))

    end do
    !$OMP END PARALLEL DO

    ! opedit: Send flux for all grids, expects sends for all
    ! nsend_fc(^D), set in connectivity.t.

    if (fix_conserve_at_step) then
      call recvflux(idim^LIM)
      call sendflux(idim^LIM)
      call fix_conserve(psb,idim^LIM,1,ncons)
      if(stagger_grid) call fix_edges(psb,idim^LIM)
    end if

    if(stagger_grid) then
      ! Now fill the cell-center values for the staggered variables
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        call phys_face_to_center(ixG^LL,psb(igrid))
      end do
      !$OMP END PARALLEL DO
    end if

    ! fixme: what is this?
    qdt = dtfactor*dt

  end subroutine advect1

  !> Prepare to advance a single grid over one partial time step
  subroutine process1_grid(method,igrid,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold)
    use mod_global_parameters
    use mod_fix_conserve

    character(len=*), intent(in) :: method
    integer, intent(in) :: igrid, ixI^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt
    type(state), target          :: sCT, s, sold

    double precision :: dx^D
    ! cell face flux
    double precision :: fC(ixI^S,1:ncons,1:ndim)
    ! cell edge flux
    double precision :: fE(ixI^S,7-2*ndim:3)

    dx^D=rnode(rpdx^D_,igrid);
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    saveigrid=igrid

    block0=>sCT
    block=>s
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    call advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D, &
         ps(igrid)%x)


    ! opedit: Obviously, flux is stored only for active grids.
    ! but we know in fix_conserve wether there is a passive neighbor
    ! via neighbor_active(i^D,igrid) thus we skip the correction for those.
    ! This violates strict conservation when the active/passive interface
    ! coincides with a coarse/fine interface.
    if (fix_conserve_at_step) then
      call store_flux(igrid,fC,idim^LIM,ncons)
      if(stagger_grid) call store_edge(igrid,ixI^L,fE,idim^LIM)
    end if

  end subroutine process1_grid

  !> Advance a single grid over one partial time step
  subroutine advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D,x)

    !  integrate one grid by one partial step
    use mod_finite_volume
    use mod_tvd
    use mod_source, only: addsource2
    use mod_global_parameters

    character(len=*), intent(in) :: method
    integer, intent(in) :: ixI^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    type(state), target          :: sCT, s, sold
    double precision :: fC(ixI^S,1:ncons,1:ndim)
    double precision :: fE(ixI^S,7-2*ndim:3)

    integer :: ixO^L

    ixO^L=ixI^L^LSUBnghostcells;

    select case (method)
    case ('tvdmu','tvdlf','hll')
       call finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D,x)
    case ('tvd')
       call tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,sCT,qt+qdt,s,fC,dx^D,x)
    case ('source')
       call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),&
            ixI^L,ixO^L,1,nprim,qtC,sCT%prim,qt,s%prim,x,.false.)
    case ('nul')
       ! There is nothing to do
    case default
       write(unitterm,*)'Error in advect1_grid:',method,' is unknown!'
       call mpistop("")
    end select

  end subroutine advect1_grid

  !> process is a user entry in time loop, before output and advance
  !>         allows to modify solution, add extra variables, etc.
  !> Warning: CFL dt already determined (and is not recomputed)! 
  subroutine process(iit,qt)
    use mod_usr_methods, only: usr_process_grid, usr_process_global
    use mod_physics, only: phys_process_grid, &
                           phys_process_global
    use mod_global_parameters
    use mod_ghostcells_update
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt
    logical :: update_bc
    integer:: iigrid, igrid,level

    update_bc = .False.

    if (associated(phys_process_global)) then
       call phys_process_global(iit,qt)
       update_bc = .True.
    end if

    if (associated(phys_process_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid,level)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         level=node(plevel_,igrid)
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         typelimiter=type_limiter(node(plevel_,igrid))
         typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

         call phys_process_grid(igrid,level,ixG^LL,ixM^LL, &
              qt,ps(igrid)%prim,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      update_bc = .True.
    end if
    if (update_bc) call getbc(qt,dt,ps,1,nprim)

    if (associated(usr_process_global)) then
       call usr_process_global(iit,qt)
       update_bc = .True.
    end if

    if (associated(usr_process_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid,level)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         level=node(plevel_,igrid)
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         typelimiter=type_limiter(node(plevel_,igrid))
         typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
         call usr_process_grid(igrid,level,ixG^LL,ixM^LL,qt,ps(igrid)%prim,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
       update_bc = .True.
    end if
    if (update_bc) call getbc(qt,dt,ps,1,nprim)

  end subroutine process

  !> process_advanced is user entry in time loop, just after advance
  !>           allows to modify solution, add extra variables, etc.
  !>           added for handling two-way coupled PIC-MHD
  !> Warning: w is now at global_time^(n+1), global time and iteration at global_time^n, it^n
  subroutine process_advanced(iit,qt)
    use mod_small_values, only: prim_NaN_checker
    use mod_usr_methods, only: usr_process_adv_grid, &
                               usr_process_adv_global
    use mod_physics, only: phys_handle_small_values, &
                           phys_process_adv_grid, &
                           phys_process_adv_global
    use mod_global_parameters
    use mod_ghostcells_update
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt
    logical :: update_bc
    integer:: iigrid, igrid,level

    update_bc = .False.

    if (associated(phys_process_adv_global)) then
       call phys_process_adv_global(iit,qt)
       update_bc = .True.
    end if

    if (associated(phys_process_adv_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid,level)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         level=node(plevel_,igrid)
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         typelimiter=type_limiter(node(plevel_,igrid))
         typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

         call phys_process_adv_grid(igrid,level,ixG^LL,ixM^LL, &
              qt,ps(igrid)%prim,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      update_bc = .True.
    end if
    if (update_bc) call getbc(qt,dt,ps,1,nprim)

    if (associated(usr_process_adv_global)) then
       call usr_process_adv_global(iit,qt)
       update_bc = .True.
    end if

    if (associated(usr_process_adv_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid,level)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         level=node(plevel_,igrid)
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         typelimiter=type_limiter(node(plevel_,igrid))
         typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

         call usr_process_adv_grid(igrid,level,ixG^LL,ixM^LL, &
              qt,ps(igrid)%prim,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      update_bc = .True.
    end if

    if (update_bc) call getbc(qt,dt,ps,1,nprim)

    ! check and optionally correct unphysical values
    ! note that NaN of rho/ eps will be treat as small value here.
    if ( fix_small_values ) then
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          call phys_handle_small_values(ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim),ixG^LL,ixM^LL, .True., 'process_advanced')
       end do
       !$OMP END PARALLEL DO
    end if

    ! check if there are any NaN_values
    if (check_NaN_values) then
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          call prim_NaN_checker(ixG^LL,ixM^LL,ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine process_advanced

end module mod_advance
