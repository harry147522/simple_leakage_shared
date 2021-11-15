module mod_gremhd_phys_implicit_update
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private
  
  logical :: debug_flag = .False.
  !logical :: debug_flag = .True.

  ! Public methods
  public :: gremhd_phys_implicit_update_init

contains

  !> Initialize the module
  subroutine gremhd_phys_implicit_update_init()
    phys_implicit_update     => gremhd_implicit_update
    phys_evaluate_implicit   => gremhd_evaluate_implicit
  end subroutine gremhd_phys_implicit_update_init

  !> Evaluate Implicit part in place, i.e. psa==>F_im(psa)
  subroutine gremhd_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks) !< Compute implicit part from this state and update it
    double precision, intent(in) :: qtC !< Both states psa and psb at this time level

    integer :: iigrid, igrid

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! primitive variables W_vel are needed in stiff source term
       call phys_to_primitive(ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
          psa(igrid)%cons(ixGlo1:ixGhi1,1:ncons),psa(igrid)%prim(ixGlo1:ixGhi1,&
          1:nprim),psa(igrid)%x(ixGlo1:ixGhi1,1:ndim))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call gremhd_get_stiff_source(qtC,ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
           psa(igrid)%cons(ixGlo1:ixGhi1,1:ncons),&
          psa(igrid)%prim(ixGlo1:ixGhi1,1:nprim),psa(igrid)%x(ixGlo1:ixGhi1,&
          1:ndim))
    end do
    !$OMP END PARALLEL DO

  end subroutine gremhd_evaluate_implicit

  !> Implicit solve of psa=psb+dtfactor*dt*F_im(psa)
  subroutine gremhd_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    type(state), target :: psa(max_blocks) !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks)   !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC !< Both states psa and psb at this time level
    double precision, intent(in) :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)

    integer :: iigrid, igrid

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%cons = psb(igrid)%cons
       call gremhd_implicit_update_grid(dtfactor,qdt,qtC,ixGlo1,ixGhi1,ixMlo1,&
          ixMhi1, psa(igrid)%cons(ixGlo1:ixGhi1,1:ncons),&
          psa(igrid)%prim(ixGlo1:ixGhi1,1:nprim),psa(igrid)%x(ixGlo1:ixGhi1,&
          1:ndim))
    end do
    !$OMP END PARALLEL DO
  end subroutine gremhd_implicit_update

  subroutine gremhd_get_stiff_source(qtC,ixImin1,ixImax1, ixOmin1,ixOmax1,&
      cons, prim, x)
    use mod_global_parameters
    use mod_usr_methods
    double precision, intent(in)    :: qtC !< Both states psa and psb at this time level
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: cons(ixImin1:ixImax1, 1:ncons)
    double precision, intent(in)    :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision                :: lfac(ixImin1:ixImax1) ! Lorentz factor
    double precision                :: psi6(ixImin1:ixImax1) !conformal factor **6
    double precision                :: eta(ixImin1:ixImax1) ! resistivity
    double precision                :: xi(ixImin1:ixImax1) !proper dynamo action term
    double precision                :: qE_dot_v(ixImin1:ixImax1),&
        qB_dot_v(ixImin1:ixImax1) 
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
    double precision                :: sqrt_gamma(ixImin1:ixImax1)
    double precision                :: qJi(ixImin1:ixImax1,1:3) !conformally transformed Current, but only the stiff part
    integer                         :: idir, jdir, kdir

    call usr_get_resistivity(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        cons(ixImin1:ixImax1, 1:ncons), x(ixImin1:ixImax1, 1:ndim),&
        eta(ixImin1:ixImax1))
    if (associated(usr_get_dynamo_coeff)) then
       call usr_get_dynamo_coeff(ixImin1,ixImax1, ixOmin1,ixOmax1,&
           cons(ixImin1:ixImax1, 1:ncons), x(ixImin1:ixImax1, 1:ndim),&
           xi(ixImin1:ixImax1))
    else
       xi(ixOmin1:ixOmax1) = 0.0d0 ! no dynamo at the moment
    end if
    call gremhd_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        prim(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3),&
        sqrt_gamma=sqrt_gamma(ixImin1:ixImax1), lfac=lfac(ixImin1:ixImax1) )

    qE_dot_v(ixOmin1:ixOmax1) = 0.0d0
    do idir = 1, ndir
       qE_dot_v(ixOmin1:ixOmax1) = qE_dot_v(ixOmin1:ixOmax1) + &
          cons(ixOmin1:ixOmax1, Econs(idir)) * prim(ixOmin1:ixOmax1,&
           W_vel(idir)) * gamma(ixOmin1:ixOmax1, idir, idir)
    end do
    qE_dot_v(ixOmin1:ixOmax1) = qE_dot_v(ixOmin1:ixOmax1) / &
       lfac(ixOmin1:ixOmax1)
    qB_dot_v(ixOmin1:ixOmax1) = 0.0d0
    do idir = 1, ndir
       qB_dot_v(ixOmin1:ixOmax1) = qB_dot_v(ixOmin1:ixOmax1) + &
          cons(ixOmin1:ixOmax1, Bcons(idir)) * prim(ixOmin1:ixOmax1,&
           W_vel(idir)) * gamma(ixOmin1:ixOmax1, idir, idir)
    end do
    qB_dot_v(ixOmin1:ixOmax1) = qB_dot_v(ixOmin1:ixOmax1) / &
       lfac(ixOmin1:ixOmax1)

    ! dynamo part
    do idir = 1, ndir
       qJi(ixOmin1:ixOmax1, idir) = lfac(ixOmin1:ixOmax1) * &
          cons(ixOmin1:ixOmax1, Bcons(idir)) - qB_dot_v(ixOmin1:ixOmax1) * &
          prim(ixOmin1:ixOmax1, W_vel(idir))
       do jdir = 1, ndir
          do kdir = 1, ndir
             qJi(ixOmin1:ixOmax1, idir) = qJi(ixOmin1:ixOmax1,&
                 idir) - lvc(idir, jdir, kdir) / gamma(ixOmin1:ixOmax1, idir,&
                 idir) * sqrt_gamma(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
                 W_vel(jdir)) * cons(ixOmin1:ixOmax1, Econs(kdir))
          end do
       end do
    end do
    do idir = 1, ndir
       qJi(ixOmin1:ixOmax1, idir) = -xi(ixOmin1:ixOmax1) * qJi(ixOmin1:ixOmax1,&
           idir)
    end do

    do idir = 1, ndir
       qJi(ixOmin1:ixOmax1, idir) = qJi(ixOmin1:ixOmax1,&
           idir) + lfac(ixOmin1:ixOmax1) * cons(ixOmin1:ixOmax1,&
           Econs(idir)) - qE_dot_v(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
           W_vel(idir))
       do jdir = 1, ndir
          do kdir = 1, ndir
             qJi(ixOmin1:ixOmax1, idir) = qJi(ixOmin1:ixOmax1,&
                 idir) + lvc(idir, jdir, kdir) / gamma(ixOmin1:ixOmax1, idir,&
                 idir) * sqrt_gamma(ixOmin1:ixOmax1) * prim(ixOmin1:ixOmax1,&
                 W_vel(jdir)) * cons(ixOmin1:ixOmax1, Bcons(kdir))
          end do
       end do
    end do

    ! finally, we store the source terms directly in the cons
    cons(ixOmin1:ixOmax1, 1:ncons) = 0.0d0
    do idir = 1, ndir
       cons(ixOmin1:ixOmax1, Econs(idir)) = - prim(ixOmin1:ixOmax1,&
           alp_) / eta(ixOmin1:ixOmax1) * qJi(ixOmin1:ixOmax1, idir)
    end do
  end subroutine gremhd_get_stiff_source

  subroutine gremhd_implicit_update_grid(dtfactor,qdt,qtC,ixImin1,ixImax1,&
      ixOmin1,ixOmax1, cons, prim, x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    use mod_eos
    use mod_eos_idealgas, only: eos_gamma
    use mod_rootfinding
    double precision, intent(in)    :: qdt      !< overall time step dt
    double precision, intent(in)    :: qtC !< Both states psa and psb at this time level
    double precision, intent(in)    :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: cons(ixImin1:ixImax1, 1:ncons)
    double precision, intent(inout) :: prim(ixImin1:ixImax1, 1:nprim)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision                :: one_over_eos_gamma_1 = 0.25d0 !< (eos_gamma-1)/eos_gamma  
    double precision                :: sigma_L = 0.0d0 !< this is defined as alp * dt
    double precision                :: eta_tilde = 0.0d0 !< this is defined as eta / ( alp * dt )

    integer                         :: idir, jdir, kdir
    double precision                :: psi6(ixImin1:ixImax1)
    integer                         :: ix1 
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3) !metric gamma_ij
    double precision                :: eta(ixImin1:ixImax1) ! resistivity
    double precision                :: xi(ixImin1:ixImax1) !proper dynamo action term
    double precision                :: sqrt_gamma(ixImin1:ixImax1)

    ! variables needed for the root founding
    logical                         :: adjustment
    integer                         :: error_code
    double precision                :: cons_tmp(ixImin1:ixImax1, 1:ncons)
    double precision                :: Ei(1:ndir)
    double precision                :: ui_hat(1:ndir)
    double precision                :: Jac(1:ndir, 1:ndir)
    double precision                :: lfac
   
    one_over_eos_gamma_1 = (eos_gamma - 1.0d0) / eos_gamma
    !one_over_eos_gamma_1 = (atmo_gamma - 1.0d0) / atmo_gamma

    call usr_get_resistivity(ixImin1,ixImax1, ixOmin1,ixOmax1,&
        cons(ixImin1:ixImax1, 1:ncons), x(ixImin1:ixImax1, 1:ndim),&
        eta(ixImin1:ixImax1))
    if(associated(usr_get_dynamo_coeff)) then
       call usr_get_dynamo_coeff(ixImin1,ixImax1, ixOmin1,ixOmax1,&
           cons(ixImin1:ixImax1, 1:ncons), x(ixImin1:ixImax1, 1:ndim),&
           xi(ixImin1:ixImax1))
    else
       xi(ixOmin1:ixOmax1) = 0.0d0 ! no dynamo at the moment
    end if

    call get_gamma_ij_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1, ixOmin1,&
       ixOmax1, gamma(ixImin1:ixImax1,1:3,1:3))
    do idir = 1, ndir
       gamma(ixOmin1:ixOmax1,idir,idir) = gamma(ixOmin1:ixOmax1,idir,&
          idir) * prim(ixOmin1:ixOmax1, psi_)**4 
    end do
    psi6(ixOmin1:ixOmax1) = prim(ixOmin1:ixOmax1, psi_)**6 
    call get_sqrt_gamma_hat(x(ixImin1:ixImax1, 1:ndim), ixImin1,ixImax1,&
        ixOmin1,ixOmax1, sqrt_gamma(ixImin1:ixImax1))
    sqrt_gamma(ixOmin1:ixOmax1) = psi6(ixOmin1:ixOmax1) * &
       sqrt_gamma(ixOmin1:ixOmax1)
 
    ! conformal transformation back to normal conserved variables
    cons_tmp(ixOmin1:ixOmax1, D_) = cons(ixOmin1:ixOmax1,&
        D_) / psi6(ixOmin1:ixOmax1)
    cons_tmp(ixOmin1:ixOmax1, tau_) = cons(ixOmin1:ixOmax1,&
        tau_) / psi6(ixOmin1:ixOmax1)
    do idir = 1, ndir
       cons_tmp(ixOmin1:ixOmax1, mom(idir)) = cons(ixOmin1:ixOmax1,&
           mom(idir)) / psi6(ixOmin1:ixOmax1)
       cons_tmp(ixOmin1:ixOmax1, Econs(idir)) = cons(ixOmin1:ixOmax1,&
           Econs(idir)) / psi6(ixOmin1:ixOmax1)
       cons_tmp(ixOmin1:ixOmax1, Bcons(idir)) = cons(ixOmin1:ixOmax1,&
           Bcons(idir)) / psi6(ixOmin1:ixOmax1)
    end do

    do ix1 = ixOmin1,ixOmax1 
      ! get sigma_L
      sigma_L = dtfactor * qdt * prim(ix1, alp_) ! sigma_L = alpha * delta_t
      eta_tilde = eta(ix1) / sigma_L
      if ( (cons_tmp(ix1, D_) <= small_D) .or. (cons_tmp(ix1,&
          tau_) <= small_tau) ) then
         ! note: EM part is contained in this tau, however this relation should be held as well.
         ! atmosphere handling 
         ! skip the primitive recovery
         ui_hat(1:ndir) = 0.0d0
      else
         ! initial guess of ui
         do idir = 1, ndir
            ui_hat(idir) = prim(ix1, W_vel(idir))
         end do
         call rootfinding_global_multid_newton_raphson(ndir, ui_hat, tol_im,&
             iter_max_im, error_code, func_of_u, jacobian_of_f, fmin_of_f)
         ! check the solution
         select case (error_code)
         !case (0) ! root is found
         !   write(*,*) "z= ", z
         case (-1) ! nothing happened
            call mpistop("have you ever attemp to find the root in con2prim?")
         case (1) ! failed to find the root
            if (debug_flag) then
               write(*,*) 'location x =',x(ix1, 1:ndim)
               write(*,*) 'rho =', prim(ix1, rho_), small_rho
               write(*,*) 'D = ', cons(ix1, D_), cons_tmp(ix1, D_),&
                   small_rho_thr
               write(*,*) 'D_atmo/D = ', small_D/cons(ix1, D_)
               write(*,*) 'Wv = ', prim(ix1, W_vel(:))
               write(*,*) 'ui = ', ui_hat(1:ndir)
               write(*,*) 'f(Wv) = ', func_of_u(prim(ix1, W_vel(:)))
               write(*,*) 'f = ', func_of_u(ui_hat)
               write(*,*) 'eta =',eta(ix1), eta_tilde
               write(*,*) 'sigma =',1.0d0/eta(ix1)
               write(*,*) 'cons = ', cons(ix1, 1:ncons)
               call mpistop&
                  ("Error: Failed to find the root in gremhd_implicit_update")
            end if
            ! treat it as atmosphere
            ui_hat = 0.0d0
         case (2) ! z is NaN
            call mpistop("Error: the root is NaN in gremhd_implicit_update")
         case (3) ! reached local minimum
            if (debug_flag) then
               write(*,*) 'location x =',x(ix1, 1:ndim)
               write(*,*) 'rho =', prim(ix1, rho_), small_rho
               write(*,*) 'D = ', cons(ix1, D_), cons_tmp(ix1, D_),&
                   small_rho_thr
               write(*,*) 'D_atmo/D = ', small_D/cons(ix1, D_)
               write(*,*) 'Wv = ', prim(ix1, W_vel(:))
               write(*,*) 'ui = ', ui_hat(1:ndir)
               write(*,*) 'f(Wv) = ', func_of_u(prim(ix1, W_vel(:)))
               write(*,*) 'f = ', func_of_u(ui_hat)
               write(*,*) 'eta =',eta(ix1), eta_tilde
               write(*,*) 'sigma =',1.0d0/eta(ix1)
               write(*,*) 'cons = ', cons(ix1, 1:ncons)
               call mpistop("local minimum in gremhd_implicit_update")
            end if
         case (4) ! singular jacobian
            call mpistop("Error: singular jacobian in gremhd_implicit_update")
         case (5) ! slope is larger or equals to zero
            call mpistop&
("Error: slope is larger or equals to zero in gremhd_implicit_update")
         end select
      end if

      ! update conserved E field from the roots
      call get_vars(ui_hat, W_out=lfac, E_out=Ei)
      ! check if the resulting rho is lower then atmo
      if ( cons_tmp(ix1, D_) / lfac <= small_rho_thr ) then
         ui_hat = 0.0d0 ! reset ui
         ! recalculate E field
         call get_vars(ui_hat, E_out=Ei)
      end if

      do idir = 1, ndir
         cons_tmp(ix1, Econs(idir)) = Ei(idir)
         ! save the solution into prim, as an initial guess for other con2prim if needed
         prim(ix1, W_vel(idir)) = ui_hat(idir)
      end do
    enddo

    do idir = 1, ndir
       cons(ixOmin1:ixOmax1, Econs(idir)) = cons_tmp(ixOmin1:ixOmax1,&
           Econs(idir)) * psi6(ixOmin1:ixOmax1) 
    end do

    contains

       !> func f(ui) for finding the root
       subroutine get_vars(ui_in, W_out, E_out, f_out, fjac_out)
          implicit none
          double precision, intent(in)    :: ui_in(1:ndir)
          double precision, intent(out), optional   :: W_out
          double precision, intent(out), optional   :: E_out(1:ndir)
          double precision, intent(out), optional   :: f_out(1:ndir)
          double precision, intent(out), optional   :: fjac_out(1:ndir,&
              1:ndir)

          double precision                :: ui(1:ndir)

          double precision                :: E_dot_u, B_dot_u
          double precision                :: coeff_A(0:5)
          double precision                :: coeff_A_dot(0:5)

          integer                         :: idir,jdir,kdir,adir,bdir
          double precision                :: dhDdu(1:ndir)
          double precision                :: E_field(1:ndir)
          double precision                :: S_i_prime(1:ndir)
          double precision                :: dEdu(1:ndir,1:ndir)
          double precision                :: fluid_energy
          double precision                :: W_hat, hD
          
          W_hat = 1.0d0
          do idir = 1, ndir
             W_hat = W_hat + ui_in(idir)**2 * gamma(ix1, idir, idir)
          end do
          W_hat = dsqrt( W_hat )
          if ( present(W_out) ) W_out = W_hat
    
          ui(:) = ui_in(:)

          ! get E field from u
          coeff_A(0) = W_hat + eta_tilde + xi(ix1)**2 * ( W_hat**2 - 1.0d0 ) / &
             ( W_hat + eta_tilde ) 
          coeff_A(1) = ( 1.0d0 + xi(ix1)**2 * W_hat / ( W_hat + eta_tilde ) ) &
             * eta_tilde / ( 1.0d0 + W_hat * eta_tilde )
          coeff_A(2) = - xi(ix1) * eta_tilde / ( W_hat + eta_tilde )
          coeff_A(3) = xi(ix1) * ( 1.0d0 + eta_tilde * W_hat ) / ( W_hat + &
             eta_tilde )
          coeff_A(4) = xi(ix1) * ( 1.0d0 - eta_tilde**2 + xi(ix1)**2 ) / ( &
             W_hat + eta_tilde ) / ( 1.0d0 + W_hat * eta_tilde )
          coeff_A(5) = - 1.0d0 - ( xi(ix1)**2 * W_hat ) / ( W_hat + eta_tilde &
             )

          E_dot_u = 0.0d0
          B_dot_u = 0.0d0
          do idir = 1, ndir
             E_dot_u = E_dot_u + cons_tmp(ix1,&
                 Econs(idir)) * ui(idir) * gamma(ix1, idir, idir)
             B_dot_u = B_dot_u + cons_tmp(ix1,&
                 Bcons(idir)) * ui(idir) * gamma(ix1, idir, idir)
          end do
          
          do idir = 1, ndir
             E_field(idir) = eta_tilde * cons_tmp(ix1,&
                 Econs(idir)) + coeff_A(1) * E_dot_u * ui(idir) + coeff_A(3) * &
                cons_tmp(ix1, Bcons(idir)) + coeff_A(4) * B_dot_u * ui(idir) 
             do jdir = 1, ndir
                do kdir = 1, ndir
                   E_field(idir) = E_field(idir) + lvc(idir,jdir,&
                      kdir) / sqrt_gamma(ix1) * ui(jdir) * gamma(ix1, jdir,&
                       jdir) * gamma(ix1, kdir,&
                       kdir) * ( coeff_A(2) * cons_tmp(ix1,&
                       Econs(kdir)) + coeff_A(5) * cons_tmp(ix1,&
                       Bcons(kdir)) ) 
                end do
             end do
          end do
          E_field = E_field / coeff_A(0)
          if (present(E_out)) E_out = E_field

          if (present(f_out).or.present(fjac_out)) then
             ! fluid_energy
             fluid_energy = 0.0d0
             do idir = 1, ndir
                fluid_energy = fluid_energy + gamma(ix1, idir,&
                    idir) * ( E_field(idir)**2 + cons_tmp(ix1,&
                    Bcons(idir))**2 )
             end do
             fluid_energy = cons_tmp(ix1,&
                 tau_) - 0.5d0 * fluid_energy + cons_tmp(ix1, D_)
   
             ! modified enthalpy
             hD = ( W_hat * fluid_energy - one_over_eos_gamma_1 * cons_tmp(ix1,&
                 D_) ) / ( W_hat**2 - one_over_eos_gamma_1 ) 

             ! S_i_prime = S_i - EM part
             do idir = 1, ndir
                S_i_prime(idir) = cons_tmp(ix1, mom(idir))
                do jdir = 1, ndir
                   do kdir = 1, ndir
                      S_i_prime(idir) = S_i_prime(idir) - sqrt_gamma(ix1) * &
                         lvc(idir,jdir,kdir) * E_field(jdir) * cons_tmp(ix1,&
                          Bcons(kdir))
                   end do
                end do
             end do
   
             if (present(f_out)) then
                do idir = 1, ndir
                   f_out(idir) = gamma(ix1, idir,&
                       idir) * ui(idir) - S_i_prime(idir) / hD
                end do
             end if
             
             if (present(fjac_out)) then 
                ! get dEdu
                coeff_A_dot(0) = 1.0d0 + xi(ix1)**2 * ( 1.0d0 + W_hat**2 + &
                   2.0d0 * eta_tilde * W_hat ) / ( W_hat + eta_tilde )**2
                coeff_A_dot(1) = - eta_tilde**2 / ( 1.0d0 + W_hat * eta_tilde &
                   )**2 * ( 1.0d0 + xi(ix1)**2 * ( W_hat**2 - 1.0d0 )**2 / ( &
                   W_hat + eta_tilde )**2 )
                coeff_A_dot(2) = xi(ix1) * eta_tilde / ( W_hat + eta_tilde &
                   )**2
                coeff_A_dot(3) = xi(ix1) * ( eta_tilde**2 - 1.0d0 ) / ( W_hat &
                   + eta_tilde )**2
                coeff_A_dot(4) = - xi(ix1) * ( 1.0d0 - eta_tilde**2 + &
                   xi(ix1)**2 ) * ( 1.0d0 + eta_tilde**2 + 2.0d0 * eta_tilde * &
                   W_hat) / ( W_hat + eta_tilde )**2 / ( 1.0d0 + W_hat * &
                   eta_tilde )**2 
                coeff_A_dot(5) = - xi(ix1)**2 * eta_tilde / ( W_hat + &
                   eta_tilde )**2
   
                dEdu(:,:) = 0.0d0
                do adir = 1, ndir
                   do bdir = 1, ndir
      
                      dEdu(adir, bdir) = ( - coeff_A_dot(0) * E_field(adir) + &
                         coeff_A_dot(3) * cons_tmp(ix1,&
                          Bcons(adir)) ) * ui(bdir) * gamma(ix1, bdir,&
                          bdir) + W_hat * ui(adir) * ( coeff_A(1) * &
                         cons_tmp(ix1, Econs(bdir)) + coeff_A(4) * &
                         cons_tmp(ix1, Bcons(bdir)) ) * gamma(ix1, bdir,&
                          bdir) + E_dot_u * ( coeff_A(1) * W_hat * kr(adir,&
                         bdir) + coeff_A_dot(1) * ui(adir) * ui(bdir) * &
                         gamma(ix1, bdir,&
                          bdir) ) + B_dot_u * ( coeff_A(4) * W_hat * kr(adir,&
                         bdir) + coeff_A_dot(4) * ui(adir) * ui(bdir) * &
                         gamma(ix1, bdir, bdir) )
                      do jdir = 1, ndir
                         do kdir = 1, ndir
                            dEdu(adir, bdir) = dEdu(adir, bdir) + lvc(adir,&
                               jdir,kdir) / sqrt_gamma(ix1) * gamma(ix1, kdir,&
                                kdir)  * ( cons_tmp(ix1,&
                                Bcons(kdir)) * ( coeff_A(5) * W_hat * &
                               gamma(ix1, jdir,&
                                bdir) + coeff_A_dot(5) * ui(jdir) * ui(adir) * &
                               gamma(ix1, jdir, jdir) * gamma(ix1, adir,&
                                adir)) + cons_tmp(ix1,&
                                Econs(kdir)) * ( coeff_A(2) * W_hat * &
                               gamma(ix1, jdir,&
                                bdir) + coeff_A_dot(2) * ui(jdir) * ui(adir) * &
                               gamma(ix1, jdir, jdir) * gamma(ix1, adir,&
                                adir)) )
                         end do
                      end do
      
                   end do
                end do
                dEdu(:,:) = dEdu(:,:) / coeff_A(0) / W_hat
       
                do jdir = 1, ndir
                   dhDdu(jdir) = ( fluid_energy / W_hat - 2.0d0 * hD ) * &
                      ui(jdir) * gamma(ix1, jdir, jdir)
                   do idir = 1, ndir
                      dhDdu(jdir) = dhDdu(jdir) - W_hat * E_field(idir) * &
                         gamma(ix1, idir, idir) * dEdu(idir, jdir)
                   end do
                end do
                dhDdu = dhDdu / ( W_hat**2 - one_over_eos_gamma_1 )
   
                ! finally, get Jac
                do idir = 1, ndir
                   do jdir = 1, ndir
   
                      fjac_out(idir, jdir) = gamma(ix1, idir,&
                          jdir) + S_i_prime(idir) / hD**2 * dhDdu(jdir)

                      do adir = 1, ndir
                         do bdir = 1, ndir
                            fjac_out(idir, jdir) = fjac_out(idir,&
                                jdir) + sqrt_gamma(ix1) / hD * lvc(idir,adir,&
                               bdir) * dEdu(adir, jdir) * cons_tmp(ix1,&
                                Bcons(bdir))
                         end do
                      end do

                   end do
                end do
             end if
          end if
       end subroutine get_vars

       !> master function f(ui) for finding the root
       function func_of_u(ui) result( f_of_u )
          implicit none
          double precision, dimension(:), intent(in)    :: ui
          double precision, dimension(size(ui))         :: f_of_u
          call get_vars(ui, f_out=f_of_u)
       end function func_of_u

       !> Jacobian of f
       function jacobian_of_f(ui) result( Jac )
          implicit none
          double precision, dimension(:), intent(in)     :: ui
          double precision, dimension(size(ui), size(ui)):: Jac
          call get_vars(ui, fjac_out=Jac(:,:))
       end function jacobian_of_f

       function fmin_of_f(ui) result (f)
          implicit none
          double precision, dimension(:), intent(in)    :: ui
          double precision                              :: f
          double precision, dimension(1:size(ui))       :: fvec
          integer                                       :: idir
          fvec(:) = func_of_u(ui(:))
          f = 0.0d0
          do idir=1, ndir
             f = f + fvec(idir)**2 / gamma(ix1, idir, idir)
          end do
          !f = dot_product(fvec(:),fvec(:))
          f = 0.5d0 * f 
       end function fmin_of_f
  end subroutine gremhd_implicit_update_grid

end module mod_gremhd_phys_implicit_update
