module mod_grhd_ccsn_leakage_phys_add_source
  use mod_physics
  use mod_grhd_ccsn_leakage_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grhd_ccsn_leakage_phys_add_source_init

contains

  !> Initialize the module
  subroutine grhd_ccsn_leakage_phys_add_source_init()
    use mod_global_parameters
    use mod_geometry

    integer :: itr, idir

    if (coordinate /= cartesian) then
       ! we need geom source terms
       phys_add_source_geom     => grhd_ccsn_leakage_add_source_geom
    end if
    phys_add_source          => grhd_ccsn_leakage_add_source
  end subroutine grhd_ccsn_leakage_phys_add_source_init

  subroutine get_g_up_munu(g, alp, beta, gamma, ixImin1,ixImax1, ixOmin1,&
     ixOmax1)
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: alp(ixImin1:ixImax1),&
        beta(ixImin1:ixImax1, 1:3)
    double precision, intent(in)    :: gamma(ixImin1:ixImax1, 1:3, 1:3)
    double precision, intent(out)   :: g(ixImin1:ixImax1, 0:3, 0:3)
    integer                         :: mu, nu
    g(ixOmin1:ixOmax1, 0:3,0:3) = 0.0d0
    g(ixOmin1:ixOmax1, 0,0) = - 1.0d0 / alp(ixOmin1:ixOmax1)**2
    do mu = 1, 3
       g(ixOmin1:ixOmax1, 0, mu) = beta(ixOmin1:ixOmax1,&
           mu) / alp(ixOmin1:ixOmax1)**2
       g(ixOmin1:ixOmax1, mu, 0) = g(ixOmin1:ixOmax1, 0, mu)  
    end do
    ! this is actually \gamma^{ij}
    do mu = 1, 3
       g(ixOmin1:ixOmax1, mu, mu) = 1.0d0 / gamma(ixOmin1:ixOmax1, mu, mu)
    end do
    ! g^{ij} = \gamma^{ij} - \beta^i \beta^j / \alpha^2
    do mu = 1, 3
       do nu = 1, 3
          g(ixOmin1:ixOmax1, mu, nu) = g(ixOmin1:ixOmax1, mu,&
              nu) - beta(ixOmin1:ixOmax1, mu) * beta(ixOmin1:ixOmax1,&
              nu) / alp(ixOmin1:ixOmax1)**2
       end do
    end do

  end subroutine get_g_up_munu

  !> Add gravitational source terms to w
  subroutine grhd_ccsn_leakage_add_source(qdt, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, primCT, cons, x, qsourcesplit, active)
    use mod_global_parameters
    use mod_geometry
    use mod_grhd_ccsn_leakage_simple_leakage, only : bounce
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(in)    :: primCT(ixImin1:ixImax1, 1:nprim)
    double precision, intent(inout) :: cons(ixImin1:ixImax1, 1:ncons)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active !< Needs to be set to true when active

    
    double precision                :: add_source(ixImin1:ixImax1,1:ncons)
    double precision                :: add_source_leakage(ixImin1:ixImax1,&
       ndir+2)
    double precision                :: h(ixImin1:ixImax1) ! enthalpy
    double precision                :: lfac(ixImin1:ixImax1) ! Lorentz factor
    ! metric variables
    double precision                :: alp(ixImin1:ixImax1)
    double precision                :: psi(ixImin1:ixImax1),&
        psi4(ixImin1:ixImax1)
    double precision                :: betai(ixImin1:ixImax1,1:3)
    ! dervitaves of the metric variables
    double precision                :: dalp(ixImin1:ixImax1,1:3),&
        dpsi(ixImin1:ixImax1,1:3)
    double precision                :: dbeta(ixImin1:ixImax1,1:3,1:3)

    ! 3-metric gamma_ij
    double precision                :: gamma(ixImin1:ixImax1,1:3,1:3)
    ! extrinsic curvature K_ij
    double precision                :: K_ij(ixImin1:ixImax1,1:3,1:3)
    ! 4-metric g^{\mu\nu}
    double precision                :: g(ixImin1:ixImax1,0:3,0:3)
    ! dervitaves of the 3-metric gamma_{jk},i
    double precision                :: Dgamma(ixImin1:ixImax1,1:3,1:3,1:3)
    ! Christoffel symbols of the reference 3-metric gamma_hat
    double precision                :: christoffel(ixImin1:ixImax1,1:3,1:3,&
       1:3)
    ! covariant dervitaves of beta D_i beta^k
    double precision                :: D_beta(ixImin1:ixImax1,1:3,1:3)
    ! energy-momentum tensor T^{\mu\nu}
    double precision                :: T(ixImin1:ixImax1,0:3,0:3)
    ! 4-velocity u^{\mu}
    double precision                :: u(ixImin1:ixImax1,0:3)
    integer                         :: iw,idir,jdir,kdir

    !-----------------------------------------------------------------------
    ! Part 0: For the split source part
    !-----------------------------------------------------------------------
    ! NO split source terms at the moment
    if (qsourcesplit) return

    add_source(ixOmin1:ixOmax1, 1:ncons) = 0.0d0
    
    ! fixme: maybe use pointers?
    ! initialize the metric
    psi(ixOmin1:ixOmax1) = primCT(ixOmin1:ixOmax1,psi_)
    psi4(ixOmin1:ixOmax1) = psi(ixOmin1:ixOmax1)**4
    alp(ixOmin1:ixOmax1) = primCT(ixOmin1:ixOmax1,alp_)
    betai(ixImin1:ixImax1,1:3) = 0.0d0
    do idir = 1, ndir
       betai(ixOmin1:ixOmax1,idir) = primCT(ixOmin1:ixOmax1, beta(idir))
    end do
    ! volume averaged Christoffel symbols
    Christoffel(ixOmin1:ixOmax1,1:3,1:3,1:3) = &
       block%christoffel(ixOmin1:ixOmax1,1:3,1:3,1:3)
    
    !-----------------------------------------------------------------------
    ! Part 1: normal source terms 
    !-----------------------------------------------------------------------
    ! nothing at the moment
   
    !-----------------------------------------------------------------------
    ! Part 2: gravitational source terms, only needed when use_GR = .true.
    !-----------------------------------------------------------------------
    if ( .not. use_GR ) return

    call grhd_ccsn_leakage_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,&
       ixOmax1, primCT(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,1:3,1:3), lfac=lfac(ixImin1:ixImax1),&
        h=h(ixImin1:ixImax1))

    ! fixme: modified primCT to normal variable names
    ! calculate derivitives of the metric variables
    dalp(ixImin1:ixImax1,1:3) = 0.0d0
    dbeta(ixImin1:ixImax1,1:3,1:3) = 0.0d0
    dpsi(ixImin1:ixImax1,1:3) = 0.0d0
    do idir = 1, ndim
       call partial_d( primCT(ixImin1:ixImax1,alp_) ,ixImin1,ixImax1,ixOmin1,&
          ixOmax1,idir,dalp(ixImin1:ixImax1,idir) )
       call partial_d( primCT(ixImin1:ixImax1,psi_) ,ixImin1,ixImax1,ixOmin1,&
          ixOmax1,idir,dpsi(ixImin1:ixImax1,idir) )
       dpsi(ixOmin1:ixOmax1,idir) = dpsi(ixOmin1:ixOmax1,&
          idir) / psi(ixOmin1:ixOmax1)
       do jdir = 1, ndir
          call partial_d( primCT(ixImin1:ixImax1,beta(jdir)) ,ixImin1,ixImax1,&
             ixOmin1,ixOmax1,idir,dbeta(ixImin1:ixImax1,jdir,idir))
       end do
    end do
 
    

    ! covariant derivative of beta: partial_i beta^k + Gamma^k_{ij} beta^j
    D_beta(ixOmin1:ixOmax1, 1:3, 1:3) = dbeta(ixOmin1:ixOmax1, 1:3, 1:3)
    if ( coordinate /= cartesian ) then
       do idir = 1, 3
          do kdir = 1, 3 
             do jdir = 1, 3 
                D_beta(ixOmin1:ixOmax1, kdir, idir) = D_beta(ixOmin1:ixOmax1,&
                    kdir, idir) + Christoffel(ixOmin1:ixOmax1, kdir,idir,&
                   jdir) * betai(ixOmin1:ixOmax1,jdir)
             end do
          end do
       end do
    end if

    ! dervitaves of metric D_i gamma_jk
    Dgamma(ixImin1:ixImax1,1:3,1:3,1:3) = 0.0d0
    do kdir = 1, ndim
       do idir = 1, 3
          Dgamma(ixOmin1:ixOmax1,idir,idir,&
             kdir) = 4.0d0 * gamma(ixOmin1:ixOmax1,idir,&
             idir) * dpsi(ixOmin1:ixOmax1,kdir)
       end do
    end do

    ! Note that g(mu,nu) here is g^{\mu\nu}
    call get_g_up_munu(g(ixImin1:ixImax1,0:3,0:3), alp(ixImin1:ixImax1),&
        betai(ixImin1:ixImax1,1:3), gamma(ixImin1:ixImax1,1:3,1:3), ixImin1,&
       ixImax1, ixOmin1,ixOmax1)

    ! Calculate the 4-velocity
    u(ixOmin1:ixOmax1,0:3) = 0.0d0
    u(ixOmin1:ixOmax1,0) = lfac(ixOmin1:ixOmax1) / alp(ixOmin1:ixOmax1)
    do idir = 1, ndir
       u(ixOmin1:ixOmax1, idir) = primCT(ixOmin1:ixOmax1,&
           W_vel(idir)) - lfac(ixOmin1:ixOmax1) * betai(ixOmin1:ixOmax1,&
          idir) / alp(ixOmin1:ixOmax1) 
    end do

    ! energy-momentum tensor T^{\mu\nu}
    do idir = 0,3
       do jdir = 0,3
          T(ixOmin1:ixOmax1,idir,jdir) = primCT(ixOmin1:ixOmax1,&
              rho_) * h(ixOmin1:ixOmax1) * u(ixOmin1:ixOmax1,&
              idir) * u(ixOmin1:ixOmax1, jdir) + primCT(ixOmin1:ixOmax1,&
              press_) * g(ixOmin1:ixOmax1,idir,jdir)
       enddo
    enddo

    ! Now calculate the source terms for HD (mom, tau) in the compact form
    if ( evolve_hydro ) then
       do idir = 1, ndir
          do jdir = 1, 3 
             do kdir = 1, 3 
                add_source(ixOmin1:ixOmax1,&
                    mom(idir)) = add_source(ixOmin1:ixOmax1,&
                    mom(idir)) + 0.5d0 * Dgamma(ixOmin1:ixOmax1, jdir, kdir,&
                    idir) * ( T(ixOmin1:ixOmax1, 0, 0) * betai(ixOmin1:ixOmax1,&
                    kdir) * betai(ixOmin1:ixOmax1,&
                    jdir) + 2.0d0 * T(ixOmin1:ixOmax1, 0,&
                    jdir) * betai(ixOmin1:ixOmax1, kdir) + T(ixOmin1:ixOmax1,&
                    jdir, kdir) )
             end do
          end do
       end do
   
       do idir = 1, ndir
          add_source(ixOmin1:ixOmax1, mom(idir)) = add_source(ixOmin1:ixOmax1,&
              mom(idir)) - T(ixOmin1:ixOmax1, 0,&
              0) * alp(ixOmin1:ixOmax1) * dalp(ixOmin1:ixOmax1, idir) 
          ! T^{0j} * gamma_jk * D_i beta^k
          do kdir = 1, 3 
             do jdir = 1, 3 
                add_source(ixOmin1:ixOmax1,&
                    mom(idir)) = add_source(ixOmin1:ixOmax1,&
                    mom(idir)) + T(ixOmin1:ixOmax1, 0,&
                    jdir) * gamma(ixOmin1:ixOmax1, jdir,&
                    kdir) * D_beta(ixOmin1:ixOmax1, kdir, idir)
             end do
          end do
       end do
   
       ! To compute the tau source term, we need to work out the extrinsic curvature K_{ij}
       K_ij(ixOmin1:ixOmax1,1:3,1:3) = 0.0d0
   
       select case (coordinate)
       case (cartesian)
          K_ij(ixOmin1:ixOmax1,1,1) = gamma(ixOmin1:ixOmax1,1,&
             1)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( 2.0d0*dbeta(ixOmin1:ixOmax1,1,&
             1)   )
      
          K_ij(ixOmin1:ixOmax1,2,2) = gamma(ixOmin1:ixOmax1,2,&
             2)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( - dbeta(ixOmin1:ixOmax1,1,&
             1)   )
      
          K_ij(ixOmin1:ixOmax1,3,3) = gamma(ixOmin1:ixOmax1,3,&
             3)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( - dbeta(ixOmin1:ixOmax1,1,&
             1)   )
       case (cylindrical)
          K_ij(ixOmin1:ixOmax1,1,1) = gamma(ixOmin1:ixOmax1,1,&
             1)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( 2.0d0*dbeta(ixOmin1:ixOmax1,1,&
             1) - betai(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,r_)   )
      
          K_ij(ixOmin1:ixOmax1,2,2) = gamma(ixOmin1:ixOmax1,2,&
             2)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( - dbeta(ixOmin1:ixOmax1,1,&
             1) - betai(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,r_)   )
      
          K_ij(ixOmin1:ixOmax1,3,3) = gamma(ixOmin1:ixOmax1,3,&
             3)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( - dbeta(ixOmin1:ixOmax1,1,&
             1) + 2.0d0 * betai(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,r_)   )
       case (spherical)
          
   
          K_ij(ixOmin1:ixOmax1,1,1) = gamma(ixOmin1:ixOmax1,1,&
             1)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( 2.0d0*dbeta(ixOmin1:ixOmax1,1,&
             1) - 2.0d0*betai(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,r_)   )
      
          K_ij(ixOmin1:ixOmax1,2,2) = gamma(ixOmin1:ixOmax1,2,&
             2)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( -dbeta(ixOmin1:ixOmax1,1,&
             1) + betai(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,r_)   )
      
          K_ij(ixOmin1:ixOmax1,3,3) = gamma(ixOmin1:ixOmax1,3,&
             3)/(3.0d0*alp(ixOmin1:ixOmax1)) * ( -dbeta(ixOmin1:ixOmax1,1,&
             1) + betai(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,r_)   )
       end select
       
   
       ! Now work out the source term of tau
       do idir=1,3
          add_source(ixOmin1:ixOmax1,tau_) = add_source(ixOmin1:ixOmax1,&
             tau_) + T(ixOmin1:ixOmax1,0,0) * ( -betai(ixOmin1:ixOmax1,&
             idir)*dalp(ixOmin1:ixOmax1,idir) ) + T(ixOmin1:ixOmax1,0,&
             idir) * ( -dalp(ixOmin1:ixOmax1,idir) ) 
       enddo

   
       do idir=1,3
          do jdir=1,3
             add_source(ixOmin1:ixOmax1,tau_) = add_source(ixOmin1:ixOmax1,&
                tau_) + T(ixOmin1:ixOmax1,0,&
                idir)*( 2.0d0*betai(ixOmin1:ixOmax1,jdir)*K_ij(ixOmin1:ixOmax1,&
                idir,jdir) ) + T(ixOmin1:ixOmax1,0,0) * ( K_ij(ixOmin1:ixOmax1,&
                idir,jdir)*betai(ixOmin1:ixOmax1,idir)*betai(ixOmin1:ixOmax1,&
                jdir) ) + T(ixOmin1:ixOmax1,idir,jdir) * K_ij(ixOmin1:ixOmax1,&
                idir,jdir)
          enddo
       enddo
   
       cons(ixOmin1:ixOmax1, tau_) = cons(ixOmin1:ixOmax1,&
           tau_) + qdt * alp(ixOmin1:ixOmax1) * psi(ixOmin1:ixOmax1)**6 * &
          add_source(ixOmin1:ixOmax1, tau_)
       do idir =1, ndir
         cons(ixOmin1:ixOmax1, mom(idir)) = cons(ixOmin1:ixOmax1,&
             mom(idir)) + qdt * alp(ixOmin1:ixOmax1) * psi(ixOmin1:ixOmax1)**6 &
            * add_source(ixOmin1:ixOmax1,mom(idir))
       end do
    end if ! end hydro source


 
!write(*,*) "passed phys_add_source"
         
!           write(*,*)  cons(5, coolingsource(1)), x(5,r_), "coolingsource"
!stop  "you arrived add_source"
    if (bounce) then
        
       
       
 ! write(*,*) "add_source"
 ! write(*,*) cons(5, coolingsource(:))
       do idir = 1, ndir+2
         add_source_leakage(ixOmin1:ixOmax1, idir) = cons(ixOmin1:ixOmax1,&
             coolingsource(idir))
       enddo
 
!  !     write(*,*) "wrong path"
       cons(ixOmin1:ixOmax1, tau_) = cons(ixOmin1:ixOmax1,&
           tau_) + qdt * add_source_leakage(ixOmin1:ixOmax1, ndir+1)
       cons(ixOmin1:ixOmax1, ye_con_) = cons(ixOmin1:ixOmax1,&
           ye_con_) + qdt * add_source_leakage(ixOmin1:ixOmax1, ndir+2)
       do idir =1, ndir
         cons(ixOmin1:ixOmax1, mom(idir)) = cons(ixOmin1:ixOmax1,&
             mom(idir)) + qdt * add_source_leakage(ixOmin1:ixOmax1, idir)
       end do
   endif


  
!write(*,*) add_source(5,tau_) , x(5,r_), qdt,"normal ource"
  end subroutine grhd_ccsn_leakage_add_source

  !> Add geometrical source terms to w
  subroutine grhd_ccsn_leakage_add_source_geom(qdt, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, consCT, primCT, cons, x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(in)    :: consCT(ixImin1:ixImax1, 1:ncons),&
        primCT(ixImin1:ixImax1, 1:nprim)
    double precision, intent(inout) :: cons(ixImin1:ixImax1, 1:ncons)

    double precision                :: source(ixImin1:ixImax1)
    double precision                :: fluxCT(ixImin1:ixImax1,1:3,1:3)
    double precision                :: alppsi6(ixImin1:ixImax1),&
        lfac(ixImin1:ixImax1)
    double precision                :: v_hat(ixImin1:ixImax1,1:ndir)
    integer                         :: iw,idir,jdir,kdir, h1xmin1,h1xmax1
    integer                         :: iphi = 3

    if ( .not. evolve_hydro ) return

    call grhd_ccsn_leakage_get_intermediate_variables(ixImin1,ixImax1, ixOmin1,&
       ixOmax1, primCT(ixImin1:ixImax1, 1:nprim), x(ixImin1:ixImax1, 1:ndim),&
        lfac=lfac(ixImin1:ixImax1) )

    ! not that v_hat here is v_hat = alp * v - beta, the alp is factored out
    do idir=1,ndir
       v_hat(ixOmin1:ixOmax1, idir) = primCT(ixOmin1:ixOmax1,&
           alp_) * primCT(ixOmin1:ixOmax1,&
           W_vel(idir)) / lfac(ixOmin1:ixOmax1) - primCT(ixOmin1:ixOmax1,&
           beta(idir))
    end do
    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1) = primCT(ixOmin1:ixOmax1,&
       alp_) * primCT(ixOmin1:ixOmax1,psi_)**6

    ! Only Momentum flux f^i_j is needed only
    fluxCT = 0.0d0
    do jdir = 1, ndir
       do idir = 1, ndir
         fluxCT(ixOmin1:ixOmax1, idir, jdir) = v_hat(ixOmin1:ixOmax1,&
            idir) * consCT(ixOmin1:ixOmax1, mom(jdir))
       end do
    end do
    ! and the pressure terms
    do jdir = 1, 3
       fluxCT(ixOmin1:ixOmax1, jdir,jdir) = ( fluxCT(ixOmin1:ixOmax1, jdir,&
           jdir) + alppsi6(ixOmin1:ixOmax1) * primCT(ixOmin1:ixOmax1,&
           press_) )
    end do

    select case (coordinate)
    case (cylindrical)
       ! geo_source[r] = Gamma^phi_{r phi}*f^phi_phi
       if ( phi_ == 2 ) iphi = 2
       source(ixOmin1:ixOmax1) = block%christoffel(ixOmin1:ixOmax1,iphi,r_,&
          iphi) * fluxCT(ixOmin1:ixOmax1,iphi,iphi)
       cons(ixOmin1:ixOmax1, mom(r_)) = cons(ixOmin1:ixOmax1,&
           mom(r_)) + qdt * source(ixOmin1:ixOmax1)

    case (spherical)
       ! geo_source[r] = Gamma^2_{12}*f^2_2 + Gamma^3_{13}*f^3_3
       source(ixOmin1:ixOmax1) = block%christoffel(ixOmin1:ixOmax1,2,1,&
          2) * fluxCT(ixOmin1:ixOmax1,2,2) +block%christoffel(ixOmin1:ixOmax1,&
          3,1,3) * fluxCT(ixOmin1:ixOmax1,3,3) 

       cons(ixOmin1:ixOmax1, mom(r_)) = cons(ixOmin1:ixOmax1,&
           mom(r_)) + qdt * source(ixOmin1:ixOmax1)

       
    end select
  end subroutine grhd_ccsn_leakage_add_source_geom

end module mod_grhd_ccsn_leakage_phys_add_source
