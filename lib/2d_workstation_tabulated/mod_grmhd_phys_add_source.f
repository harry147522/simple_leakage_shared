module mod_grmhd_phys_add_source
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_add_source_init

contains

  !> Initialize the module
  subroutine grmhd_phys_add_source_init()
    use mod_global_parameters
    use mod_geometry

    integer :: itr, idir

    if (coordinate /= cartesian) then
       ! we need geom source terms
       phys_add_source_geom     => grmhd_add_source_geom
    end if

    phys_add_source          => grmhd_add_source
  end subroutine grmhd_phys_add_source_init

  subroutine get_g_up_munu(g, alp, beta, gamma, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: alp(ixImin1:ixImax1,ixImin2:ixImax2),&
        beta(ixImin1:ixImax1,ixImin2:ixImax2, 1:3)
    double precision, intent(in)    :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:3, 1:3)
    double precision, intent(out)   :: g(ixImin1:ixImax1,ixImin2:ixImax2, 0:3,&
        0:3)

    integer                         :: mu, nu

    g(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 0:3,0:3) = 0.0d0

    g(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 0,0) = - 1.0d0 / alp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2
    do mu = 1, 3
       g(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 0, mu) = beta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mu) / alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2
       g(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mu, 0) = g(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, 0, mu)  
    end do

    ! this is actually \gamma^{ij}
    do mu = 1, 3
       g(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mu,&
           mu) = 1.0d0 / gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mu, mu)
    end do

    ! g^{ij} = \gamma^{ij} - \beta^i \beta^j / \alpha^2
    do mu = 1, 3
       do nu = 1, 3
          g(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mu, nu) = g(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, mu, nu) - beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mu) * beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              nu) / alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2
       end do
    end do

  end subroutine get_g_up_munu

  !> Add gravitational source terms to w
  subroutine grmhd_add_source(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, primCT, cons, x, qsourcesplit, active)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(in)    :: primCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nprim)
    double precision, intent(inout) :: cons(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active !< Needs to be set to true when active

    ! MHD variables needed when using GLM
    double precision                :: grad_Bphi(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndir)

    
    double precision                :: cot_theta(ixImin1:ixImax1,&
       ixImin2:ixImax2)
   
    double precision                :: alppsi6(ixImin1:ixImax1,&
       ixImin2:ixImax2) 
    double precision                :: add_source(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ncons)
    double precision                :: htot(ixImin1:ixImax1,ixImin2:ixImax2) !enthalpy
    double precision                :: Ptot(ixImin1:ixImax1,ixImin2:ixImax2) !enthalpy
    double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2) !Lorentz factor
    ! metric variables
    double precision                :: alp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: psi(ixImin1:ixImax1,ixImin2:ixImax2),&
        psi4(ixImin1:ixImax1,ixImin2:ixImax2), psi6(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                :: betai(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3)
    ! dervitaves of the metric variables
    double precision                :: dalp(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3), dpsi(ixImin1:ixImax1,ixImin2:ixImax2,1:3)
    double precision                :: dbeta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)

    ! 3-metric gamma_ij
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    ! extrinsic curvature K_ij
    double precision                :: K_ij(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    ! 4-metric g^{\mu\nu}
    double precision                :: g(ixImin1:ixImax1,ixImin2:ixImax2,0:3,&
       0:3)
    ! dervitaves of the 3-metric gamma_{jk},i
    double precision                :: Dgamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3,1:3)
    ! Christoffel symbols of the reference 3-metric gamma_hat
    double precision                :: christoffel(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3,1:3,1:3)
    ! covariant dervitaves of beta D_i beta^k
    double precision                :: D_beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    ! energy-momentum tensor T^{\mu\nu}
    double precision                :: T(ixImin1:ixImax1,ixImin2:ixImax2,0:3,&
       0:3)
    ! 4-velocity u^{\mu}
    double precision                :: u(ixImin1:ixImax1,ixImin2:ixImax2,0:3)
    double precision                :: bmu(ixImin1:ixImax1,ixImin2:ixImax2,&
       0:3)
    integer                         :: iw,idir,jdir,kdir

    !-----------------------------------------------------------------------
    ! Part 0: For the split source part
    !-----------------------------------------------------------------------
    ! NO split source terms at the moment
    if (qsourcesplit) return
    
    add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1:ncons) = 0.0d0

    ! fixme: maybe use pointers?
    ! initialize the metric
    psi(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = primCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_)
    psi4(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = psi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**4
    psi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = psi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**6
    alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = primCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,alp_)
    betai(ixImin1:ixImax1,ixImin2:ixImax2,1:3) = 0.0d0
    do idir = 1, ndir
       betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = primCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, beta(idir))
    end do
    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = alp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * psi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**6
    ! volume averaged Christoffel symbols
    Christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3,1:3,&
       1:3) = block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3,1:3,1:3)
    
    !-----------------------------------------------------------------------
    ! Part 1: normal source terms 
    !-----------------------------------------------------------------------
   
    if ( type_divb == divb_glm ) then
       add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bphi_cons_) = alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bphi_) * divB_glm_kappa
       do idir = 1, ndir
          call gradient( primCT(ixImin1:ixImax1,ixImin2:ixImax2, Bphi_),&
             ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
             idir,grad_Bphi(ixImin1:ixImax1,ixImin2:ixImax2,idir) )
          add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(idir)) = - alp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) * grad_Bphi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              idir)
       end do
    end if
    
    !-----------------------------------------------------------------------
    ! Part 2: gravitational source terms, only needed when use_GR = .true.
    !-----------------------------------------------------------------------
    if ( .not. use_GR ) return

    ! fixme: modified primCT to normal variable names
    ! calculate derivitives of the metric variables
    dalp(ixImin1:ixImax1,ixImin2:ixImax2,1:3) = 0.0d0
    dbeta(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3) = 0.0d0
    dpsi(ixImin1:ixImax1,ixImin2:ixImax2,1:3) = 0.0d0
    do idir = 1, ndim
       call partial_d( primCT(ixImin1:ixImax1,ixImin2:ixImax2,alp_) ,ixImin1,&
          ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idir,&
          dalp(ixImin1:ixImax1,ixImin2:ixImax2,idir) )
       call partial_d( primCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_) ,ixImin1,&
          ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idir,&
          dpsi(ixImin1:ixImax1,ixImin2:ixImax2,idir) )
       dpsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = dpsi(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,idir) / psi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       do jdir = 1, ndir
          call partial_d( primCT(ixImin1:ixImax1,ixImin2:ixImax2,beta(jdir)) ,&
             ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
             idir,dbeta(ixImin1:ixImax1,ixImin2:ixImax2,jdir,idir))
       end do
    end do

    call grmhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, primCT(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        bmu=bmu(ixImin1:ixImax1,ixImin2:ixImax2,0:ndir),&
        lfac=lfac(ixImin1:ixImax1,ixImin2:ixImax2), htot=htot(ixImin1:ixImax1,&
       ixImin2:ixImax2), Ptot=Ptot(ixImin1:ixImax1,ixImin2:ixImax2) )

    ! partial dervitaves of metric D_i gamma_jk
    Dgamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3,1:3) = 0.0d0
    do kdir = 1, ndim
       do idir = 1, 3
          Dgamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,idir,&
             kdir) = 4.0d0 * gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idir) * dpsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
       end do
    end do

    ! Note that g(mu,nu) here is g^{\mu\nu}
    call get_g_up_munu(g(ixImin1:ixImax1,ixImin2:ixImax2,0:3,0:3),&
        alp(ixImin1:ixImax1,ixImin2:ixImax2), betai(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3), gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)

    ! Calculate the 4-velocity
    u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0:3) = 0.0d0
    u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0) = lfac(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) / alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    do idir = 1, ndir
       u(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir) = ( primCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, W_vel(idir)) - lfac(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          idir) / alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) )
    end do

    ! energy-momentum tensor T^{\mu\nu}
    do idir = 0,3
       do jdir = 0,3
          T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             jdir) = primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_) * htot(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) * u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              idir) * u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              jdir) + Ptot(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) * g(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             jdir) - bmu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              idir) * bmu(ixOmin1:ixOmax1,ixOmin2:ixOmax2, jdir)
       enddo
    enddo

    ! covariant derivative of beta: partial_i beta^k + Gamma^k_{ij} beta^j
    D_beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1:3, 1:3) = dbeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, 1:3, 1:3)
    if ( coordinate /= cartesian ) then
       do idir = 1, 3
          do kdir = 1, 3 
             do jdir = 1, 3 
                D_beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2, kdir,&
                    idir) = D_beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2, kdir,&
                    idir) + Christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2, kdir,&
                   idir,jdir) * betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,jdir)
             end do
          end do
       end do
    end if

    ! Now calculate the source terms for HD (mom, tau) in the compact form
    do idir = 1, ndir
       do jdir = 1, 3 
          do kdir = 1, 3 
             add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) + 0.5d0 * Dgamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 jdir,kdir,idir) * ( T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
                0)*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                kdir)*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                jdir) + 2.0d0*T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
                jdir)*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                kdir) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2, jdir,kdir) )
          end do
       end do
    end do

    do idir = 1, ndir
       add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir)) - T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
          0)*alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dalp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, idir) 
       ! T^{0j} * gamma_jk * D_i beta^k
       do kdir = 1, 3 
          do jdir = 1, 3 
             add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(idir)) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
                jdir) * gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,jdir,&
                kdir) * D_beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir,idir)
          end do
       end do
    end do

    ! To compute the tau source term ( add_source(:,:,:,tau_) ), we need to work out the extrinsic curvature K_{ij}
    ! fixme: make get_K_ij as a stand alone subroutine
    K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3,1:3) = 0.0d0

    select case (coordinate)
    case (cartesian)
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,1) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1,1)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( 2.0d0*dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1)  - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2)  )
   
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,2,2)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1)  + 2.0d0*dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2)  )
   
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,3,3)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1)  - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2)  )
    case (cylindrical)
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,1) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1,1)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( 2.0d0*dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1) - betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,r_)  - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2)  )
   
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,2,2)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1) - betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,r_)  + 2.0d0 * dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2,2)  )
   
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,3,3)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1) + 2.0d0 * betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_)  - dbeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,2,2)  )
    case (spherical)
       
       cot_theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dcos(x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,theta_))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          theta_))
      

       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,1) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1,1)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( 2.0d0*dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1) - 2.0d0*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_)  - dbeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,2,2) - betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)*cot_theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  )
   
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,2,2)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( -dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1) + betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,r_)  + 2.0d0*dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
          2) - betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)*cot_theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  )
   
       K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) = gamma(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,3,3)/(3.0d0*alp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)) * ( -dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          1) + betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,r_)  - dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
          2) + 2.0d0*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)*cot_theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  )
    end select
    
    K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,2) = ( dbeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1,2)*gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
       1) + dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,1)*gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2,2) ) / (2.0d0*alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,3) = ( dbeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,3,1)*gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
       3) + dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,3)*gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1,1) ) / (2.0d0*alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,3) = ( dbeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,3,2)*gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
       3) + dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,3)*gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2,2) ) / (2.0d0*alp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    ! K_ji=K_ij
    do idir=1,2
       do jdir=idir+1,3
          K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,jdir,&
             idir) = K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,jdir)
       end do
    end do
   


    ! Now work out the source term of tau
    do idir=1,3
       add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          tau_) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          tau_) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
          0)*( -betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          idir)*dalp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          idir) ) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
          idir)*( -dalp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) )

       do jdir=1,3
          add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             tau_) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             tau_) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
             idir)*( 2.0d0*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             jdir)*K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,jdir) )
       enddo

    enddo

    do idir=1,3
       do jdir=1,3
          add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             tau_) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             tau_) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,0,&
             0) * ( K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             jdir)*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idir)*betai(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             jdir) ) + T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             jdir) * K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,jdir)
       enddo
    enddo

    cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_) = cons(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, tau_) + qdt * alppsi6(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tau_)
    do idir =1, ndir
      cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = cons(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(idir)) + qdt * alppsi6(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idir))
    end do

    ! Now work out the source term of Bphi if glm is being used
    if ( type_divb == divb_glm ) then
       do idir=1, ndir ! becareful! 
          add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             Bphi_cons_) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             Bphi_cons_) - alp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) * K_ij(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idir) / gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idir) + primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bvec(idir)) * dalp(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir)
       enddo
       do idir = 1, ndir
          do jdir = 1, ndir
             add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 Bcons(idir)) = add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 Bcons(idir)) - primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 Bvec(jdir)) * dbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,&
                 jdir) 
          end do
       end do
       ! add the source terms on cons
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bphi_cons_) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bphi_cons_) +  qdt * psi6(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bphi_cons_)
       do idir = 1, ndir
          cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(idir)) =  cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(idir)) + qdt * psi6(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) * add_source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(idir))
       end do
    end if

  end subroutine grmhd_add_source

  !> Add geometrical source terms to w
  subroutine grmhd_add_source_geom(qdt, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, consCT, primCT, cons, x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(in)    :: consCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons), primCT(ixImin1:ixImax1,ixImin2:ixImax2, 1:nprim)
    double precision, intent(inout) :: cons(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ncons)

    integer                         :: iw,idir,jdir,kdir, h1xmin1,h1xmin2,&
       h1xmax1,h1xmax2, h2xmin1,h2xmin2,h2xmax1,h2xmax2
    double precision                :: source(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: fluxCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3)
    double precision                :: alppsi6(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                :: B_dot_v(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                :: v_hat(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)
    double precision                :: gamma(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3,1:3) !metric gamma_ij
    double precision                :: B_over_W(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:ndir) 
    double precision                :: Ptot(ixImin1:ixImax1,ixImin2:ixImax2) !total pressure
    double precision                :: lfac(ixImin1:ixImax1,ixImin2:ixImax2) !Lorentz factor

    call grmhd_get_intermediate_variables(ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, primCT(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nprim), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim),&
        gamma=gamma(ixImin1:ixImax1,ixImin2:ixImax2,1:3,1:3),&
        lfac=lfac(ixImin1:ixImax1,ixImin2:ixImax2),&
        v_hat=v_hat(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
        B_dot_v=B_dot_v(ixImin1:ixImax1,ixImin2:ixImax2),&
        Ptot=Ptot(ixImin1:ixImax1,ixImin2:ixImax2) )

    ! alp * psi**6
    alppsi6(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = primCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,alp_) * primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)**6

    ! B^i / W
    do idir = 1, ndir
       B_over_W(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idir) = primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Bvec(idir)) / lfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! Only Momentum flux f^i_j is needed only
    fluxCT = 0.0d0
    do jdir = 1, ndir
       do idir = 1, ndir
         fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,&
             jdir) = v_hat(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            idir) * consCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(jdir)) - alppsi6(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) * ( B_over_W(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             jdir) + B_dot_v(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) * primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             W_vel(jdir)) ) * gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2, jdir,&
             jdir) * B_over_W(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir)
       end do
    end do
    ! and the pressure terms
    do jdir = 1, 3
       fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, jdir,&
          jdir) = fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, jdir,&
           jdir) + alppsi6(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * Ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    select case (coordinate)
    case (cylindrical)
       ! geo_source[r] = Gamma^phi_{r phi}*f^phi_phi
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,1,&
          3) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3)
       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(r_)) = cons(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(r_)) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)

    case (spherical)
       ! geo_source[r] = Gamma^2_{12}*f^2_2 + Gamma^3_{13}*f^3_3
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,1,&
          2) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
          2) +block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,1,&
          3) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) 

       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(r_)) = cons(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(r_)) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)

       
 !geo_source[theta] = Gamma1_22,Gamma1_22,22*f2_1,22*f2_1 + Gamma2_21,Gamma2_21,21*f1_2,21*f1_2 + Gamma3_23,Gamma3_23,23*f3_3,23*f3_3
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,2,&
          2) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
          1) +block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2,&
          1) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
          2) +block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,2,&
          3) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) 

       cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(2)) = cons(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(2)) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       ! geo_source[phi] = 0
      
    end select

    ! now work on the geo_source for B field,
    ! note that this is only for the case when using GLM as it introduce an
    ! extra term
    if (type_divb == divb_glm) then
       ! B field flux f^ij, with the extra term only
       fluxCT = 0.0d0
       do jdir = 1, ndir
          do idir = 1, ndir
            fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,&
                jdir) = primCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                beta(jdir)) * consCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                Bcons(idir))
          end do
       end do
       select case (coordinate)
       case (cylindrical)
          call mpistop("no cylindrical yet")
       case (spherical)
          ! geo_source_B[r] = Gamma^1_{22}*f^22 + Gamma^1_{33}*f^33
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
             block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,2,&
             2) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
             2) +block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,3,&
             3) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) 
          cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(r_)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(r_)) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   
          
 !geo_source_B[theta] = Gamma2_12,Gamma2_12,12*f12,12*f12 + Gamma2_21,Gamma2_21,21*f21,21*f21 + Gamma2_33,Gamma2_33,33*f33,33*f33
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
             block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,1,&
             2) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
             2) +block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,2,&
             1) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
             1) +block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,3,&
             3) * fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,3) 
          cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(2)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(2)) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   
 !geo_source_B[phi] = Gamma3_13,Gamma3_13,13*f13,13*f13 + Gamma3_31,Gamma3_31,31*f31,31*f31 + Gamma3_12,Gamma3_12,12*f21,12*f21 + Gamma3_21,Gamma3_21,21*f12,21*f12
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
             block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,1,&
             3) * ( fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
             3) + fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
             1) )+block%christoffel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,2,&
             1) * ( fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
             1) + fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,2) ) 

          cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(3)) = cons(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              Bcons(3)) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
         
       end select
    end if
  end subroutine grmhd_add_source_geom

end module mod_grmhd_phys_add_source
