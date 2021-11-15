#include "cpp_macros.h"
!> Module which contains multigrid procedures for the L operator for (psi - 1)
module m_cfc_psi
  use m_data_structures
  use m_finite_difference

  implicit none
  private

  public :: cfc_psi_set_methods

contains

  subroutine cfc_psi_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    mg%vector_equation = .False.

    select case (mg%geometry_type)
#if NDIM == 3
    case (mg_cartesian)
       mg%box_op => box_lpsi

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_lpsi
       case default
          error stop "cfc_psi_set_methods: unsupported smoother type"
       end select
#endif

#if NDIM != 1
    case (mg_cylindrical)
       mg%box_op => box_clpsi

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_clpsi
       case default
          error stop "cfc_psi_set_methods: unsupported smoother type"
       end select
#endif

    case (mg_spherical)
       mg%box_op => box_slpsi

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_slpsi
       case default
          error stop "cfc_psi_set_methods: unsupported smoother type"
       end select

    case default
       error stop "cfc_psi_set_methods: unsupported geometry"
    end select

  end subroutine cfc_psi_set_methods

  !> Perform L operator on a box cartesian geometry
  subroutine box_lpsi(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: u(-1:1, NDIM)
    real(dp)                  :: dr(NDIM), idr2(NDIM)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      do i = 1, nc
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         call Laplacian(cc(i, i_out), u, idr2)

         ! nonlinear source terms
         cc(i, i_out) = cc(i, i_out) &
                      +( cc(i, f1) / ( 1.0d0 + cc(i, n) ) &
                       + cc(i, f2) / ( 1.0d0 + cc(i, n) )**7 )
      end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc
            ! Laplacian
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            call Laplacian(cc(i, j, i_out), u, idr2)

            ! nonlinear source terms
            cc(i, j, i_out) = cc(i, j, i_out) + &
                         ( cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) ) &
                         + cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**7)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               ! Laplacian
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               call Laplacian(cc(i, j, k, i_out), u, idr2)
   
               ! nonlinear source terms
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + &
                            ( cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) ) &
                            + cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**7)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_lpsi

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cartesian geometry.
  subroutine box_gs_lpsi(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di

    real(dp)                  :: u(-1:1, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM)
    real(dp)                  :: Lop, dLdu

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2

    i0  = 1
    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         call Laplacian(Lop, u, idr2, dLdu = dLdu)

         ! nonlinear source terms
         Lop = Lop +   ( cc(i, f1) / ( 1.0d0 + cc(i, n) ) &
                         + cc(i, f2) / ( 1.0d0 + cc(i, n) )**7 )

         dLdu = dLdu + (- cc(i, f1) / ( 1.0d0 + cc(i, n) )**2 &
                     - 7.0d0 * cc(i, f2) / ( 1.0d0 + cc(i, n) )**8 )
 
         cc(i, n) = cc(i, n) - ( Lop - cc(i, mg_irhs) ) / dLdu
      end do

#elif NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di

            ! Laplacian 
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            call Laplacian(Lop, u, idr2, dLdu = dLdu)

            ! nonlinear source terms
            Lop = Lop + ( cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) ) &
                         + cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**7)
            dLdu = dLdu &
                       + (- cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) )**2 &
                         - 7.0d0 * cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**8 )
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, j), 1)
   
            do i = i0, nc, di
   
               ! Laplacian 
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               call Laplacian(Lop, u, idr2, dLdu = dLdu)
   
               ! nonlinear source terms
               Lop = Lop + ( cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) ) &
                            + cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**7)
               dLdu = dLdu + (- cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) )**2 &
                            - 7.0d0 * cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**8 )
       
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_lpsi

  !> Perform L operator on a box in cylindrical geometry, using (r,z,phi)
  subroutine box_clpsi(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    real(dp)                  :: dr(NDIM), idr2(NDIM)
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            ! Laplacian
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)
            pre_fac(1) = r_cc(i)
            pre_fac(2) = r_cc(i)**2
            call Laplacian(cc(i, j, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)

            ! nonlinear source terms
            cc(i, j, i_out) = cc(i, j, i_out) &
                       + r_cc(i)**2 * &
                         ( cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) ) &
                         + cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**7)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               ! Laplacian
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(cc(i, j, k, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                          + r_cc(i)**2 * &
                            ( cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) ) &
                            + cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**7)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_clpsi

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cylindrical geometry.
  subroutine box_gs_clpsi(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM)
    real(dp)                  :: Lop, dLdu
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    i0  = 1
    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di

            ! Laplacian 
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)
            pre_fac(1) = r_cc(i)
            pre_fac(2) = r_cc(i)**2
            call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)

            ! nonlinear source terms
            Lop = Lop + r_cc(i)**2 * &
                         ( cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) ) &
                         + cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**7)
            dLdu = dLdu &
                       + r_cc(i)**2 * &
                        (- cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) )**2 &
                         - 7.0d0 * cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**8 )
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, j), 1)
   
            do i = i0, nc, di
   
               ! Laplacian 
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               Lop = Lop + r_cc(i)**2 * &
                           ( cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) ) &
                            + cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**7)
               dLdu = dLdu + r_cc(i)**2 * &
                           (- cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) )**2 &
                            - 7.0d0 * cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**8 )
       
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_clpsi

  !> Perform L operator on a box in spherical geometry, using (r,theta,phi)
  subroutine box_slpsi(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    real(dp)                  :: dr(NDIM), idr2(NDIM)
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)
#if NDIM != 1
    real(dp)                  :: sin_theta_face(1:nc+1)
    real(dp)                  :: sin_theta_cc(nc)
#endif

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])
#if NDIM != 1
    sin_theta_face  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-1.0_dp, i=1,nc+1)])
    sin_theta_cc  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=1,nc)])
#endif

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      do i = 1, nc
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(cc(i, i_out), u, idr2, face_coeff_in = face_coeff)

         ! nonlinear source terms
         cc(i, i_out) = cc(i, i_out) &
                       + r_cc(i)**2 *  &
                       ( cc(i, f1) / ( 1.0d0 + cc(i, n) ) &
                       + cc(i, f2) / ( 1.0d0 + cc(i, n) )**7 )
      end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc
            ! Laplacian
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            face_coeff(0:1, 2) = sin_theta_face(j:j+1)
            pre_fac(1) = sin_theta_cc(j)
            call Laplacian(cc(i, j, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)

            ! nonlinear source terms
            cc(i, j, i_out) = cc(i, j, i_out) &
                       + r_cc(i)**2 * sin_theta_cc(j) *  &
                         ( cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) ) &
                         + cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**7)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               ! Laplacian
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(cc(i, j, k, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                       + r_cc(i)**2 * sin_theta_cc(j)**2 *  &
                            ( cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) ) &
                            + cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**7)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_slpsi

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> spherical geometry.
  subroutine box_gs_slpsi(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM)
    real(dp)                  :: Lop, dLdu
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)
#if NDIM != 1
    real(dp)                  :: sin_theta_face(1:nc+1)
    real(dp)                  :: sin_theta_cc(nc)
#endif

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])
#if NDIM != 1
    sin_theta_face  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-1.0_dp, i=1,nc+1)])
    sin_theta_cc  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=1,nc)])
#endif

    i0  = 1
    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff)

         ! nonlinear source terms
         Lop = Lop + r_cc(i)**2 *  &
                         ( cc(i, f1) / ( 1.0d0 + cc(i, n) ) &
                         + cc(i, f2) / ( 1.0d0 + cc(i, n) )**7 )

         dLdu = dLdu + r_cc(i)**2 *   &
                    (- cc(i, f1) / ( 1.0d0 + cc(i, n) )**2 &
                     - 7.0d0 * cc(i, f2) / ( 1.0d0 + cc(i, n) )**8 )
 
         cc(i, n) = cc(i, n) - ( Lop - cc(i, mg_irhs) ) / dLdu
      end do

#elif NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di

            ! Laplacian 
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            face_coeff(0:1, 2) = sin_theta_face(j:j+1)
            pre_fac(1) = sin_theta_cc(j)
            call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)

            ! nonlinear source terms
            Lop = Lop + r_cc(i)**2 * sin_theta_cc(j) * &
                         ( cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) ) &
                         + cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**7)
            dLdu = dLdu &
                       + r_cc(i)**2 * sin_theta_cc(j) * &
                        (- cc(i, j, f1) / ( 1.0d0 + cc(i, j, n) )**2 &
                         - 7.0d0 * cc(i, j, f2) / ( 1.0d0 + cc(i, j, n) )**8 )
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, j), 1)
   
            do i = i0, nc, di
   
               ! Laplacian 
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               Lop = Lop + r_cc(i)**2 * sin_theta_cc(j)**2 * &
                           ( cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) ) &
                            + cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**7)
               dLdu = dLdu + r_cc(i)**2 * sin_theta_cc(j)**2 * &
                           (- cc(i, j, k, f1) / ( 1.0d0 + cc(i, j, k, n) )**2 &
                            - 7.0d0 * cc(i, j, k, f2) / ( 1.0d0 + cc(i, j, k, n) )**8 )
       
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_slpsi

end module m_cfc_psi
