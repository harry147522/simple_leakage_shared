#include "cpp_macros.h"
!> Module which contains multigrid procedures for a Laplacian operator
module m_laplacian
  use m_data_structures
  use m_finite_difference

  implicit none
  private

  public :: laplacian_set_methods

contains

  subroutine laplacian_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (all(mg%periodic)) then
       ! For a fully periodic Laplacian, remove the mean from the rhs and phi so
       ! that a unique and periodic solution can be found
       mg%subtract_mean = .true.
    end if

    mg%vector_equation = .False.

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_lpl

       select case (mg%smoother_type)
       ! case (smoother_jacobi)
       !    mg%box_smoother => box_jacobi_lpl
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_lpl
       case default
          error stop "laplacian_set_methods: unsupported smoother type"
       end select
#if NDIM != 1
    case (mg_cylindrical)
       mg%box_op => box_clpl

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_clpl
       case default
          error stop "laplacian_set_methods: unsupported smoother type"
       end select
#endif
    case (mg_spherical)
       mg%box_op => box_slpl

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_slpl
       case default
          error stop "laplacian_set_methods: unsupported smoother type"
       end select
    case default
       error stop "laplacian_set_methods: unsupported geometry"
    end select

  end subroutine laplacian_set_methods

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine box_gs_lpl(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    real(dp)                  :: idr2(NDIM), fac
    logical                   :: redblack
#if NDIM == 3
    real(dp), parameter       :: sixth = 1/6.0_dp
#endif

    idr2 = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    fac = 0.5_dp / sum(idr2)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 1
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         cc(i, n) = fac * ( &
              idr2(1) * (cc(i+1, n) + cc(i-1, n)) - &
              cc(i, mg_irhs))
      end do
#elif NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di
            cc(i, j, n) = fac * ( &
                 idr2(1) * (cc(i+1, j, n) + cc(i-1, j, n)) + &
                 idr2(2) * (cc(i, j+1, n) + cc(i, j-1, n)) - &
                 cc(i, j, mg_irhs))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, di
               cc(i, j, k, n) = fac * ( &
                    idr2(1) * (cc(i+1, j, k, n) + cc(i-1, j, k, n)) + &
                    idr2(2) * (cc(i, j+1, k, n) + cc(i, j-1, k, n)) + &
                    idr2(3) * (cc(i, j, k+1, n) + cc(i, j, k-1, n)) - &
                    cc(i, j, k, mg_irhs))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_lpl

!   !> Perform Jacobi relaxation on box for a Laplacian operator
!   subroutine box_jacobi_lpl(mg, id, nc, cntr)
!     type(mg_t), intent(inout) :: mg
!     integer, intent(in)       :: id
!     integer, intent(in)       :: nc
!     integer, intent(in)       :: cntr !< Not used
!     integer                   :: IJK
!     real(dp), parameter       :: w     = 2.0_dp / 3
!     real(dp)                  :: tmp(DTIMES(0:nc+1))
!     real(dp)                  :: dr2
! #if NDIM == 3
!     real(dp), parameter       :: sixth = 1/6.0_dp
! #endif

!     dr2   = mg%dr(mg%boxes(id)%lvl)**2

!     associate (box => mg%boxes(id))
!       tmp = box%cc(DTIMES(:), mg_iphi)
!       do KJI_DO(1, nc)
! #if NDIM == 2
!          box%cc(i, j, mg_iphi) = (1-w) * box%cc(i, j, mg_iphi) + &
!               0.25_dp * w * ( &
!               tmp(i+1, j) + tmp(i-1, j) + &
!               tmp(i, j+1) + tmp(i, j-1) - &
!               dr2 * box%cc(i, j, mg_irhs))
! #elif NDIM == 3
!          box%cc(i, j, k, mg_iphi) = (1-w) * &
!               box%cc(i, j, k, mg_iphi) + &
!               sixth * w * ( &
!               tmp(i+1, j, k) + tmp(i-1, j, k) + &
!               tmp(i, j+1, k) + tmp(i, j-1, k) + &
!               tmp(i, j, k+1) + tmp(i, j, k-1) - &
!               dr2 * box%cc(i, j, k, mg_irhs))
! #endif
!       end do; CLOSE_DO
!     end associate
!   end subroutine box_jacobi_lpl

  !> Perform Laplacian operator on a box
  subroutine box_lpl(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK
    real(dp)                  :: idr2(NDIM)

    idr2 = 1 / mg%dr(:, mg%boxes(id)%lvl)**2

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 1
      do i = 1, nc
            cc(i, i_out) = &
                 idr2(1) * (cc(i-1, n) + cc(i+1, n) - 2 * cc(i, n))
         end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc
            cc(i, j, i_out) = &
                 idr2(1) * (cc(i-1, j, n) + cc(i+1, j, n) - 2 * cc(i, j, n)) + &
                 idr2(2) * (cc(i, j-1, n) + cc(i, j+1, n) - 2 * cc(i, j, n))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               cc(i, j, k, i_out) = &
                    idr2(1) * (cc(i-1, j, k, n) + cc(i+1, j, k, n) &
                    - 2 * cc(i, j, k, n)) &
                    + idr2(2) * (cc(i, j-1, k, n) + cc(i, j+1, k, n) &
                    - 2 * cc(i, j, k, n)) &
                    + idr2(3) * (cc(i, j, k-1, n) + cc(i, j, k+1, n) &
                    - 2 * cc(i, j, k, n))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_lpl

  !> Perform Laplacian operator on a box in cylindrical geometry, using (r,z)
  !> and (r,phi,z) coordinates in 2D/3D.
  subroutine box_clpl(mg, id, nc, i_out)
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

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)
            pre_fac(1) = r_cc(i)
            pre_fac(2) = r_cc(i)**2
            call Laplacian(cc(i, j, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(cc(i, j, k, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_clpl

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cylindrical geometry. 
  subroutine box_gs_clpl(mg, id, nc, redblack_cntr)
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
    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)
            pre_fac(1) = r_cc(i)
            pre_fac(2) = r_cc(i)**2
            call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, j), 1)
   
            do i = i0, nc, di
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)

               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_clpl

  !> Perform Laplacian operator on a box in spherical geometry, using (r,theta,phi)
  subroutine box_slpl(mg, id, nc, i_out)
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
    sin_theta_cc    = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=1,nc)])
#endif

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 1
      do i = 1, nc
         u(-1:1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(cc(i, i_out), u, idr2, face_coeff_in = face_coeff)
      end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            face_coeff(0:1, 2) = sin_theta_face(j:j+1)
            pre_fac(1) = sin_theta_cc(j)
            call Laplacian(cc(i, j, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(cc(i, j, k, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_slpl

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> spherical geometry.
  subroutine box_gs_slpl(mg, id, nc, redblack_cntr)
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
    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 1
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         u(-1:1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff)
         cc(i, n) = cc(i, n) - ( Lop - cc(i, mg_irhs) ) / dLdu
      end do

#elif NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            face_coeff(0:1, 2) = sin_theta_face(j:j+1)
            pre_fac(1) = sin_theta_cc(j)
            call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, j), 1)
   
            do i = i0, nc, di
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_slpl

end module m_laplacian
