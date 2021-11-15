!> Module to set boundary conditions from user data
module mod_bc_data
  use mod_lookup_table

  implicit none
  private

  integer, parameter :: max_boundary_conds = 10

  type bc_data_t
     integer                        :: n_variables
     double precision               :: origin(3)
     double precision               :: dx(3)
     integer                        :: n_points(3)
     character(len=40), allocatable :: names(:)
     double precision, allocatable  :: values(:, :, :, :)
  end type bc_data_t

  type(LT_t)  :: lt_1d(max_boundary_conds)
  type(LT2_t) :: lt_2d(max_boundary_conds)
  type(LT3_t) :: lt_3d(max_boundary_conds)

  !> Whether boundary condition data is time varying
  logical, public, protected :: bc_data_time_varying = .false.

  !> Integer array for indexing lookup tables per variable per direction
  integer, public, protected, allocatable :: bc_data_ix(:, :)

  public :: bc_data_init
  public :: bc_data_set
  public :: bc_data_get_2d
  public :: bc_data_get_3d

contains

  subroutine bc_data_init()
    use mod_global_parameters

    integer                :: i, iw, ib, n_files, n_bc
    character(len=std_len) :: bc_name, fname
    double precision       :: xmax(3)
    type(bc_data_t)        :: bc

    allocate(bc_data_ix(nprimbc, 2*ndim))

    bc_data_ix(:, :) = -1
    n_bc             = 0

    do ib = 1, 2 * ndim
       do iw = 1, nprimbc
          bc_name = typeboundary(iw, ib)
          if (bc_name(1:4) == "vtk:") then

             n_bc               = n_bc + 1
             fname              = bc_name(5:)
             bc_data_ix(iw, ib) = n_bc

             ! Other routines don't have to parse the full name
             typeboundary(iw, ib)   = "bc_data"

             call read_vtk_structured_points(trim(fname), bc)
             xmax = bc%origin + (bc%n_points-1) * bc%dx

             if (n_bc == 1) then
                bc_data_time_varying = (bc%n_points(ndim) > 1)
             else if (bc_data_time_varying .neqv. (bc%n_points(ndim) > 1)) &
                then
                call mpistop("bc_data_init: only some files are time varying")
             end if

             
             call mpistop("bc_data_init: 1D case not supported")
            
             
             
          end if
       end do
    end do

  end subroutine bc_data_init

  elemental function bc_data_get_3d(n_bc, x1, x2, qt) result(val)
    integer, intent(in)          :: n_bc
    double precision, intent(in) :: x1, x2, qt
    double precision             :: val

    if (bc_data_time_varying) then
       val = LT3_get_col(lt_3d(n_bc), 1, x1, x2, qt)
    else
       val = LT2_get_col(lt_2d(n_bc), 1, x1, x2)
    end if
  end function bc_data_get_3d

  elemental function bc_data_get_2d(n_bc, x1, qt) result(val)
    integer, intent(in)          :: n_bc
    double precision, intent(in) :: x1, qt
    double precision             :: val

    if (bc_data_time_varying) then
       val = LT2_get_col(lt_2d(n_bc), 1, x1, qt)
    else
       val = LT_get_col(lt_1d(n_bc), 1, x1)
    end if
  end function bc_data_get_2d

  !> Set boundary conditions according to user data
  subroutine bc_data_set(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nprim)
    double precision                :: tmp(ixOmin1:ixOmax1)
    integer                         :: i, ix, iw, n_bc

    

    

  end subroutine bc_data_set

  subroutine read_vtk_structured_points(fname, bc)
    character(len=*), intent(in)  :: fname
    type(bc_data_t), intent(out)  :: bc
    double precision, allocatable :: tmp_data(:, :, :, :)
    integer, parameter            :: max_variables = 10
    character(len=40)             :: tmp_names(max_variables)
    character(len=256)            :: line
    character(len=40)             :: word, typename
    integer, parameter            :: my_unit = 123
    integer                       :: n, n_points_total
    integer                       :: n_components

    open(my_unit, file=trim(fname), status="old", action="read")

    ! Header, e.g. # vtk DataFile Version 2.0
    read(my_unit, "(A)") line

    ! Dataset name
    read(my_unit, "(A)") line

    ! ASCII / BINARY
    read(my_unit, "(A)") line

    if (line /= "ASCII") then
       print *, "line: ", trim(line)
       error stop "read_vtk: not an ASCII file"
    end if

    ! DATASET STRUCTURED_POINTS
    read(my_unit, "(A)") line

    if (line /= "DATASET STRUCTURED_POINTS") then
       print *, "line: ", trim(line)
       error stop "read_vtk must have: DATASET STRUCTURED_POINTS"
    end if

    ! DIMENSIONS NX NY NZ
    read(my_unit, "(A)") line
    read(line, *) word, bc%n_points

    if (word /= "DIMENSIONS") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: DIMENSIONS"
    end if

    ! SPACING DX DY DZ
    read(my_unit, *) word, bc%dx

    if (word /= "SPACING") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: SPACING"
    end if

    ! ORIGIN 0 0 0
    read(my_unit, *) word, bc%origin
    if (word /= "ORIGIN") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: ORIGIN"
    end if

    ! POINT_DATA NPOINTS
    read(my_unit, *) word, n_points_total

    if (word /= "POINT_DATA") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: POINT_DATA n_points"
    end if

    if (n_points_total /= product(bc%n_points)) error stop "read_vtk: n_points not equal to NX*NY*NZ"

    allocate(tmp_data(bc%n_points(1), bc%n_points(2), bc%n_points(3),&
        max_variables))

    ! Read all scalar variables
    do n = 1, max_variables

       ! SCALARS name type ncomponents
       read(my_unit, *, end=900) word, tmp_names(n), typename, n_components

       if (word /= "SCALARS") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: SCALARS name type ncomponents"
       end if

       if (n_components /= 1) error stop "read_vtk: ncomponents should be 1"

       ! LOOKUP_TABLE default
       read(my_unit, *) word, typename

       if (word /= "LOOKUP_TABLE") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: LOOKUP_TABLE name"
       end if

       ! Read list of data values
       read(my_unit, *) tmp_data(:, :, :, n)
    end do

    ! Done reading variables
900 continue

    close(my_unit)

    if (n == max_variables + 1) error stop "read_vtk: increase max_variables"

    ! Loop index is one higher than number of variables
    bc%n_variables = n-1
    bc%values      = tmp_data(:, :, :, 1:n-1)
    bc%names       = tmp_names(1:n-1)
  end subroutine read_vtk_structured_points

end module mod_bc_data
