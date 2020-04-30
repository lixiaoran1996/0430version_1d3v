module m_lookup_table
  use m_types
  use m_error

  implicit none
  private

  type LT_type
     private
     logical               :: linear_table
     integer               :: n_rows, n_vars
     real(dp)              :: scaling, inv_scaling
     real(dp)              :: x_min_invpow, inv_xrange_inv
     real(dp), allocatable :: x_values(:), y_values(:,:)
  end type LT_type

  type LT_loc_type
     private
     integer  :: low_ix
     real(dp) :: low_frac
  end type LT_loc_type

  interface LT_get
     module procedure LT_get_at_col, LT_get_at_loc_at_col, LT_get_normal, LT_get_at_loc, &
          LT_get_at_cols, LT_get_at_loc_at_cols, &
          LT_get_many_at_col, LT_get_many_at_locs_at_col, LT_get_many_normal, &
          LT_get_many_at_locs, LT_get_many_at_cols, LT_get_many_at_locs_at_cols
  end interface LT_get

  public :: LT_type
  public :: LT_loc_type

  public :: LT_create
  public :: LT_destroy
  public :: LT_add_column
  public :: LT_add_to_column
  public :: LT_get
  public :: LT_get_location
  public :: LT_get_locations
  public :: LT_get_num_vars
  public :: LT_get_num_rows
  public :: LT_get_all_data

contains

  function LT_create(x_min, x_max, n_rows, scale_fac, n_storage) result(my_lt)
    use m_utils
    real(dp), intent(in)           :: x_min, x_max
    real(dp), intent(in), optional :: scale_fac
    integer, intent(in)            :: n_rows
    integer, intent(in), optional  :: n_storage
    type(LT_type), pointer         :: my_lt

    integer                        :: n
    real(dp)                       :: frac, temp, scaling, inv_scaling

    if (present(scale_fac)) then
       scaling = scale_fac
    else
       scaling = 1.0_dp
    end if
    inv_scaling = 1.0_dp / scaling

    if (x_max <= x_min) call ERR_show("LT_create: x_max should be > x_min")
    if (n_rows <= 1)    call ERR_show("LT_create: n_rows should be bigger than 1")
    if (scaling <= 0)   call ERR_show("LT_create: scaling should be > 0")

    ! Create lookup table
    allocate(my_lt)
    allocate(my_lt%x_values(n_rows))

    if (present(n_storage)) then
       allocate(my_lt%y_values(n_storage, n_rows))
    else
       allocate(my_lt%y_values(1, n_rows))
    end if

    my_lt%n_rows         = n_rows
    my_lt%n_vars         = 0
    my_lt%x_min_invpow   = x_min**inv_scaling
    my_lt%inv_xrange_inv = 1 / (x_max**inv_scaling - x_min**inv_scaling)
    my_lt%scaling        = scaling
    my_lt%inv_scaling    = inv_scaling
    if (my_lt%scaling == 1.0_dp) my_lt%linear_table = .true.

    do n = 1, n_rows
       frac = (n-1) / (n_rows-1.0_dp)
       temp = frac * (x_max**inv_scaling - x_min**inv_scaling) + x_min**inv_scaling
       my_lt%x_values(n) = temp**scaling
    end do
  end function LT_create

  subroutine LT_destroy(my_lt)
    type(LT_type), pointer, intent(inout) :: my_lt
    deallocate(my_lt)
  end subroutine LT_destroy

  subroutine LT_add_column(my_lt, xx, yy)
    use m_utils
    type(LT_type), intent(inout) :: my_lt
    real(dp), intent(in)         :: xx(:), yy(:)
    integer                      :: n
    real(dp), allocatable        :: copy_y_values(:,:)

    if (size(xx) /= size(yy)) call ERR_show("LT_add_column: incompatible argument sizes")
    if (any(xx(2:size(xx)) < xx(1:size(xx)-1))) call ERR_show("LT_add_column: unsorted x-values")

    if (my_lt%n_vars >= size(my_lt%y_values, 1)) then
       allocate(copy_y_values(my_lt%n_vars, my_lt%n_rows))
       copy_y_values = my_lt%y_values
       deallocate(my_lt%y_values)
       allocate(my_lt%y_values(my_lt%n_vars+1, my_lt%n_rows))
       my_lt%y_values(1:my_lt%n_vars, :) = copy_y_values
       deallocate(copy_y_values)
    end if

    do n = 1, my_lt%n_rows
       call UT_linInterpList(xx, yy, my_lt%x_values(n), my_lt%y_values(my_lt%n_vars+1, n))
    end do

    my_lt%n_vars = my_lt%n_vars + 1
  end subroutine LT_add_column

  subroutine LT_add_to_column(my_lt, col_ix, xx, yy)
    use m_utils
    type(LT_type), intent(inout) :: my_lt
    integer, intent(in)          :: col_ix
    real(dp), intent(in)         :: xx(:), yy(:)
    integer                      :: n
    real(dp)                     :: temp

    if (size(xx) /= size(yy)) call ERR_show("LT_add_to_column: incompatible argument sizes")

    do n = 1, size(my_lt%x_values)
       call UT_linInterpList(xx, yy, my_lt%x_values(n), temp)
       my_lt%y_values(col_ix, n) = my_lt%y_values(col_ix, n) + temp
    end do
  end subroutine LT_add_to_column

  !> Sets a location object to point to the right place in the data
  subroutine LT_get_location(my_lt, x, my_loc)
    type(LT_type), intent(in)      :: my_lt
    real(dp), intent(in)           :: x
    type(LT_loc_type), intent(out) :: my_loc
    integer                        :: n_rows
    real(dp)                       :: frac, temp

    n_rows = my_lt%n_rows

    if (my_lt%linear_table) then
       frac   = (x - my_lt%x_min_invpow) * my_lt%inv_xrange_inv
    else
       frac   = (x**my_lt%inv_scaling - my_lt%x_min_invpow) * my_lt%inv_xrange_inv
    end if

    if (frac <= 0) then
       my_loc%low_ix = 1
       my_loc%low_frac = 1.0_dp
    else if (frac >= 1) then
       my_loc%low_ix = n_rows-1
       my_loc%low_frac = 0.0_dp
    else
       temp = 1 + frac * (n_rows-1)
       my_loc%low_ix = int(temp)
       my_loc%low_frac = 1.0_dp - temp + my_loc%low_ix
    end if

  end subroutine LT_get_location

    !> Sets a location object to point to the right place in the data
  subroutine LT_get_locations(my_lt, x_values, my_locs)
    type(LT_type), intent(in)      :: my_lt
    real(dp), intent(in)           :: x_values(:)
    type(LT_loc_type), intent(out) :: my_locs(:)
    integer                        :: n

    ! if (size(x_values) /= size(my_locs)) call ERR_show("LT_get*: incompatible argument sizes")

    do n = 1, size(x_values)
       call LT_get_location(my_lt, x_values(n), my_locs(n))
    end do
  end subroutine LT_get_locations

  subroutine LT_get_normal(my_lt, x, y_values)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(in)      :: x
    real(dp), intent(out)     :: y_values(:)
    type(LT_loc_type)         :: loc

    call LT_get_location(my_lt, x, loc)
    call LT_get_at_loc(my_lt, loc, y_values)
  end subroutine LT_get_normal

    subroutine LT_get_many_normal(my_lt, x_values, y_values)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(in)      :: x_values(:)
    real(dp), intent(out)     :: y_values(:,:)
    type(LT_loc_type)         :: locs(size(x_values))

    call LT_get_locations(my_lt, x_values, locs)
    call LT_get_many_at_locs(my_lt, locs, y_values)
  end subroutine LT_get_many_normal

  subroutine LT_get_at_loc(my_lt, loc, y_values)
    type(LT_type), intent(in)     :: my_lt
    type(LT_loc_type), intent(in) :: loc
    real(dp), intent(out)         :: y_values(:)

    if (size(y_values) /= my_lt%n_vars) call ERR_show("LT_get*: incompatible argument sizes")
    y_values = loc%low_frac * my_lt%y_values(1:my_lt%n_vars, loc%low_ix) &
         + (1-loc%low_frac) * my_lt%y_values(1:my_lt%n_vars, loc%low_ix+1)
  end subroutine LT_get_at_loc

  subroutine LT_get_many_at_locs(my_lt, locs, y_values)
    type(LT_type), intent(in)     :: my_lt
    type(LT_loc_type), intent(in) :: locs(:)
    real(dp), intent(out)         :: y_values(:,:)
    integer :: n

    if (size(y_values, 2) /= size(locs)) call ERR_show("LT_get*: incompatible argument sizes")

    do n = 1, size(locs)
       call LT_get_at_loc(my_lt, locs(n), y_values(:, n))
    end do
  end subroutine LT_get_many_at_locs

  subroutine LT_get_at_col(my_lt, x, col_ix, y_value)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(in)      :: x
    integer, intent(in)       :: col_ix
    real(dp), intent(out)     :: y_value
    type(LT_loc_type)         :: loc

    call LT_get_location(my_lt, x, loc)
    call LT_get_at_loc_at_col(my_lt, loc, col_ix, y_value)
  end subroutine LT_get_at_col

  subroutine LT_get_many_at_col(my_lt, x_values, col_ix, y_values)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(in)      :: x_values(:)
    integer, intent(in)       :: col_ix
    real(dp), intent(out)     :: y_values(:)
    type(LT_loc_type)         :: locs(size(x_values))

    call LT_get_locations(my_lt, x_values, locs)
    call LT_get_many_at_locs_at_col(my_lt, locs, col_ix, y_values)
  end subroutine LT_get_many_at_col

  subroutine LT_get_at_cols(my_lt, x, col_ixs, y_values)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(in)      :: x
    integer, intent(in)       :: col_ixs(:)
    real(dp), intent(out)     :: y_values(:)
    type(LT_loc_type)         :: loc

    call LT_get_location(my_lt, x, loc)
    call LT_get_at_loc_at_cols(my_lt, loc, col_ixs, y_values)
  end subroutine LT_get_at_cols

  subroutine LT_get_many_at_cols(my_lt, x_values, col_ixs, y_values)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(in)      :: x_values(:)
    integer, intent(in)       :: col_ixs(:)
    real(dp), intent(out)     :: y_values(:,:)
    type(LT_loc_type)         :: locs(size(x_values))

    call LT_get_locations(my_lt, x_values, locs)
    call LT_get_many_at_locs_at_cols(my_lt, locs, col_ixs, y_values)
  end subroutine LT_get_many_at_cols

  subroutine LT_get_at_loc_at_col(my_lt, loc, col_ix, y_value)
    type(LT_type), intent(in)     :: my_lt
    type(LT_loc_type), intent(in) :: loc
    integer, intent(in)           :: col_ix
    real(dp), intent(out)         :: y_value

    y_value = loc%low_frac * my_lt%y_values(col_ix, loc%low_ix) &
         + (1-loc%low_frac) * my_lt%y_values(col_ix, loc%low_ix+1)
  end subroutine LT_get_at_loc_at_col

  subroutine LT_get_many_at_locs_at_col(my_lt, locs, col_ix, y_values)
    type(LT_type), intent(in)     :: my_lt
    type(LT_loc_type), intent(in) :: locs(:)
    integer, intent(in)           :: col_ix
    real(dp), intent(out)         :: y_values(:)
    integer                       :: n

    if (size(y_values) /= size(locs)) call ERR_show("LT_get*: incompatible argument sizes")

    do n = 1, size(locs)
       call LT_get_at_loc_at_col(my_lt, locs(n), col_ix, y_values(n))
    end do
  end subroutine LT_get_many_at_locs_at_col

  subroutine LT_get_at_loc_at_cols(my_lt, loc, col_ixs, y_values)
    type(LT_type), intent(in)     :: my_lt
    type(LT_loc_type), intent(in) :: loc
    integer, intent(in)           :: col_ixs(:)
    real(dp), intent(out)         :: y_values(:)

    if (size(y_values) /= size(col_ixs)) call ERR_show("LT_get*: incompatible argument sizes")

    y_values = loc%low_frac * my_lt%y_values(col_ixs, loc%low_ix) &
         + (1-loc%low_frac) * my_lt%y_values(col_ixs, loc%low_ix+1)
  end subroutine LT_get_at_loc_at_cols

  subroutine LT_get_many_at_locs_at_cols(my_lt, locs, col_ixs, y_values)
    type(LT_type), intent(in)     :: my_lt
    type(LT_loc_type), intent(in) :: locs(:)
    integer, intent(in)           :: col_ixs(:)
    real(dp), intent(out)         :: y_values(:,:)
    integer                       :: n

    if (size(y_values, 2) /= size(locs)) call ERR_show("LT_get*: incompatible argument sizes")

    do n = 1, size(locs)
       call LT_get_at_loc_at_cols(my_lt, locs(n), col_ixs, y_values(:, n))
    end do
  end subroutine LT_get_many_at_locs_at_cols

  integer function LT_get_num_vars(my_lt)
    type(LT_type), intent(in) :: my_lt
    LT_get_num_vars = my_lt%n_vars
  end function LT_get_num_vars

  integer function LT_get_num_rows(my_lt)
    type(LT_type), intent(in) :: my_lt
    LT_get_num_rows = my_lt%n_rows
  end function LT_get_num_rows

  subroutine LT_get_all_data(my_lt, x_data, y_data)
    type(LT_type), intent(in) :: my_lt
    real(dp), intent(out) :: x_data(:), y_data(:,:)

    if (any(shape(y_data) /= (/my_lt%n_vars, my_lt%n_rows/)) .or. size(x_data) /= size(my_lt%x_values)) &
       call ERR_show("LT_get_all_data: incompatible argument sizes")

    x_data = my_lt%x_values
    y_data = my_lt%y_values(1:my_lt%n_vars, :)
  end subroutine LT_get_all_data

end module m_lookup_table
