module m_write_grid
   use m_types

   implicit none
   private

   include 'silo.inc'
   integer, parameter :: WG_DB_TYPE = DB_PDB
   integer, parameter :: WG_cartesian = 1
   integer, parameter :: WG_cylindrical = 2

   type WG_file_t
      private
      logical                 :: is_open
      integer                 :: db_ix
      character(len=line_len) :: name
   end type WG_file_t

   ! Public types
   public :: WG_file_t
   public :: WG_cartesian
   public :: WG_cylindrical

   ! Public methods
   public :: WG_create_file
   public :: WG_close_file
   public :: WG_create_dir
   public :: WG_change_dir
   public :: WG_add_grid_2d
   public :: WG_add_data_2d

contains

   function WG_create_file(filename) result(my_file)
      use m_error
      character(len=*), intent(in) :: filename
      type(WG_file_t)              :: my_file
      integer                      :: ierr, db_ix
      character(len=line_len)      :: fileinfo

      fileinfo = "A silo file created by WG_create_file"

      ! Create the Silo file
      ierr = dbcreate(trim(filename), len_trim(filename), DB_CLOBBER, DB_LOCAL, &
           fileinfo, len_trim(fileinfo), WG_DB_TYPE, db_ix)

      if (ierr /= 0) call ERR_show("WG_create_file error: cannot create " // trim(filename))

      my_file%name    = filename
      my_file%is_open = .true.
      my_file%db_ix   = db_ix

   end function WG_create_file

   subroutine WG_close_file(my_file)
      use m_error
      type(WG_file_t), intent(inout) :: my_file
      integer                        :: ierr

      if (my_file%is_open) then
         ierr = dbclose(my_file%db_ix)
         if (ierr /= 0) call ERR_show("WG_close_file error: cannot close " // trim(my_file%name))
         my_file%is_open = .false.
      else
         call ERR_warn("WG_close_file warning: file not open: " // trim(my_file%name))
      end if
   end subroutine WG_close_file

   subroutine WG_create_dir(my_file, dirname)
      use m_error
      type(WG_file_t), intent(in) :: my_file
      character(len=*), intent(in) :: dirname
      integer :: ierr, iostat

      ierr = dbmkdir(my_file%db_ix, trim(dirname), len_trim(dirname), iostat)
      if (ierr /= 0) call ERR_show("WG_create_dir error: cannot create " // trim(dirname))
   end subroutine WG_create_dir

   subroutine WG_change_dir(my_file, dirname)
      use m_error
      type(WG_file_t), intent(in)  :: my_file
      character(len=*), intent(in) :: dirname
      integer                      :: ierr
      ierr = dbsetdir(my_file%db_ix, trim(dirname), len_trim(dirname))
      if (ierr /= 0) call ERR_show("WG_change_dir: unable to go to directory " // trim(dirname))
   end subroutine WG_change_dir

   subroutine WG_add_grid_2d(my_file, grid_name, grid_type, grid_shape, pos_min, pos_max, curr_time)
      use m_error
      type(WG_file_t), intent(in)    :: my_file
      character(len=*), intent(in)   :: grid_name
      integer, intent(in)            :: grid_shape(2), grid_type
      real(dp), intent(in)           :: pos_min(2), pos_max(2)
      real(dp), intent(in), optional :: curr_time

      real(dp)                       :: time, deltas(2), dummy(1)
      real(dp)                       :: x_coords(grid_shape(1)), y_coords(grid_shape(2))
      integer                        :: i, ierr, iostat, dboptIx, db_ix

      interface
         function dbputqm(dbid, name, lname, xname, lxname, yname, &
                                  lyname, zname, lzname, x, y, z, dims, ndims, &
                                  datatype, coordtype, optlist_id, status)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbid, lname, lxname, lyname, lzname, dims(*), ndims
            integer(c_int) :: datatype, coordtype, optlist_id, status, dbputqm
            real(c_double) :: x(*), y(*), z(*)
            character(kind=c_char) :: name(*), xname(*), yname(*), zname(*)
         end function dbputqm
      end interface

      if (.not. my_file%is_open) call ERR_show("WG_add_grid_2d error: file not open " // trim(my_file%name))

      if ( present(curr_time) ) then
         time = curr_time
      else
         time = 0.0_dp
      end if

      deltas = (pos_max - pos_min) / dble(grid_shape-1)

      do i = 1, grid_shape(1)
         x_coords(i) = pos_min(1) + (i-1) * deltas(1)
      end do

      do i = 1, grid_shape(2)
         y_coords(i) = pos_min(2) + (i-1) * deltas(2)
      end do

      ! Make option list
      ierr = dbmkoptlist(20, dboptIx)

      select case (grid_type)
      case (WG_cartesian)
         ierr = dbaddcopt(dboptIx, DBOPT_XLABEL, 'x', 1)
         ierr = dbaddcopt(dboptIx, DBOPT_YLABEL, 'y', 1)
         ierr = dbaddiopt(dboptIx, DBOPT_COORDSYS, DB_CARTESIAN)
      case (WG_cylindrical)
         ierr = dbaddcopt(dboptIx, DBOPT_XLABEL, 'r', 1)
         ierr = dbaddcopt(dboptIx, DBOPT_YLABEL, 'z', 1)
         ierr = dbaddiopt(dboptIx, DBOPT_COORDSYS, DB_CYLINDRICAL)
      end select

      ierr  = dbaddiopt(dboptIx, DBOPT_NSPACE, 2)
      ierr  = dbaddiopt(dboptIx, DBOPT_HIDE_FROM_GUI, 0)
      ierr  = dbadddopt(dboptIx, DBOPT_DTIME, time)
      db_ix = my_file%db_ix

      ! Write the grid structure
      select case (grid_type)
      case (WG_cartesian)
         ierr = dbputqm(db_ix, trim(grid_name), len_trim(grid_name), 'x', 1, 'y', 1, '*', 1, x_coords, y_coords, dummy, &
              grid_shape, 2, DB_DOUBLE, DB_COLLINEAR, dboptIx, iostat)
      case (WG_cylindrical)
         ierr = dbputqm(db_ix, trim(grid_name), len_trim(grid_name), 'r', 1, 'z', 1, '*', 1, x_coords, y_coords, dummy, &
              grid_shape, 2, DB_DOUBLE, DB_COLLINEAR, dboptIx, iostat)
      end select

      ierr = dbfreeoptlist(dboptIx)
      if (ierr /= 0) call ERR_show("WG_add_grid_2d error: cannot free option list")
   end subroutine WG_add_grid_2d

   subroutine WG_add_data_2d(my_file, data_name, grid_name, my_data_2d, data_unit)
      use m_error
      type(WG_file_t), intent(in)  :: my_file
      character(len=*), intent(in) :: grid_name, data_name, data_unit
      real(dp), intent(in)         :: my_data_2d(:,:)
      integer                      :: dboptIx, ierr, iostat, grid_shape(2), db_ix
      real(dp)                     :: dummy(1)

      interface
         function dbputqv1(dbid, name, lname, meshname, lmeshname, &
                                   var, dims, ndims, mixvar, mixlen, datatype, &
                                   centering, optlist_id, status)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbid, lname, lmeshname, dims(*), ndims, mixlen
            integer(c_int) :: centering, optlist_id, status, datatype, dbputqv1
            real(c_double) :: var(*), mixvar(*)
            character(kind=c_char) :: name(*), meshname(*)
         end function dbputqv1
      end interface

      if (.not. my_file%is_open) call ERR_show("WG_add_data_2d error: file not open " // trim(my_file%name))

      grid_shape = shape(my_data_2d)

      ! Set options
      ierr = dbmkoptlist(10, dboptIx)
      ierr = dbaddcopt(dboptIx, DBOPT_UNITS, trim(data_unit), len_trim(data_unit))
      db_ix = my_file%db_ix

      ! Write the data to the grid
      ierr = dbputqv1(db_ix, trim(data_name), len_trim(data_name), trim(grid_name), len_trim(grid_name), &
           my_data_2d, grid_shape, 2, dummy, 0, DB_DOUBLE, DB_NODECENT, dboptIx, iostat)

      ierr = dbfreeoptlist(dboptIx)
      if (ierr /= 0) call ERR_show("WG_add_data_2d error: cannot free option list")
   end subroutine WG_add_data_2d

end module m_write_grid
