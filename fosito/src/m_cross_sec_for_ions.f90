! Copyright 2005-2012, Chao Li, Anbang Sun, Jannis Teunissen
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Module that contains routines to read in cross section data from textfiles.
!!!*****************************************************!!!!!!!!!!!!!
!!!! Anbang wrote this code for reading ions' cross section
module m_cross_sec_for_ions
   use m_types

   implicit none
   private

   !> The type of cross section table
   type CS_ion_type
      real(dp), allocatable   :: en_cs(:,:)    ! Stores the energy vs cross sec table
      integer                 :: n_rows        ! Number of rows in the table
      integer                 :: col_type      ! Type of collision, defined below
      real(dp)                :: spec_value    ! Specific value that can be used
      real(dp)                :: min_energy    ! Minimum energy for the collision
      real(dp)                :: max_energy    ! Maximum energy for the collision
      character(LEN=tiny_len) :: gas_name      ! Name of the colliding neutral molecule
      character(LEN=line_len) :: description   ! Description of the collision
      character(LEN=line_len) :: comment       ! Additional comments
   end type CS_ion_type

   type(CS_ion_type), allocatable :: CS_table_for_ions(:)

   integer, parameter :: CS_ion_elastic_t = 1, CS_ion_backwardscattering_t = 2
   integer, parameter :: max_num_cols_per_gas = 50
   integer, parameter :: max_num_rows = 200
   integer            :: CS_num_colls_found = 0

   ! Public variables
   public :: CS_ion_type
   public :: CS_ion_elastic_t, CS_ion_backwardscattering_t

   ! Methods
   public :: CS_reset_for_ions
   public :: CS_get_cross_secs_for_ions
   public :: CS_read_file_for_ions
   public :: CS_write_summary_for_ions
   public :: CS_write_all_for_ions

contains

   ! Reset all the cross sections that have been found
   subroutine CS_reset_for_ions()
      deallocate(CS_table_for_ions)
      CS_num_colls_found = 0
   end subroutine CS_reset_for_ions

   ! Return a copy of the cross sections that have been found
   subroutine CS_get_cross_secs_for_ions(cross_secs)
      type(CS_ion_type), allocatable :: cross_secs(:)
      allocate(cross_secs(CS_num_colls_found))
      cross_secs = CS_table_for_ions(1:CS_num_colls_found)
   end subroutine CS_get_cross_secs_for_ions

   ! Search 'filename' for cross section data concerning 'gas_name'
   subroutine CS_read_file_for_ions(filename, gas_name, x_normalization, y_normalization, req_energy)
      use m_error

      character(LEN=*), intent(IN) :: gas_name, filename
      real(dp), intent(in)         :: x_normalization, y_normalization, req_energy
      integer                      :: n, cIx, nL, n_rows, col_type
      integer                      :: my_unit, io_state, len_gas_name
      character(LEN=name_len)      :: lineFMT
      character(LEN=line_len)      :: line, prev_line, err_string
      real(dp)                     :: tempArray(2, max_num_rows)
      real(dp)                     :: x_scaling, y_scaling

      my_unit      = 333
      nL           = 0 ! Set the number of lines to 0
      len_gas_name = len(trim(gas_name))

      ! Set the line format to read, only depends on line_len currently
      write(lineFMT, FMT = "(A,I0,A)") "(A", line_len, ")"

      ! Open 'filename' (with error checking)
      open(my_unit, FILE = filename, STATUS = "OLD", ACTION = "READ", ERR = 999, IOSTAT = io_state)

      ! Look for collision processes with the correct gas name in the file,
      ! which should look for example like:

      !     ATTACHMENT                    [description of the type of process, always in CAPS]
      !     H2O -> H2O^-                  [the gas name possibly followed by the result of the process]
      !     COMMENT: total attachment     [possibly comments]
      !     UPDATED: 2010-06-24 15:04:36
      !     SCALING: 1.0 1.0              [optionally scale factors for the columns]
      !     ------------------            [at least 5 dashes]
      !     xxx   xxx                     [cross section data in two column format]
      !     ...   ...
      !     xxx   xxx
      !     ------------------

      ! So if we find the gas name the previous line holds the type of collision, then
      ! there is possibly some extra information and between the dashes the actual cross
      ! sections are found.

      ! The outer DO loop, running until the end of the file is reached
      do
         ! Search for 'gas_name' in the file
         line = ' '
         do
            prev_line = line
            read(my_unit, FMT = lineFMT, ERR = 999, end = 666) line; nL = nL+1
            line = adjustl(line)
            if (line(1:len_gas_name) == gas_name) exit
         end do

         ! Check prev_line for the type of collision
         select case (prev_line)
         case ("ELASTIC", "MOMENTUM")
            col_type = CS_ion_elastic_t
         case ("BACKWARDSCATTERING")
            col_type = CS_ion_backwardscattering_t
         case ("COMMENT")
            cycle
         case DEFAULT
            write(err_string, *) "CS_read_file_for_ions warning: ignoring unknown process type for ", &
                 gas_name, " in ", filename, " at line ", nL
            call ERR_warn(err_string)
            cycle
         end select

         ! Update the number of processes and set the gas name and collision type
         CS_num_colls_found = CS_num_colls_found + 1
         cIx = CS_num_colls_found
         call ensure_free_storage_for_ions(cIx)

         ! Add the reaction description to the table
         CS_table_for_ions(cIx)%description = adjustl(line)

!          ! For all collisions except attachment, there is a specific value on the next line
!          if (col_type /= CS_attach_t) then
!             read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
!             read(line, FMT = *, ERR = 999, end = 555) CS_table_for_ions(cIx)%spec_value
!          else
!             CS_table_for_ions(cIx)%spec_value = 0.0_dp
!          end if

         ! Read the spec_value underline
          read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
          read(line, FMT = *, ERR = 999, end = 555) CS_table_for_ions(cIx)%spec_value

         CS_table_for_ions(cIx)%gas_name = gas_name
         CS_table_for_ions(cIx)%col_type = col_type
         CS_table_for_ions(cIx)%comment = "COMMENT: (empty)"
         x_scaling = 1.0_dp
         y_scaling = 1.0_dp

         ! Now we can check whether there is a ZDPLASKIN statement and a reaction description,
         ! while scanning lines until dashes are found, which indicate the start of the cross section data
         do
            read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:9) == "ZDPLASKIN" ) then
               CS_table_for_ions(cIx)%description = trim(gas_name) // " [" // trim(adjustl(line(11:))) // "]"
            else if ( line(1:7) == "COMMENT") then
               CS_table_for_ions(cIx)%comment = line
            else if ( line(1:7) == "SCALING") then
               read(line(9:), *) x_scaling, y_scaling
            else if ( line(1:5) == "-----" ) then
               exit
            end if
         end do

         ! Read the cross section data into a temporary array
         n_rows = 0
         do
            read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:5) == "-----" ) then
               exit  ! Dashes mark the end of the data
            else if (trim(line) == "" .or. line(1:1) == "#") then
               cycle ! Ignore whitespace or comments
            else if (n_rows < max_num_rows) then
               n_rows = n_rows + 1
               read(line, FMT = *, ERR = 999, end = 555) tempArray(:, n_rows)
            else
               write(err_string, *) "CS_read_file_for_ions error: too many rows in ", filename, " at line ", nL
               call ERR_show(err_string)
            end if
         end do

         if (n_rows < 2) then
            write(err_string, *) "CS_read_file_for_ions error: need at least two values in ", &
              filename, " at line number ", nL
            call ERR_show(err_string)
         end if

         ! Store the data in the actual table
         allocate( CS_table_for_ions(cIx)%en_cs(2, n_rows) )
         CS_table_for_ions(cIx)%n_rows = n_rows
         CS_table_for_ions(cIx)%en_cs(1,:) = tempArray(1, 1:n_rows) * x_normalization * x_scaling
         CS_table_for_ions(cIx)%en_cs(2,:) = tempArray(2, 1:n_rows) * y_normalization * y_scaling
         CS_table_for_ions(cIx)%max_energy = CS_table_for_ions(cIx)%en_cs(1, n_rows)

         ! Locate minimum energy (first value followed by non-zero cross sec)
         do n = 1, n_rows-1
            if (CS_table_for_ions(cIx)%en_cs(2, n+1) > 0.0_dp) then
               CS_table_for_ions(cIx)%min_energy = CS_table_for_ions(cIx)%en_cs(1, n)
               exit
            end if
         end do

         ! Locate maximum energy (last value preceded by non-zero)
         do n = n_rows, 2, -1
            if (CS_table_for_ions(cIx)%en_cs(2, n-1) > 0.0_dp) then
               CS_table_for_ions(cIx)%max_energy = CS_table_for_ions(cIx)%en_cs(1, n)
               exit
            end if
         end do

         ! Check whether the tables that have been read in go up to high enough energies for our
         ! simulation. They can also have 0.0 as their highest listed cross section, in which
         ! case we assume the cross section is 0.0 for all higher energies

         if ( CS_table_for_ions(cIx)%en_cs(1, n_rows) < req_energy .and. &
              & CS_table_for_ions(cIx)%en_cs(2, n_rows) > 0.0D0 ) then
            write(err_string, *) "CS_read_file_for_ions error: cross section data at line ", nL, &
                 " does not go up to high enough x-values (energy)"
            call ERR_show(err_string)
         end if

      end do

555   continue ! Routine ends here if the end of "filename" is reached erroneously
      close(my_unit, ERR = 999, IOSTAT = io_state)
      write(err_string, *) "CS_read_file_for_ions error, reached end of file while searching. ", &
           "io_state = ", io_state, " while reading from [", filename, "] at line ", nL
      call ERR_show(err_string)
      return

666   continue ! Routine ends here if the end of "filename" is reached correctly
      close(my_unit, ERR = 999, IOSTAT = io_state)
      return

999   continue ! If there was an error, the routine will end here
      write(err_string, *) "CS_read_file_for_ions error at line ", nL, &
           " io_state = ", io_state, " while searching [", gas_name, "] in [", filename, "]"
      call ERR_show(err_string)

   end subroutine CS_read_file_for_ions

   subroutine ensure_free_storage_for_ions(req_size)
      use m_error

      integer, intent(in) :: req_size
      type(CS_ion_type), allocatable :: tbl_copy(:)
      integer :: curr_size

      if (allocated(CS_table_for_ions)) then
         curr_size = size(CS_table_for_ions)
         if (curr_size < req_size) then
            allocate(tbl_copy(curr_size))
            tbl_copy = CS_table_for_ions
            deallocate(CS_table_for_ions)
            allocate(CS_table_for_ions(max(req_size, 2 * curr_size)))
            CS_table_for_ions(1:curr_size) = tbl_copy
         end if
      else
         allocate(CS_table_for_ions(req_size))
      end if
   end subroutine ensure_free_storage_for_ions

   subroutine CS_write_all_for_ions(filename)
      use m_error
      character(LEN=*), intent(in) :: filename
      character(LEN=line_len)      :: err_string
      integer                      :: ics, n, io_state, my_unit

      my_unit = 333
      open(my_unit, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = io_state)

      write(my_unit, *, ERR = 999) "# A list of all the cross sections that have been read in by the"
      write(my_unit, *, ERR = 999) "# m_cross_sec_for_ions.f90 module. You can use this file as input again."
      write(my_unit, *, ERR = 999) ""

      do ics = 1, CS_num_colls_found
         select case (CS_table_for_ions(ics)%col_type)
         case (CS_ion_elastic_t)
            write(my_unit, *, ERR = 999) "ELASTIC"
            write(my_unit, *, ERR = 999) trim(CS_table_for_ions(ics)%description)
            write(my_unit, *, ERR = 999) CS_table_for_ions(ics)%spec_value, " / mass ratio"
         case (CS_ion_backwardscattering_t)
            write(my_unit, *, ERR = 999) "BACKWARDSCATTERING"
            write(my_unit, *, ERR = 999) trim(CS_table_for_ions(ics)%description)
            write(my_unit, *, ERR = 999) CS_table_for_ions(ics)%spec_value, " / mass ratio"
         end select

         write(my_unit, *, ERR = 999) trim(CS_table_for_ions(ics)%comment)
         write(my_unit, *, ERR = 999) "------------------------"
         do n = 1, CS_table_for_ions(ics)%n_rows
            write(my_unit, *, ERR = 999) CS_table_for_ions(ics)%en_cs(:, n)
         end do
         write(my_unit, *, ERR = 999) "------------------------"
         write(my_unit, *, ERR = 999) ""
      end do

      close(my_unit, ERR = 999, IOSTAT = io_state)
      return

999   continue ! If there was an error, the routine will end here
      write(err_string, *) "CS_write_all_for_ions error, io_state = ", io_state, " while writing to ", filename
      call ERR_show(err_string)

   end subroutine CS_write_all_for_ions

   subroutine CS_write_summary_for_ions(filename)
      use m_error
      character(LEN=*), intent(in) :: filename
      character(LEN=name_len)      :: col_name
      character(LEN=line_len)      :: err_string
      integer                      :: n, io_state, my_unit
      my_unit = 333

      open(my_unit, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = io_state)

      write(my_unit, ERR = 999, FMT = "(A)") "# List of collision processes"
      write(my_unit, ERR = 999, FMT = "(A)") "Index      Gasname            Coltype     Description"
      write(my_unit, ERR = 999, FMT = "(A)") "-----------------------------------------------------"

      do n = 1, CS_num_colls_found
         select case (CS_table_for_ions(n)%col_type)
         case (CS_ion_elastic_t)
            col_name = "Elastic"
         case (CS_ion_backwardscattering_t)
            col_name = "Backwardscattering"
         end select

         write(my_unit, ERR = 999, FMT = "((I4),(A),(A12),(A),(A15),(A),(A25))") &
              n, "    ", trim(CS_table_for_ions(n)%gas_name), "  ", trim(col_name), "     ", CS_table_for_ions(n)%description
      end do

      close(my_unit, STATUS = "KEEP", ERR = 999, IOSTAT = io_state)
      return

999   continue ! If there was an error, the routine will end here
      write(err_string, *) "CS_write_summary_for_ions error, io_state = ", io_state, " while writing to ", filename
      call ERR_show(err_string)

   end subroutine CS_write_summary_for_ions

end module m_cross_sec_for_ions
