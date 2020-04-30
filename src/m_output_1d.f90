module m_output_1d
  use m_types
  use m_model_choice

  implicit none
  private

  integer  :: OUT_num_bins
  real(dp), allocatable :: OUT_eedf_eV_range(:)
  real(dp), allocatable :: OUT_eedf_min_max_fields(:)
  real(dp), allocatable :: ave_EEDF(:,:)

  public :: OUT_write_vars
  public :: OUT_write_coeffs
  public :: OUT_init
  public :: Cal_aver_EEDF_elec
  public :: OUT_write_aver_EEDF_elec
  
  public :: OUT_save_mid_sim_info, OUT_read_mid_sim_info

contains

    subroutine OUT_init()
        use m_config
        use m_error
        integer :: number, varSize

        OUT_num_bins = CFG_get_int("output_num_bins")

        varSize = CFG_get_size("output_eedf_min_max_fields")
        allocate(OUT_eedf_min_max_fields(varSize))
        do number = 1, varSize
            OUT_eedf_min_max_fields(number) = CFG_get_real("output_eedf_min_max_fields",number)
        end do

        varSize = CFG_get_size("output_eedf_eV_range")
        allocate(OUT_eedf_eV_range(varSize))
        do number = 1, varSize
            OUT_eedf_eV_range(number)  = CFG_get_real("output_eedf_eV_range",number)
        end do

        allocate(ave_EEDF(3, OUT_num_bins))
        ave_EEDF = 0.d0
    end subroutine OUT_init

  subroutine OUT_write_vars(sim_name, sim_type, cntr, time, myrank, root, ntasks)
    use m_fluid_dd_1d, only: FL_get_output
    use m_particle_1d, only: PM_get_output, PM_get_eedf,PM_collectEEDFAtRoot,PM_shareEEDF
    use m_error

    character(len=*), intent(in) :: sim_name
    integer, intent(in) :: sim_type, cntr
    real(dp), intent(in) :: time

    integer :: n_pos, n_sca
    character(len=name_len) :: filename,  eedf_names(3)

    real(dp), allocatable :: pos_data(:,:), sca_data(:), eedf(:,:)
    character(len=name_len), allocatable :: data_names(:)

    integer, intent(in)         :: myrank, root, ntasks

    write(filename, fmt="(I0,A)"), cntr, ".txt"

    select case (sim_type)
    case (MODEL_part)
       filename = "output/" // trim(sim_name) // "_part_" // trim(filename)
       call PM_get_output(pos_data, sca_data, data_names, n_pos, n_sca,time, myrank, root)
       if (n_pos > 0 .and. myrank == root) then
          call write_data_2d(filename, pos_data, data_names(1:n_pos), 30)
       end if

       ! write EEDF for electrons
       eedf_names(1) = "eV"
       eedf_names(2) = "eedf"
       eedf_names(3) = "eepf"

        call PM_get_eedf(eedf, OUT_eedf_eV_range, &
            OUT_eedf_min_max_fields, OUT_num_bins)

        !Collect and share
        call PM_collectEEDFAtRoot(eedf, myrank, root, ntasks, OUT_num_bins)
        call PM_shareEEDF(eedf, myrank, root, OUT_num_bins)
        write(filename, fmt="(A,I0,A)") &
            "output/" // trim(sim_name) // "_part_", cntr, "_eedf"//".txt"
        if (myrank == root) call write_data_2d(filename, eedf, eedf_names, &
            30, do_transpose = .true.)
        deallocate(eedf)
 

    case (MODEL_fluid)
       filename = "output/" // trim(sim_name) // "_fluid_" // trim(filename)
       call FL_get_output(pos_data, sca_data, data_names, n_pos, n_sca,time)

       if (n_pos > 0) then
          call write_data_2d(filename, pos_data, data_names(1:n_pos), 30)
       end if

    case default
       call ERR_show("output not implemented yet")
    end select

    if (myrank == root) print *, "Written output " // trim(filename) // " at t = ", time
  end subroutine OUT_write_vars

  ! here we calculate and output the average EEDF during a time, according to turner' code
  subroutine Cal_aver_EEDF_elec(stepsAver, myrank, root, ntasks)
      use m_particle_1d, only: PM_get_eedf,PM_collectEEDFAtRoot,PM_shareEEDF
      use m_error

      integer, intent(in)   :: stepsAver
      real(dp), allocatable :: eedf(:,:)
      integer, intent(in)   :: myrank, root, ntasks
          
       call PM_get_eedf(eedf, OUT_eedf_eV_range, &
            OUT_eedf_min_max_fields, OUT_num_bins)
        !Collect and share
        call PM_collectEEDFAtRoot(eedf, myrank, root, ntasks, OUT_num_bins)
        call PM_shareEEDF(eedf, myrank, root, OUT_num_bins)

       if (size(ave_EEDF) /= size(eedf)) then
            call ERR_show("Error in calculating the average EEDF for electrons in m_Output_1.f90!")
       end if

       ave_EEDF(1,:) = eedf(1,:)
       ave_EEDF(2,:) = ave_EEDF(2,:) + eedf(2,:) / stepsAver
       ave_EEDF(3,:) = ave_EEDF(3,:) + eedf(3,:) / stepsAver

       deallocate(eedf)

  end subroutine Cal_aver_EEDF_elec

  subroutine OUT_write_aver_EEDF_elec(sim_name, time)
        character(len=*), intent(in) :: sim_name
        real(dp), intent(in) :: time
        character(len=name_len) :: filename
        integer              :: myunit, nn
        
        myunit = 78956
        filename = "output/" // trim(sim_name) // "_part_" // "aver_EEDF.txt"
        open(unit = myunit, file = filename , status = 'unknown')
        write(myunit, *) "eV / eedf"
        do nn = 1, OUT_num_bins
            write(myunit, *) ave_EEDF(:,nn)
        end do
        close(myunit)

        print *, "output the average EEDF of electrons " // trim(filename) // " at t = ", time       
  end subroutine OUT_write_aver_EEDF_elec

  subroutine OUT_write_coeffs(sim_name, sim_type)
    use m_fluid_dd_1d, only: FL_get_coeffs
    use m_particle_1d, only: PM_get_coeffs
    use m_error

    character(len=*), intent(in) :: sim_name
    integer, intent(in) :: sim_type

    integer :: n_coeffs
    character(len=name_len) :: filename

    real(dp), allocatable :: coeff_data(:,:)
    character(len=name_len), allocatable :: coeff_names(:)

    select case (sim_type)
    case (MODEL_part)
       filename = "output/" // trim(sim_name) // "_part_coeffs.txt"
       call PM_get_coeffs(coeff_data, coeff_names, n_coeffs)

       if (n_coeffs > 0) then
          call write_data_2d(filename, coeff_data, coeff_names, 30, do_transpose = .true.)
       end if
    case (MODEL_fluid)
       filename = "output/" // trim(sim_name) // "_fluid_coeffs.txt"
       call FL_get_coeffs(coeff_data, coeff_names, n_coeffs)

       if (n_coeffs > 0) then
          call write_data_2d(filename, coeff_data, coeff_names, 30, do_transpose = .true.)
       end if

    case default
       call ERR_show("output not implemented yet")
    end select

    print *, "Written output " // trim(filename)
  end subroutine OUT_write_coeffs

  subroutine write_data_2d(filename, data_2d, col_names, col_width, do_transpose)
    use m_error
    character(len=*), intent(in)                           :: filename
    real(dp), intent(in) :: data_2d(:,:)
    integer, intent(in) :: col_width
    character(len=name_len), intent(in) :: col_names(:)
    logical, intent(in), optional :: do_transpose

    integer :: my_unit, n, n_rows, n_cols, io_state
    character(len=tiny_len) :: fmt_string
    real(dp), allocatable :: copy_of_data(:,:)
    logical :: transpose_data

    my_unit = 333
    if (present(do_transpose)) then
       transpose_data = do_transpose
    else
       transpose_data = .false.
    end if

    if (transpose_data) then
       n_rows = size(data_2d, 2)
       n_cols = size(data_2d, 1)
       allocate(copy_of_data(n_rows, n_cols))
       copy_of_data = transpose(data_2d)
    else
       n_rows = size(data_2d, 1)
       n_cols = size(data_2d, 2)
       allocate(copy_of_data(n_rows, n_cols))
       copy_of_data = data_2d
    end if

    if (size(col_names) /= n_cols) call ERR_show("write_data_2d: incompatible argument sizes")

    open(unit=my_unit, file=filename, iostat=io_state)
    if (io_state /= 0) call ERR_show("write_data_2d: error writing " // filename)

    ! Create format string for the header
    write(fmt_string, fmt="(A,I0,A)") "(A", col_width, ")"

    ! Write header
    do n = 1, n_cols
       write(my_unit, advance="NO", FMT=fmt_string) "  # " // col_names(n)
    end do

    write(my_unit, *) ""

    ! Create format string for data
    write(fmt_string, fmt="(A,I0,A,I0,A,I0,A)") "(", n_cols, "E", col_width, ".", col_width - 9, "e3)"

    ! Write data
    do n = 1, n_rows
       write(my_unit, fmt_string) copy_of_data(n, :)
    end do

    close(my_unit)

  end subroutine write_data_2d

! save and read the mid info during simulations
  subroutine OUT_save_mid_sim_info(sim_time, steps, stepIons, output_cnrt,  sim_name)
 
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer, intent(in)         :: steps, output_cnrt, stepIons
        real(dp), intent(in)         :: sim_time
        integer                      :: myunit

        myunit = 399
        filename = "midstore/" // trim(sim_name) // "_info_sim"// ".txt"
        open(unit = myunit, file = filename , status = 'unknown')
        write(myunit, *) sim_time, steps, stepIons, output_cnrt
        close(myunit)
  
  end subroutine OUT_save_mid_sim_info

  subroutine OUT_read_mid_sim_info(sim_time, steps, stepIons, output_cnrt, sim_name)
 
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer, intent(out)         ::  steps, output_cnrt, stepIons
        real(dp), intent(out)         :: sim_time
        integer                      :: myunit

        myunit = 399
        filename = "midstore/" // trim(sim_name) // "_info_sim"// ".txt"
        open(unit = myunit, file = filename , status = 'old')
        read(myunit, *) sim_time, steps, stepIons, output_cnrt
        close(myunit)
  
  end subroutine OUT_read_mid_sim_info

end module m_output_1d
