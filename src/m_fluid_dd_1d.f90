!> A not so accurate implementation of a drift-diffusion-reaction plasma fluid model

module m_fluid_dd_1d
  use m_types
  use m_lookup_table
  use m_transport_schemes

  implicit none
  private

  integer               :: FL_num_vars
  integer, parameter    :: FL_iv_elec = 1, FL_iv_ion = 2, FL_iv_en = 3

  integer               :: FL_if_mob, FL_if_dif, FL_if_src, FL_if_en
  integer               :: FL_ie_mob, FL_ie_dif, FL_ie_src, FL_ie_loss
  logical               :: FL_use_en, FL_use_en_mob, FL_use_en_dif, FL_use_en_src

  integer               :: FL_grid_size
  real(dp), allocatable :: FL_vars(:,:)
  real(dp)              :: FL_delta_x
  real(dp)              :: FL_small_dens

  type(LT_type), pointer            :: FL_lkp_fld, FL_lkp_en
  procedure(TS_dd_1d_type), pointer :: FL_transport_scheme

  public :: FL_init_cfg
  public :: FL_advance
  public :: FL_get_output
  public :: FL_get_coeffs

contains

  subroutine FL_init_cfg()
    use m_transport_data
    use m_init_cond_1d
    use m_config
    use m_error

    integer, parameter     :: nameLen = 80
    integer                :: n, table_size
    real(dp)               :: max_energy, max_efield, xx
    real(dp), allocatable  :: x_data(:), y_data(:)

    FL_transport_scheme => TS_dd_1d_up_fl

    table_size   = CFG_get_int("fluid_lkptbl_size")
    max_energy   = CFG_get_real("fluid_lkptbl_max_energy")
    max_efield   = CFG_get_real("fluid_lkptbl_max_efield")
    FL_grid_size = CFG_get_int("grid_num_points")
    FL_delta_x   = CFG_get_real("grid_delta_x")
    FL_small_dens = CFG_get_real("fluid_small_density")
    FL_use_en = CFG_get_logic("fluid_use_energy")
    FL_use_en_mob = CFG_get_logic("fluid_use_en_mob")
    FL_use_en_dif = CFG_get_logic("fluid_use_en_dif")
    FL_use_en_src = CFG_get_logic("fluid_use_en_src")

    if (.not. FL_use_en .and. (FL_use_en_mob .or. FL_use_en_dif .or. FL_use_en_src)) &
         call ERR_show("FL_init_cfg() error: if fluid_use_energy is false, you cannot&
         & use energy coefficients: fluid_use_en_(mob/dif/src)!")

    if (allocated(FL_vars)) deallocate(FL_vars)
    if (FL_use_en) then
       FL_num_vars = 3
       allocate(FL_vars(FL_grid_size, FL_num_vars))
    else
       FL_num_vars = 2
       allocate(FL_vars(FL_grid_size, FL_num_vars))
    end if

    ! Initialization of electron and ion density
    do n = 1, FL_grid_size
       xx = (n-1) * FL_delta_x
       call INIT_get_elec_dens(xx, FL_vars(n, FL_iv_elec))
       call INIT_get_ion_dens(xx, FL_vars(n, FL_iv_ion))
       if (FL_use_en) call INIT_get_en_dens(xx, FL_vars(n, FL_iv_en))
    end do

    ! Create a lookup table for the model coefficients
    FL_lkp_fld => LT_create(0.0_dp, max_efield, table_size, scale_fac = 1.0_dp)
    if (FL_use_en) FL_lkp_en => LT_create(0.0_dp, max_energy, table_size, scale_fac = 1.0_dp)

    if (FL_use_en_mob) then
       call TD_get_data(trim(CFG_get_string("fluid_en_mob")), x_data, y_data)
       call LT_add_column(FL_lkp_en, x_data, y_data)
       FL_ie_mob = LT_get_num_vars(FL_lkp_en)
    else
       call TD_get_data(trim(CFG_get_string("fluid_fld_mob")), x_data, y_data)
       call LT_add_column(FL_lkp_fld, x_data, y_data)
       FL_if_mob = LT_get_num_vars(FL_lkp_fld)
    end if

    if (FL_use_en_dif) then
       call TD_get_data(trim(CFG_get_string("fluid_en_dif")), x_data, y_data)
       call LT_add_column(FL_lkp_en, x_data, y_data)
       FL_ie_dif = LT_get_num_vars(FL_lkp_en)
    else
       call TD_get_data(trim(CFG_get_string("fluid_fld_dif")), x_data, y_data)
       call LT_add_column(FL_lkp_fld, x_data, y_data)
       FL_if_dif = LT_get_num_vars(FL_lkp_fld)
    end if

    if (FL_use_en_src) then
       call TD_get_data(trim(CFG_get_string("fluid_en_alpha")), x_data, y_data)
       call LT_add_column(FL_lkp_en, x_data, y_data)
       FL_ie_src = LT_get_num_vars(FL_lkp_en)
       call TD_get_data(trim(CFG_get_string("fluid_en_eta")), x_data, y_data)
       call LT_add_to_column(FL_lkp_en, FL_ie_src, x_data, y_data)
    else
       call TD_get_data(trim(CFG_get_string("fluid_fld_alpha")), x_data, y_data)
       call LT_add_column(FL_lkp_fld, x_data, y_data)
       FL_if_src = LT_get_num_vars(FL_lkp_fld)
       call TD_get_data(trim(CFG_get_string("fluid_fld_eta")), x_data, y_data)
       call LT_add_to_column(FL_lkp_fld, FL_if_src, x_data, y_data)
    end if

    call TD_get_data(trim(CFG_get_string("fluid_fld_en")), x_data, y_data)
    call LT_add_column(FL_lkp_fld, x_data, y_data)
    FL_if_en = LT_get_num_vars(FL_lkp_fld)

    if (FL_use_en) then
       call TD_get_data(trim(CFG_get_string("fluid_en_loss")), x_data, y_data)
       call LT_add_column(FL_lkp_en, x_data, y_data)
       FL_ie_loss = LT_get_num_vars(FL_lkp_en)
    end if

  end subroutine FL_init_cfg

  subroutine FL_time_derivs(vars, time, time_derivs)
    use m_efield_1d
    use m_units_constants
    real(dp), intent(in)  :: vars(:,:), time
    real(dp), intent(out) :: time_derivs(:,:)

    integer               :: n_cc
    real(dp), parameter   :: five_third = 5 / 3.0_dp
    real(dp), parameter   :: mu_smooth  = 1.0_dp
    real(dp)              :: inv_delta_x
    real(dp)              :: mob_c(   size(vars, 1)-1)
    real(dp)              :: dif_c(   size(vars, 1)-1)
    real(dp)              :: src_c(   size(vars, 1)-1)
    real(dp)              :: flux(    size(vars, 1)-1)
    real(dp)              :: fld(     size(vars, 1)-1)
    real(dp)              :: fld_en(  size(vars, 1)-1)
    real(dp)              :: en_loss( size(vars, 1)-1)
    real(dp)              :: source(  size(vars, 1))
    real(dp)              :: mean_en( size(vars, 1)-1)
    type(LT_loc_type)     :: lt_locs( size(vars, 1)-1)

    n_cc = size(vars, 1)
    inv_delta_x = 1.0_dp / FL_delta_x

    ! Get electric field
    source = (vars(:, FL_iv_ion) - vars(:, FL_iv_elec)) * UC_elem_charge
    call EF_compute_and_get_st(source, fld,time, [0.d0,0.d0], 0, 0)

    ! Get efield-coefficients from the lookup table
    call LT_get_locations(FL_lkp_fld, abs(fld), lt_locs(1:n_cc-1))
    if (.not. FL_use_en_mob) call LT_get(FL_lkp_fld, lt_locs, FL_if_mob, mob_c)
    if (.not. FL_use_en_dif) call LT_get(FL_lkp_fld, lt_locs, FL_if_dif, dif_c)
    if (.not. FL_use_en_src) call LT_get(FL_lkp_fld, lt_locs, FL_if_src, src_c)

    if (FL_use_en) then
       ! There is a regularization: at density zero, we use the energy corresponding to the efield
       call LT_get(FL_lkp_fld, lt_locs, FL_if_en, fld_en)
       mean_en = (vars(1:n_cc-1, FL_iv_en) + vars(2:n_cc, FL_iv_en) + 2 * FL_small_dens * fld_en) &
            / (vars(1:n_cc-1, FL_iv_elec) + vars(2:n_cc, FL_iv_elec) + 2 * FL_small_dens)

       ! Get energy-coefficients from the lookup table
       call LT_get_locations(FL_lkp_en, mean_en, lt_locs)
       call LT_get(FL_lkp_en, lt_locs, FL_ie_loss, en_loss)
       if (FL_use_en_mob) call LT_get(FL_lkp_en, lt_locs, FL_ie_mob, mob_c)
       if (FL_use_en_dif) call LT_get(FL_lkp_en, lt_locs, FL_ie_dif, dif_c)
       if (FL_use_en_src) call LT_get(FL_lkp_en, lt_locs, FL_ie_src, src_c)
    end if

    ! ~~~ Electron transport ~~~
    call FL_transport_scheme(vars(:, FL_iv_elec), -mob_c * fld, dif_c * inv_delta_x, flux)

    ! ~~~ Electron source ~~~ TODO: try different formulations
    source(n_cc)     = 0.0_dp
    source(1:n_cc-1) = 0.5_dp * src_c * abs(flux)
    source(2:n_cc)   = source(2:n_cc) + source(1:n_cc-1)

    ! Set time derivatives
    time_derivs(:, FL_iv_elec) = source
    time_derivs(:, FL_iv_ion) = source
    call add_grad_flux_1d(time_derivs(:, FL_iv_elec), flux * inv_delta_x)

    if (FL_use_en) then
       ! ~~~ Energy source ~~~
       source(n_cc)     = 0.0_dp
       source(1:n_cc-1) = -0.5_dp * (fld * flux + en_loss * vars(1:n_cc-1, FL_iv_elec))
       source(2:n_cc)   = source(2:n_cc) - 0.5_dp * (fld * flux + en_loss * vars(2:n_cc, FL_iv_elec))

       ! ~~~ Energy transport ~~~
       call FL_transport_scheme(vars(:, FL_iv_en), -five_third * mob_c * fld, &
            five_third * dif_c * inv_delta_x, flux)

       ! Set time derivatives
       time_derivs(:, FL_iv_en) = source
       call add_grad_flux_1d(time_derivs(:, FL_iv_en), flux * inv_delta_x)
    end if

    ! Fix values at endpoints
    time_derivs(1, :) = 0.0_dp
    time_derivs(n_cc, :) = 0.0_dp
  end subroutine FL_time_derivs

  subroutine add_grad_flux_1d(dens_c, flux_f)
    real(dp), intent(inout) :: dens_c(:)
    real(dp), intent(in)    :: flux_f(:)
    integer                 :: n_cc

    n_cc             = size(dens_c)
    dens_c(1:n_cc-1) = dens_c(1:n_cc-1) - flux_f
    dens_c(2:n_cc)   = dens_c(2:n_cc) + flux_f
  end subroutine add_grad_flux_1d

  subroutine FL_advance(time, dt, max_dt, new_dt, abs_err)
    use m_time_steppers
    real(dp), intent(in)    :: dt, max_dt, abs_err
    real(dp), intent(out)   :: new_dt
    real(dp), intent(inout) :: time
    integer                 :: n
    real(dp)                :: max_errs(FL_grid_size, FL_num_vars)

    do n = 1, FL_num_vars
       max_errs(:, n) = abs_err * max(epsilon(1.0_dp), maxval(abs(FL_vars(:, n))))
    end do

    call STEP_rk4a_2d(FL_vars, max_errs, time, dt, max_dt, new_dt, FL_time_derivs)
  end subroutine FL_advance

  subroutine FL_get_output(pos_data, sca_data, data_names, n_pos, n_sca, time)
    use m_efield_1d
    use m_units_constants
    real(dp), intent(out), allocatable                :: pos_data(:,:), sca_data(:)
    character(len=name_len), intent(out), allocatable :: data_names(:)
    integer, intent(out)                              :: n_pos, n_sca
    integer                                           :: n
    real(dp)                                          :: temp_data(FL_grid_size), fld_en(FL_grid_size)
    real(dp)                                          :: potential(FL_grid_size)
    real(dp),intent(in)                               :: time

    n_pos = 6
    n_sca = 0
    allocate(pos_data(FL_grid_size, n_pos))
    allocate(sca_data(n_sca))
    allocate(data_names(n_pos+n_sca))

    do n = 1, FL_grid_size
       temp_data(n) = (n-1) * FL_delta_x
    end do

    data_names(1) = "position (m)"
    pos_data(:,1) = temp_data

    call EF_compute_and_get((FL_vars(:, FL_iv_ion) - FL_vars(:, FL_iv_elec)) * UC_elem_charge, potential, temp_data, &
                & time, [0.d0,0.d0], 0, 0)
    call LT_get(FL_lkp_fld, abs(temp_data), FL_if_en, fld_en)

    data_names(2) = "electric field (V/m)"
    pos_data(:,2) = temp_data

    data_names(3) = "electron density (1/m3)"
    pos_data(:,3) = FL_vars(:, FL_iv_elec)

    data_names(4) = "ion density (1/m3)"
    pos_data(:,4) = FL_vars(:, FL_iv_ion)

    data_names(5) = "energy density (1/m3)"
    data_names(6) = "mean energy (eV/m3)"

    if (FL_use_en) then
       pos_data(:,5) = FL_vars(:, FL_iv_en)
       pos_data(:,6) = (FL_vars(:, FL_iv_en) + FL_small_dens * fld_en) &
            / (FL_small_dens + FL_vars(:, FL_iv_elec))
    else
       pos_data(:,5) = FL_vars(:, FL_iv_elec) * fld_en
       pos_data(:,6) = fld_en
    end if

  end subroutine FL_get_output

  subroutine FL_get_coeffs(coeff_data, coeff_names, n_coeffs)
    real(dp), intent(out), allocatable :: coeff_data(:,:)
    character(len=name_len), intent(out), allocatable :: coeff_names(:)
    integer, intent(out) :: n_coeffs
    integer :: ix, n_rows, n_fld_coeffs, n_en_coeffs

    if (FL_use_en) then
       n_en_coeffs = 2 + count((/FL_use_en_mob, FL_use_en_dif, FL_use_en_src/))
       n_fld_coeffs = 7 - n_en_coeffs
    else
       n_en_coeffs = 0
       n_fld_coeffs = 5
    end if

    n_coeffs = n_fld_coeffs + n_en_coeffs
    n_rows = LT_get_num_rows(FL_lkp_fld)
    allocate(coeff_data(n_coeffs, n_rows))
    allocate(coeff_names(n_coeffs))

    call LT_get_all_data(FL_lkp_fld, coeff_data(1, :), coeff_data(2:n_fld_coeffs, :))
    coeff_names(1) = "efield (V/m)"
    coeff_names(1+FL_if_en) = "fld_en (eV)"
    if (.not. FL_use_en_mob) coeff_names(1+FL_if_mob) = "fld_mob (m2/Vs)"
    if (.not. FL_use_en_dif) coeff_names(1+FL_if_dif) = "fld_dif (m2/s)"
    if (.not. FL_use_en_src) coeff_names(1+FL_if_src) = "fld_src (1/m)"

    if (FL_use_en) then
       ix = n_fld_coeffs + 1
       call LT_get_all_data(FL_lkp_en, coeff_data(ix, :), coeff_data(ix+1:, :))
       coeff_names(ix) = "energy (eV/s)"
       coeff_names(ix+FL_ie_loss) = "en_loss (eV/s)"
       if (FL_use_en_mob) coeff_names(ix+FL_ie_mob) = "en_mob (m2/Vs)"
       if (FL_use_en_dif) coeff_names(ix+FL_ie_dif) = "en_dif (m2/s)"
       if (FL_use_en_src) coeff_names(ix+FL_ie_src) = "en_src (1/m)"
    end if

  end subroutine FL_get_coeffs

end module m_fluid_dd_1d
