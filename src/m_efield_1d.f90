!> Module to compute electric fields in 1D
! In 1D the electric field is simply given by
! E_i = E_0 - Sum of charge before current point / epsilon0,

module m_efield_1d
  use m_types
  use m_utils

  implicit none
  private

  integer               :: EF_grid_size
  real(dp)              :: EF_delta_x, EF_inv_delta_x
  real(dp),allocatable  :: EF_applied_field(:),EF_efield_time(:),EF_applied_voltage(:)
  real(dp), allocatable :: EF_values(:), EP_values(:)
  logical               :: EF_is_constant
  logical               :: EF_use_voltage
  integer               :: voltage_type   ! decide AC or DC
  integer, parameter    :: voltage_DC = 0, voltage_AC = 1, voltage_pulse = 2
  real(dp)              :: voltage_AC_max, voltage_AC_fre   ! the amplitude and frequency of the capactive voltage 
  ! the amplitude,time  and frequency of the pulse voltage
  real(dp)              :: voltage_pulse_falling_edge,voltage_pulse_fre,voltage_pulse_max
  real(dp)              :: voltage_pulse_platform,voltage_pulse_rising_edge
  integer               :: EF_DBD_type
  real(dp)              :: EF_DBD_len(2)
  real(dp)              :: EF_epsDBD(2)   ! the permitivity of two dielectrics
  logical               :: EF_useDBD  !Whether we use DBD
  integer               :: EF_DBD_index(2)
  real(dp)              :: EF_old_field(4)

  public :: EF_initialize
  public :: EF_compute
  public :: EF_compute_and_get
  public :: EF_compute_and_get_st
  public :: EF_get_at_pos
  public :: EF_get_values
  public :: EF_get_values_st
  public :: EF_get_potential_at_bound
  public :: EF_get_DBD_index_and_type

  public :: EF_displace_currtent_at
  public :: EF_index_for_current
  public :: EF_get_old_efield
  public :: EF_get_boundary_potential
  public :: EF_get_efield_index
  public :: EF_get_values_dbd
  public :: EF_get_Gap_and_dielectric_potential


contains

  subroutine EF_get_values_dbd(out_efield)
    real(dp), intent(out) :: out_efield(2)
    out_efield(1) = EF_values(EF_DBD_index(1) + 1)
    out_efield(2) = EF_values(EF_DBD_index(2))
  end subroutine EF_get_values_dbd

  subroutine EF_initialize()
    use m_config
    use m_error
    use m_units_constants
    
    integer :: varSize
    integer :: number

    EF_grid_size      = CFG_get_int("grid_num_points")
    EF_delta_x     = CFG_get_real("grid_delta_x")
    EF_inv_delta_x = 1.0_dp / EF_delta_x
    
    EF_is_constant = CFG_get_logic("sim_constant_efield")
    EF_use_voltage = CFG_get_logic("sim_use_voltage")
    voltage_type   = CFG_get_int("sim_voltage_type")
    voltage_AC_max = CFG_get_real("sim_voltage_AC_max")
    voltage_AC_fre = CFG_get_real("sim_voltage_AC_fre") 
	voltage_pulse_max = CFG_get_real("sim_voltage_pulse_max ")
    voltage_pulse_fre = CFG_get_real("sim_voltage_pulse_fre ")
	voltage_pulse_rising_edge = CFG_get_real("sim_voltage_rising_edge ")
    voltage_pulse_falling_edge = CFG_get_real("sim_voltage_falling_edge ")
	voltage_pulse_platform = CFG_get_real("sim_voltage_platform ")
    EF_useDBD      = CFG_get_logic("sim_DBD")
 
    ! permitivity of dielectrics
    do number = 1, 2
        EF_epsDBD(number)  = CFG_get_real("sim_relPer_del",number) * UC_eps0
    end do

    !Anbang: check the length of dbd on both sides, then we can decide the boundary condition for possion solver with dielectrics
    do number = 1, 2
        EF_DBD_len(number)  = CFG_get_real("sim_len_del",number)
    end do

    ! decide the type of dbd: 1- only left; 2-only right, 3- both sides
    if (EF_DBD_len(1) > 0.d0 .and. EF_DBD_len(2) == 0.d0) then
        EF_DBD_type = 1
    else if (EF_DBD_len(1) == 0.d0 .and. EF_DBD_len(2) > 0.d0) then
        EF_DBD_type = 2
    else if (EF_DBD_len(1) > 0.d0 .and. EF_DBD_len(2) > 0.d0) then
        EF_DBD_type = 3
    else
        if (EF_useDBD) call ERR_show("EF_DBD_type is wrong!")
    end if

    

    !Anbang: get the index od the dbd position
    EF_DBD_index(1) = get_index_pos(EF_DBD_len(1))
    EF_DBD_index(2) = get_index_pos((EF_grid_size -1) * EF_delta_x - EF_DBD_len(2))
    
    if (EF_useDBD) then
        print *, "EF_DBD_index =", EF_DBD_index
        print *, "EF_DBD_type = ", EF_DBD_type
    end if

   ! initial old field at some positions
    EF_old_field = 0.d0
    

    ! Anbang: Here we decide to use voltage/efield as boundary conditions
    if ( .not. EF_use_voltage) then  ! if we use electric field as input data
        varSize = CFG_get_size("sim_applied_efield")
        allocate(EF_efield_time(varSize))
        allocate(EF_applied_field(varSize))
        if (size(EF_efield_time) /= size(EF_applied_field)) call ERR_show("EF_initialize(efield): argument has wrong size")

        do number = 1, varSize
            EF_applied_field(number)  = CFG_get_real("sim_applied_efield",number)
            EF_efield_time(number)    = CFG_get_real("sim_elec_times",number)
        end do

        ! The electric field is defined at cell faces, so it includes one extra point.
        ! e.g., if the charge density is defined at 0, dx, 2*dx, then the electric field is
        ! defined at -0.5*dx, 0.5*dx, 1.5*dx, 2.5*dx, with the first value equal to the applied field.
        ! These extra values are mostly useful as a buffer for EF_get_at_pos
        allocate(EF_values(EF_grid_size+1))
        EF_values = EF_applied_field(1)

        ! if we use efield boundary conditions, we simply set potential as zero
        allocate(EP_values(EF_grid_size))
        EP_values = 0.d0
     else  ! if we use voltage as boundary
        select case (voltage_type)
        case (voltage_DC)
            varSize = CFG_get_size("sim_applied_voltage")
            allocate(EF_efield_time(varSize))
            allocate(EF_applied_voltage(varSize))
            if (size(EF_efield_time) /= size(EF_applied_voltage)) call ERR_show("EF_initialize(voltage): argument has wrong size")

            do number = 1, varSize
                EF_applied_voltage(number)  = CFG_get_real("sim_applied_voltage",number)
                EF_efield_time(number)    = CFG_get_real("sim_elec_times",number)
            end do

            ! The electric field is defined at cell faces, so it includes one extra point.
            ! e.g., if the charge density is defined at 0, dx, 2*dx, then the electric field is
            ! defined at -0.5*dx, 0.5*dx, 1.5*dx, 2.5*dx, with the first value equal to the applied field.
            ! These extra values are mostly useful as a buffer for EF_get_at_pos
            allocate(EF_values(EF_grid_size+1))
            EF_values = -EF_applied_voltage(1) / (EF_delta_x * (EF_grid_size - 1))
            
            !left and right boundary: efield is zero
            EF_values(1) = 0.d0
            EF_values(EF_grid_size+1) = 0.d0

            ! if we use voltage boundary conditions, we only set potential as zero initially
            allocate(EP_values(EF_grid_size))
            EP_values = 0.d0
         case (voltage_AC)  ! here we simulate capactive discharges
            
            allocate(EF_values(EF_grid_size+1)) 
            ! here the voltage is: U = U0 * sin( 2 * pi * f * t)
            EF_values = - voltage_AC_max * sin (2.d0 * UC_pi * voltage_AC_fre * 0.d0) / (EF_delta_x * (EF_grid_size - 1))

            !left and right boundary: efield is zero
            EF_values(1) = 0.d0
            EF_values(EF_grid_size+1) = 0.d0

            ! if we use voltage boundary conditions, we only set potential as zero initially
            allocate(EP_values(EF_grid_size))
            EP_values = 0.d0
		case (voltage_pulse)
			allocate(EF_values(EF_grid_size+1)) 
            ! here the voltage is: U = U0 * sin( 2 * pi * f * t)
            EF_values = 0.d0

            !left and right boundary: efield is zero
            EF_values(1) = 0.d0
            EF_values(EF_grid_size+1) = 0.d0

            ! if we use voltage boundary conditions, we only set potential as zero initially
            allocate(EP_values(EF_grid_size))
            EP_values = 0.d0
            
         end select    
    end if

  end subroutine EF_initialize

  subroutine EF_compute(net_charge,time, surfCharge, myrank, root)
    use m_error
    use m_units_constants
    include "mpif.h"

    real(dp), intent(in) :: net_charge(:),surfCharge(:)
    real(dp)             :: conv_fac
    integer              :: iz
    real(dp),intent(in)  :: time
    integer, intent(in)  :: myrank, root
    integer              :: ierr
    real(dp), allocatable:: tempV(:), tempE(:)

    if (myrank == root) then   !only calculate on root, and then broadcast to all nodes
        if (size(net_charge) /= EF_grid_size) call ERR_show("EF_compute: argument has wrong size")

        if (.not. EF_use_voltage) then
            !>Anbang: Here we need to change the boundary electric field every time step
            !> Be careful here is only valid when the left boundary is far from the plasma region, 
            !>    later if the electron arrives at the boundary or there is charge accumulation,
            !  this boundary condition should be reconsidered
            EF_values(1) = EL_getfield(time)

            conv_fac = EF_delta_x / UC_eps0
            if (.not. EF_is_constant) then
                do iz = 2, EF_grid_size+1
                    EF_values(iz) = EF_values(iz-1) + net_charge(iz-1) * conv_fac
                end do
            end if
        else
            !> Here we use Dirichalet boundary conditions: left:V=V0, right: V=0
            if (EF_useDBD) then
                call PSol_DBD(net_charge, time, surfCharge, EP_values, EF_values)
            else
                call PSol_without_DBD(net_charge,time,EP_values,EF_values)
            end if
        end if
    end if

    allocate(tempV(EF_grid_size))
    allocate(tempE(EF_grid_size + 1))

    if (myrank == root) then
        tempV = EP_values
        tempE = EF_values
    end if
    ! broadcast the electric field and electric potential to all nodes
    call MPI_BCAST(tempV, EF_grid_size, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)
    call MPI_BCAST(tempE, EF_grid_size+1, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)

    if (myrank /= root) then
        EP_values = tempV
        EF_values = tempE
    end if

    deallocate(tempV, tempE)

  end subroutine EF_compute

  real(dp) function EL_getfield(time)
    real(dp),intent(in)  :: time

    call UT_linInterpList(EF_efield_time, EF_applied_field, time, EL_getfield)
  end function EL_getfield

  real(dp) function EL_getvoltage(time)
    real(dp),intent(in)  :: time

    call UT_linInterpList(EF_efield_time, EF_applied_voltage, time, EL_getvoltage)
  end function EL_getvoltage

  subroutine EF_compute_and_get_st(net_charge, out_efield,time,surfCharge, myrank, root)
    real(dp), intent(in) :: net_charge(:), surfCharge(:)
    real(dp), intent(out) :: out_efield(:)
    real(dp),intent(in)  :: time
    integer, intent(in)  :: myrank, root
    call EF_compute(net_charge,time, surfCharge, myrank, root)
    call EF_get_values_st(out_efield)
  end subroutine EF_compute_and_get_st

  subroutine EF_compute_and_get(net_charge, out_voltage, out_efield,time, surfCharge,myrank, root)
    real(dp), intent(in) :: net_charge(:), surfCharge(:)
    real(dp), intent(in) :: time
    real(dp), intent(out) :: out_efield(:), out_voltage(:)
    integer, intent(in)  :: myrank, root

    call EF_compute(net_charge,time, surfCharge,myrank, root)
    call EF_get_values(out_efield)
    out_voltage = EP_values
  end subroutine EF_compute_and_get

! Anbang: computer and get the displacement current at a certain position
 ! I_dis = (E_t+1 - E_t) / dt * sigma
  subroutine EF_displace_currtent_at(dt, dis_current, aera)
    use m_units_constants
    real(dp), intent(in) :: dt, aera
    real(dp)             :: out_efield(EF_grid_size)
    real(dp), intent(out),allocatable :: dis_current(:)

    integer               :: indx(4)
    integer               :: i

    call EF_index_for_current(indx)
    allocate(dis_current(size(indx)))
    call EF_get_values(out_efield)
    do i = 1, size(indx)
        dis_current(i) = (out_efield(indx(i)) - EF_old_field(i))/dt * UC_eps0 * aera
    end do
  end subroutine EF_displace_currtent_at

  ! anbang: pass efield to the old value
  subroutine EF_get_old_efield()

    real(dp)             :: out_efield(EF_grid_size)

    integer               :: indx(4)
    integer               :: i

    call EF_index_for_current(indx)
    call EF_get_values(out_efield)
    do i = 1, size(indx)
        EF_old_field(i) = out_efield(indx(i))
    end do
  end subroutine EF_get_old_efield

!get efield at indx points
  subroutine EF_get_efield_index(field)

    real(dp)             :: out_efield(EF_grid_size)

    integer               :: indx(4)
    integer               :: i
    real(dp),intent(out), allocatable  :: field(:)

    allocate(field(size(indx)))
    call EF_index_for_current(indx)
    
    call EF_get_values(out_efield)
    do i = 1, size(indx)
        field(i) = out_efield(indx(i))
    end do
  end subroutine

 ! anbang: index for positions where we output current 
  subroutine EF_index_for_current(indx)
        integer,intent(out)  :: indx(4)
    ! for testing, at different position   
        if (EF_useDBD) then
            indx(1) = EF_DBD_index(1) + 1
            indx(4) = EF_DBD_index(2) - 1
        else
            indx(1) =  1 + 1
            indx(4) = EF_grid_size - 1
        end if
        indx(2) = int(EF_grid_size/2.d0) + 1
        indx(3) = int(EF_grid_size * 2.d0/3.d0)
   end subroutine EF_index_for_current

  !> Get the electric field at a position in the domain (useful for the particle model)
  real(dp) function EF_get_at_pos(pos)
    real(dp), intent(in) :: pos
    real(dp) :: Efield_pos, temp
    integer :: lowIx

    ! EF_values(1) is defined at -0.5 * EF_delta_x
    lowIx = nint(pos * EF_inv_delta_x) + 1
    lowIx = min(EF_grid_size, max(1, lowIx))

    Efield_pos = (lowIx - 1.5_dp) * EF_delta_x
    temp = (pos - Efield_pos) * EF_inv_delta_x

    ! Do linear interpolation between lowIx and lowIx + 1 in the Efield array, given the position
    EF_get_at_pos = (1.0_dp - temp) * EF_values(lowIx) + temp * EF_values(lowIx+1)

  end function EF_get_at_pos

  !> Get a copy of the electric field at cell centers
  subroutine EF_get_values(out_efield)
    real(dp), intent(out) :: out_efield(:)
    out_efield(:) = 0.5_dp * (EF_values(1:EF_grid_size) + EF_values(2:EF_grid_size+1))
  end subroutine EF_get_values

  !> Get a copy of the electric field at cell faces (interior ones)
  subroutine EF_get_values_st(out_efield)
    real(dp), intent(out) :: out_efield(:)
    out_efield(:) = EF_values(2:EF_grid_size) ! Return only the interior points
  end subroutine EF_get_values_st

  real(dp) function EF_get_potential_at_bound () ! anbang: give the potential value of the left bound
      EF_get_potential_at_bound = EP_values(1)
  end function EF_get_potential_at_bound
	
  subroutine EF_get_Gap_and_dielectric_potential(EF_gap_potential , EF_dielectric_potential)
	real(dp), intent(out) :: EF_gap_potential , EF_dielectric_potential
	if (EF_useDBD) then 
		select case(EF_DBD_type)
		case(1)
			EF_dielectric_potential=EP_values(1)-EP_values(EF_DBD_index(1))
			EF_gap_potential=EP_values(EF_DBD_index(1))-EP_values(EF_grid_size)
		case(2)
			EF_dielectric_potential=EP_values(EF_DBD_index(2))-EP_values(EF_grid_size)
			EF_gap_potential=EP_values(1)-EP_values(EF_DBD_index(2))
		case(3)
			EF_gap_potential=EP_values(EF_DBD_index(1))-EP_values(EF_DBD_index(2))
			EF_dielectric_potential = EP_values(1) - EP_values(EF_grid_size)-EF_gap_potential
		end select
	end if	
  end subroutine EF_get_Gap_and_dielectric_potential

    !-------------------------------------------------------------------------------
    ! Subroutine for solving the Poisson equation without dielectric
    !-------------------------------------------------------------------------------
  SUBROUTINE PSol_without_DBD(net_charge,time,potential, E)
        use m_units_constants

        IMPLICIT NONE
    
        REAL(dp),DIMENSION(:),INTENT(in)    :: net_charge
        real(dp),intent(in)                 :: time 
        REAL(dp),DIMENSION(:),INTENT(OUT)   :: E, potential
        REAL(dp),DIMENSION(:),ALLOCATABLE   :: V,Qb,Qd,Qa,Qc
        INTEGER                             :: i
        real(dp)                            :: boundPot(2)
        
        allocate(Qa(EF_grid_size),Qb(EF_grid_size),Qc(EF_grid_size),Qd(EF_grid_size),V(EF_grid_size))

        Qa = 0.d0; Qb = 0.d0; Qc = 0.d0; Qd = 0.d0; V = 0.d0
        
    ! coefficient matrix, Anbang: Here we use uniform grids  
        do i= 2, EF_grid_size-1
        
            Qa(i) = -1.d0 / EF_delta_x
            Qb(i) = 2.d0 / EF_delta_x
            Qc(i) = -1.d0 / EF_delta_x
            Qd(i) = EF_delta_x * net_charge(i) / UC_eps0  !> Check net charge is right or not !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do

         ! set boundary conditions
         call EF_get_boundary_potential(boundPot, time)
         Qd(1) = boundPot(1)
         Qd(EF_grid_size) = boundPot(2)
        

         Qb(1) = 1.d0
         Qb(EF_grid_size) = 1.d0

        ! solve AES
        V = SolveLESTridiag(Qa,Qb,Qc,Qd)

        ! computation of axial electric field E = -grad(V)
        potential = V
        do i= 2,EF_grid_size
            E(i) = (V(i-1) - V(i)) *  EF_inv_delta_x
        end do
        ! Anbang: Here maybe needs further checking, how to set up boundary for efield
        E(1) = E(2)
        E(EF_grid_size+1) = E(EF_grid_size)
      
        deallocate(Qa,Qb,Qc,Qd,V)

   END SUBROUTINE PSol_without_DBD  

    !-------------------------------------------------------------------------------
    ! Subroutine for solving the Poisson equation for dielectric covered electrodes
    !-------------------------------------------------------------------------------
    !> Anbang: here is a little bit different from markus's code
    !> Here i include the dielectrics as well in the domain,
    !> so, we have to be careful when we calculate the potential
    SUBROUTINE PSol_DBD(net_charge, time, surChar, potential, E)
        use m_units_constants

        IMPLICIT NONE

        REAL(dp),DIMENSION(:),INTENT(in)    :: net_charge
        real(dp),intent(in)                 :: time 
        REAL(dp),DIMENSION(:),INTENT(OUT)   :: E, potential
        REAL(dp),DIMENSION(:),ALLOCATABLE   :: V,Qb,Qd,Qa,Qc
        INTEGER                             :: i
        real(dp),dimension(:),intent(in)    :: surChar    ! surface charge at two dielectrics
        real(dp)                            :: boundPot(2) ! the boundaty potential on the electrodes
        real(dp),allocatable                :: Vbulk(:)  ! the potential in the bulk plasma region
        integer                             :: varSize      
        
        ! here we find the size of bulk plasma region
        varSize = EF_DBD_index(2) - EF_DBD_index(1) + 1
        
        allocate(Qa(varSize),Qb(varSize),Qc(varSize),Qd(varSize),Vbulk(varSize),V(EF_grid_size))

        Qa = 0.d0; Qb = 0.d0; Qc = 0.d0; Qd = 0.d0; V = 0.d0
        
        ! coefficient matrix, Anbang: Here we use uniform grids  
        do i= 2, varSize-1
        
            Qa(i) = -1.d0 / EF_delta_x
            Qb(i) = 2.d0 / EF_delta_x
            Qc(i) = -1.d0 / EF_delta_x
            Qd(i) = EF_delta_x * net_charge(i + EF_DBD_index(1) - 1) / UC_eps0   !attention: the net charge has to shift to right 
        end do      

         ! set boundary conditions at the left electrode: powered electrode, right: ground electrode
         call EF_get_boundary_potential(boundPot, time)
        
        ! boundary conditions for the electric potential
        select case(EF_DBD_type)

            case(1) ! Dielectric at x = 0 only

                Qb(1) = UC_eps0 * EF_inv_delta_x + EF_epsDBD(1)/EF_DBD_len(1)
                Qc(1) = - UC_eps0 * EF_inv_delta_x
                Qd(1) = surChar(1) + EF_epsDBD(1) / EF_DBD_len(1) * boundPot(1)+net_charge(EF_DBD_index(1))*EF_delta_x/2.d0!
                
                Qb(varSize) = 1.d0
                Qd(varSize) = boundPot(2)

            case(2) ! Dielectric at x = d only

                Qb(1) = 1.d0
                Qd(1) = boundPot(1)
                
                Qa(varSize) = -UC_eps0 * EF_inv_delta_x
                Qb(varSize) = UC_eps0 / EF_delta_x + EF_epsDBD(2)/EF_DBD_len(2)
                Qd(varSize) = surChar(2) + EF_epsDBD(2) / EF_DBD_len(2) * boundPot(2)+net_charge(EF_DBD_index(2))*EF_delta_x/2.d0!

            case(3) ! Dielectric at x = 0 and x = d

                Qb(1) = UC_eps0 / EF_delta_x + EF_epsDBD(1)/ EF_DBD_len(1)
                Qc(1) = -UC_eps0 / EF_delta_x
                Qd(1) = surChar(1)+ EF_epsDBD(1)/ EF_DBD_len(1) * boundPot(1)+net_charge(EF_DBD_index(1))*EF_delta_x/2.d0!
                
                Qa(varSize) = -UC_eps0 / EF_delta_x
                Qb(varSize) = UC_eps0/ EF_delta_x + EF_epsDBD(2)/ EF_DBD_len(2)
                Qd(varSize) = surChar(2) + EF_epsDBD(2) / EF_DBD_len(2) * boundPot(2)+net_charge(EF_DBD_index(2))*EF_delta_x/2.d0!
        end select

    ! solve AES, here we get the potential in the bulk plasma region
        Vbulk = SolveLESTridiag(Qa,Qb,Qc,Qd)
       
     ! Anbang: Here we include the dielectric region
     select case(EF_DBD_type)
        
            case(1)
                V(1) = boundPot(1)
                V(EF_grid_size) = boundPot(2)
                do i = 2, EF_DBD_index(1) - 1
                    V(i) = V(1) + (Vbulk(1) - V(1)) / (EF_DBD_index(1) - 1) * (i - 1)
                end do

                do i = EF_DBD_index(1), EF_grid_size - 1
                    V(i) = Vbulk(i - EF_DBD_index(1) + 1)
                end do
                
            case(2)
                V(1) = boundPot(1)
                V(EF_grid_size) = boundPot(2)

                do i = 2, EF_DBD_index(2) 
                    V(i) = Vbulk(i)
                end do

                do i = EF_DBD_index(2) + 1, EF_grid_size - 1
                    V(i) = V(EF_DBD_index(2)) + (V(EF_grid_size) - V(EF_DBD_index(2))) / &
                            & (EF_grid_size - EF_DBD_index(2)) * (i - EF_DBD_index(2))
                end do

            case(3)
                 V(1) = boundPot(1)
                 V(EF_grid_size) = boundPot(2)

                 do i = 2, EF_DBD_index(1) - 1
                    V(i) = V(1) + (Vbulk(1) - V(1)) / (EF_DBD_index(1) - 1) * (i - 1)
                 end do  

                do i = EF_DBD_index(1), EF_DBD_index(2)
                     V(i) = Vbulk(i - EF_DBD_index(1) + 1)
                end do

                do i = EF_DBD_index(2) + 1, EF_grid_size - 1
                    V(i) = V(EF_DBD_index(2)) + (V(EF_grid_size) - V(EF_DBD_index(2))) / &
                            & (EF_grid_size - EF_DBD_index(2)) * (i - EF_DBD_index(2))
                end do

        end select
                

        ! computation of axial electric field E = -grad(V)
        potential = V
        do i= 2,EF_grid_size
            E(i) = (V(i-1) - V(i)) *  EF_inv_delta_x
        end do
    
        !Anbang: Here we need to pay attention to the electric field at the dielectric surface
        select case (EF_DBD_type)

            case(1)
                E(EF_DBD_index(1)) = E(EF_DBD_index(1) + 1)
            case(2)
                E(EF_DBD_index(2) + 1) = E(EF_DBD_index(2))
            case(3)
                E(EF_DBD_index(1)) = E(EF_DBD_index(1) + 1)
                E(EF_DBD_index(2) + 1) = E(EF_DBD_index(2))
        end select

        ! Anbang: Here maybe needs further checking, how to set up boundary for efield
        E(1) = E(2)
        E(EF_grid_size+1) = E(EF_grid_size)
        !print *, "here we finish the poisson solver!", E
        deallocate(Qa,Qb,Qc,Qd,V,Vbulk)
        
    END SUBROUTINE PSol_DBD    

   !solve tridiagonal matrix A * U = F, A is the matrix
   Function SolveLESTridiag(A,B,C,D) Result(X)
        !****f* ModMathESSolver/SolveLESTridiag
        ! NAME
        ! SolveLESTridiag -- solves tridiagonal systems of linear equations

        ! DESCRIPTION
        ! Solves linear equations in tridiagonal form using Gaussian
        ! elimination. The method is also known as Thomas algorithm.
        
        ! NOTES
        ! The main diagonal 'B' may not contain zero elements.

        ! USES
        ! * ModPrecision
        ! * ModSystem

        ! SYNOPSIS
        ! X = SolveLESTridiag(A,B,C,D)

        ! INPUTS
        ! * A   rank 1 array of reals: lower secondary diagonal
        ! * B   rank 1 array of reals: main diagonal
        ! * C   rank 1 array of reals: upper secondary diagonal
        ! * D   rank 1 array of reals: right hand side
        
        ! OUTPUT
        ! * X   solution vector of the linear equation system
        
        ! AUTHOR
        ! M.M. Becker

        ! SOURCE
        
        !USE ModSystem
        IMPLICIT NONE

        REAL(dp),INTENT(in),DIMENSION(:)  :: A,B,C,D
        REAL(dp),DIMENSION(size(D))       :: X
        REAL(dp),DIMENSION(:),ALLOCATABLE :: C_loc,D_loc
        REAL(dp)                          :: tmp
        INTEGER                           :: N,j

        N=size(D)

        allocate(C_loc(N),D_loc(N))
        
        if(any(B == 0.d0)) then
            print *, 'ERROR: In SolveLESTridiag:  &
            & the main diagonal may not contain zeros when not using pivoting.'
            stop
        end if

        C_loc(1) = C(1)/B(1)
        D_loc(1) = D(1)/B(1)
        do j=2,N
            tmp = 1.d0/(B(j)-A(j)*C_loc(j-1))
            C_loc(j) = tmp*C(j)
            D_loc(j) = tmp*(D(j)-A(j)*D_loc(j-1))
        end do

        X(N) = D_loc(N)
        do j=N-1,1,-1
            X(j) = D_loc(j)-C_loc(j)*X(j+1)
        end do

        deallocate(C_loc,D_loc)
        !****
    END Function SolveLESTridiag

   !> Anbang: boundary conditions for potential at both sides
   !> Left: powered electrde, right: ground electrode
    subroutine EF_get_boundary_potential(bd_pot, time)
         use m_units_constants

         real(dp), intent(out) :: bd_pot(2)
         real(dp), intent(in)  :: time
		 real(dp)              :: sim_time
		 integer               :: n
		 
         ! set boundary conditions at the left electrode: powered electrode
         select case (voltage_type)
         case (voltage_DC)
              bd_pot(1) = EL_getvoltage(time) 
         case (voltage_AC)
              bd_pot(1) = voltage_AC_max * sin(2.d0 * UC_pi * voltage_AC_fre * time)
		 case (voltage_pulse)
			  n=int(time * voltage_pulse_fre)
			  sim_time = time
			  sim_time = sim_time - real(n,dp) / voltage_pulse_fre
			  if(sim_time < voltage_pulse_rising_edge ) then
				bd_pot(1) = sim_time * voltage_pulse_max / voltage_pulse_rising_edge
			  else if(sim_time < (voltage_pulse_rising_edge + voltage_pulse_platform)) then
				bd_pot(1) = voltage_pulse_max
			  else if(sim_time < (voltage_pulse_rising_edge + voltage_pulse_platform + voltage_pulse_falling_edge)) then
				bd_pot(1) = voltage_pulse_max - (sim_time - (voltage_pulse_rising_edge + voltage_pulse_platform))*&
				(voltage_pulse_max / voltage_pulse_falling_edge)
			  else 
			  bd_pot(1) = 0.d0
			  end if 
         end select 
         
         ! right electrode: grounded electrode
         bd_pot(2) = 0.d0
    end subroutine EF_get_boundary_potential

    !anbang: get the index of positions of DBDtype
    integer function get_index_pos(pos)
        real(dp), intent(in)    :: pos
        
        get_index_pos = floor(pos * EF_inv_delta_x) + 1
        get_index_pos = min(EF_grid_size, max(1, get_index_pos))
    end function get_index_pos

   !Anbang: a routine to read the type and index
    subroutine EF_get_DBD_index_and_type(idx, ty)
        integer, intent(out) :: idx(2)
        integer, intent(out) :: ty

        idx = EF_DBD_index
        ty  = EF_DBD_type
    end subroutine EF_get_DBD_index_and_type



end module m_efield_1d
