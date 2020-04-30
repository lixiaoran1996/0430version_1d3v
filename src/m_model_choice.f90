module m_model_choice
  use m_types
  use m_error

  implicit none
  private

  ! Enumerate the different simulation types
  integer, parameter :: MODEL_part = 1
  integer, parameter :: MODEL_fluid = 2

  integer :: MODEL_type = -1

  ! Public types
  public :: MODEL_part, MODEL_fluid

  ! Public routines
  public :: MODEL_initialize
  public :: MODEL_get_type

contains

  subroutine MODEL_initialize(sim_type_name)
    character(len=*), intent(in) :: sim_type_name

    select case (sim_type_name)
    case ("particle")
       MODEL_type = MODEL_part
    case ("fluid")
       MODEL_type = MODEL_fluid
    case default
       call ERR_show("MODEL_initialize: Invalid simulation type given: " // sim_type_name)
    end select

  end subroutine MODEL_initialize

  integer function MODEL_get_type()
    MODEL_get_type = MODEL_type
    if (MODEL_get_type == -1) call ERR_show("MODEL_get_type: call MODEL_initialize first")
  end function MODEL_get_type

end module m_model_choice
