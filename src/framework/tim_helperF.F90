#include "timH.h"
module tim_helperF

  use MOM_string_functions, only : uppercase
  use iso_c_binding, only : c_int

  implicit none

  logical, parameter :: use_AMREX = .TRUE.
  integer, parameter :: TIMH_runAMREX   = TIMH_RUNAMREX_, &
                        TIMH_capture    = TIMH_CAPTURE_, &
                        TIMH_runFORTRAN = TIMH_RUNFORTRAN_
  integer(c_int), parameter :: TIMH_CAPTURE_INPUT = TIMH_CAPTURE_INPUT_,  &
                        TIMH_CAPTURE_OUTPUT = TIMH_CAPTURE_OUTPUT_, &
                        TIMH_RUN = TIMH_RUN_

   public :: TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN
   public :: TIMH_CAPTURE_INPUT, TIMH_CAPTURE_OUTPUT, TIMH_RUN
   public :: getenv_mode

contains

  function getenv_mode(name, default) result(mode)
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: default
    integer :: mode

    character(len=:), allocatable :: str
    integer :: length, status

    ! Get length first
    call get_environment_variable(name, length=length, status=status)

    if (status /= 0 .or. length == 0) then
      if (present(default)) then
        mode = default
      else
        mode = TIMH_runFORTRAN  ! sensible fallback
      end if
      return
    end if

    allocate(character(len=length) :: str)
    call get_environment_variable(name, str)

    ! Normalize (optional but recommended)
    str= uppercase(str)

    select case (trim(str))
    case ("AMREX")
      mode = TIMH_runAMREX
    case ("FORTRAN")
      mode = TIMH_runFORTRAN
    case ("CAPTURE")
      mode = TIMH_capture
    case default
      print *, "ERROR: Invalid value for ", name, " = ", trim(str)
      print *, "Allowed values: AMREX, FORTRAN, CAPTURE"
      stop 1
    end select

  end function getenv_mode

end module tim_helperF
