#include "timH.h"
module tim_helperF

  use MOM_string_functions, only : uppercase
  use MOM_error_infra, only : MOM_err, FATAL
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use iso_c_binding, only : c_int
  use array_mod, only : RealArray_t
  use box_mod, only : Box_t

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

   public :: already_recorded, mark_recorded

   integer, parameter :: max_kernels = 1000
   character(len=128), save :: recorded(max_kernels)
   integer, save :: n_recorded = 0
   integer, parameter :: type_read = 1, type_write = 2

 type :: io_entry
  character(len=:), allocatable :: name
  character(len=:), allocatable :: type_name
  integer(kind=int64) :: offset
end type io_entry

type :: io_recorder
  integer :: unit_bin      ! binary file
  integer :: unit_meta     ! metadata file
  type(io_entry), allocatable :: entries(:)
  integer :: n = 0
  integer :: type 
contains
  procedure :: open_write
  procedure :: open_read
  procedure :: add_entry
  procedure :: close

  ! Write side 
  procedure :: add_realarray
  procedure :: add_box
  procedure :: add_real
  procedure :: add_integer
  procedure :: add_logical

  ! Read side
  procedure :: get_realarray
  procedure :: get_box
  procedure :: get_real
  procedure :: get_integer
  procedure :: get_logical

  procedure :: load_metadata
  procedure :: find_entry
  generic   :: add => add_realarray, add_box, add_real, add_integer, add_logical
  generic   :: get => get_realarray, get_box, get_real, get_integer, get_logical
end type io_recorder

contains

  function getenv_mode(name, default) result(mode)
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: default
    integer :: mode

    character(len=:), allocatable :: str
    integer :: length, status
    character(len=256) :: mesg

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
      write(mesg,'("tim_helperF::getenv_mod called with a ",A, &
               & " allowed values: AMREX, FORTRAN, CAPTURE")') TRIM(str)
      call MOM_err(FATAL,mesg)
    end select

  end function getenv_mode


  logical function already_recorded(name)
    character(len=*), intent(in) :: name
    integer :: i

    already_recorded = .false.
    do i = 1, n_recorded
      if (trim(recorded(i)) == trim(name)) then
        already_recorded = .true.
        return
      end if
    end do
  end function already_recorded

  subroutine mark_recorded(name)
    character(len=*), intent(in) :: name

    if (.not. already_recorded(name)) then
      if (n_recorded < max_kernels) then
        n_recorded = n_recorded + 1
        recorded(n_recorded) = name
      else
        print *, "Recorder registry full!"
      end if
    end if
  end subroutine mark_recorded

subroutine add_real(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  real(kind=real64), intent(in) :: val

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'real64', pos)

  write(this%unit_bin) val
end subroutine

subroutine get_real(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  real(kind=real64), intent(out) :: val

  integer :: idx
  integer(kind=int64) :: pos

  idx = this%find_entry(name)
  if (idx < 0) stop "get_real: not found"

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos) val
end subroutine

subroutine add_integer(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  integer, intent(in) :: val

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'integer', pos)

  write(this%unit_bin) val
end subroutine

subroutine get_integer(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  integer, intent(out) :: val

  integer :: idx
  integer(kind=int64) :: pos

  idx = this%find_entry(name)
  if (idx < 0) stop "get_integer: not found"

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos) val
end subroutine

subroutine add_logical(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  logical, intent(in) :: val

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'logical', pos)

  write(this%unit_bin) val
end subroutine

subroutine get_logical(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  logical, intent(out) :: val

  integer :: idx
  integer(kind=int64) :: pos

  idx = this%find_entry(name)
  if (idx < 0) stop "get_logical: not found"

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos) val
end subroutine

subroutine open_write(this, binfile, metafile)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: binfile, metafile

  open(newunit=this%unit_bin, file=binfile, &
       access='stream', form='unformatted', status='replace')

  open(newunit=this%unit_meta, file=metafile, &
       form='formatted', status='replace')

  this%type = type_write
  this%n = 0
end subroutine open_write

subroutine open_read(this, binfile, metafile)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: binfile, metafile

  open(newunit=this%unit_bin, file=binfile, &
       access='stream', form='unformatted', status='old')

  open(newunit=this%unit_meta, file=metafile, &
       form='formatted', status='old')

  call this%load_metadata()

  this%type = type_read

end subroutine open_read

subroutine load_metadata(this)
  class(io_recorder), intent(inout) :: this

  character(len=128) :: name, type_name
  integer(kind=int64) :: offset

  this%n = 0
  if (allocated(this%entries)) deallocate(this%entries)

  do
    read(this%unit_meta, *, end=100) name, type_name, offset

    this%n = this%n + 1
    if (.not. allocated(this%entries)) then
      allocate(this%entries(1))
    else
      this%entries = [this%entries, io_entry(name, type_name, offset)]
      cycle
    endif

    this%entries(1)%name = name
    this%entries(1)%type_name = type_name
    this%entries(1)%offset = offset
  enddo

  100 continue
end subroutine load_metadata

function find_entry(this, name) result(idx)
  class(io_recorder), intent(in) :: this
  character(*), intent(in) :: name
  integer :: idx
  integer :: i

  idx = -1
  do i = 1, this%n
    if (trim(this%entries(i)%name) == trim(name)) then
      idx = i
      return
    endif
  enddo
end function find_entry

subroutine add_entry(this, name, type_name, offset)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name, type_name
  integer(kind=int64), intent(in) :: offset

  this%n = this%n + 1

  if (.not. allocated(this%entries)) then
    allocate(this%entries(1))
  else
    this%entries = [this%entries, io_entry(name, type_name, offset)]
    return
  endif

  this%entries(1)%name = name
  this%entries(1)%type_name = type_name
  this%entries(1)%offset = offset
end subroutine add_entry

subroutine add_realarray(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  type(RealArray_t), intent(in) :: val

  integer(kind=int64) :: pos

  ! --- Get current file position ---
  inquire(unit=this%unit_bin, pos=pos)

  ! --- Register metadata ---
  call this%add_entry(name, 'RealArray_t', pos)

  ! --- Write binary payload ---
  call val%write_binary(this%unit_bin)

end subroutine add_realarray

subroutine get_realarray(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  type(RealArray_t), intent(inout) :: val

  integer :: idx
  integer(kind=int64) :: pos

  ! --- Find metadata entry ---
  idx = this%find_entry(name)
  if (idx < 0) then
    stop "get_realarray: variable not found: "//trim(name)
  endif

  ! --- Optional type check ---
  if (trim(this%entries(idx)%type_name) /= 'RealArray_t') then
    stop "get_realarray: type mismatch"
  endif

  pos = this%entries(idx)%offset

  ! --- Seek to position ---
  read(this%unit_bin, pos=pos)

  ! --- Delegate to type ---
  call val%read_binary(this%unit_bin)

end subroutine get_realarray

subroutine add_box(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  type(Box_T), intent(in) :: val

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'Box_t', pos)
  call val%write_binary(this%unit_bin)
end subroutine add_box

subroutine get_box(this, name, val)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: name
  type(Box_t), intent(inout) :: val

  integer :: idx
  integer(kind=int64) :: pos

  idx = this%find_entry(name)
  if (idx < 0) then
    stop "get_box: variable not found: "//trim(name)
  endif

  if (trim(this%entries(idx)%type_name) /= 'Box_t') then
    stop "get_box: type mismatch"
  endif

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos)

  call val%read_binary(this%unit_bin)

end subroutine

subroutine close(this)
  class(io_recorder), intent(inout) :: this
  integer :: i

  if(this%type.eq.type_write) then 
    do i = 1, this%n
      write(this%unit_meta,'(A,1X,A,1X,I0)') &
        this%entries(i)%name, &
        this%entries(i)%type_name, &
        this%entries(i)%offset
    enddo
  endif

  close(this%unit_bin)
  close(this%unit_meta)
end subroutine close

end module tim_helperF
