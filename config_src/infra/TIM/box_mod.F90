module box_mod
  use iso_c_binding, only : c_ptr, c_loc
  implicit none
  private

  public :: Box_t, Box_c

  !< Box struct for C bridge
  type, bind(C) :: Box_C
     type(c_ptr) :: idxS   !< start index of a box
     type(c_ptr) :: idxE   !< end index of a box
  end type Box_C

  !< Class box_t defines an iteration index range
  type :: Box_T
     integer, pointer :: idxS(:) => NULL()  !< Start index of a box
     integer, pointer :: idxE(:) => NULL()  !< End index of a box
  contains
     procedure   :: allocBox           !< allocate index box
     procedure   :: setBox             !< Sets the index range for the bocx
     procedure   :: freeBox            !< deallocates index box
     procedure   :: to_c_Box           !< Converts Box to C
     procedure   :: expandBox          !< Increase the bounds of a box in one dimension
                                       !! both extents of box are increased by a fixed amount
     procedure   :: contractBox        !< Decrease the bounds of a box in one dimension
                                       !! both extents of box are decreased by a fixed amount 
     generic     :: alloc => allocBox  !< Allocate memory for a box
     generic     :: free => freeBox    !< Deallocates memory used by box
     generic     :: set => setBox      !< Set the extent of a box
     generic     :: to_c => to_c_Box   !< Convert box to C
     generic     :: expand => expandBox !< Expand the extent of a box
     generic     :: contract => contractBox !< Contract the exxtent of a box
  end type Box_T

contains

!< Allocates an iteration box
subroutine allocBox(this,ndims)
  class(Box_t), intent(inout) :: this   !< The box to be allocated
  integer, intent(in) :: ndims          !< The number of dimension in the box

  ! If already associated deallocate
  if(associated(this%idxS)) deallocate(this%idxS) 
  if(associated(this%idxE)) deallocate(this%idxS)

  allocate(this%idxS(ndims), source=0)
  allocate(this%idxE(ndims), source=0)
end subroutine allocBox

!< Allocates an iteration box
subroutine freeBox(this)
  class(Box_t), intent(inout) :: this   !< The box to be deallocated

  if(associated(this%idxS)) deallocate(this%idxS)
  if(associated(this%idxE)) deallocate(this%idxE)

end subroutine freeBox

!< Set the extents of the iteration box
subroutine setBox(this,idxS,idxE)
  class(Box_t), intent(inout) :: this        !< The box to set
  integer, dimension(:), intent(in) :: idxS  !< The starting indices
  integer, dimension(:), intent(in) :: idxE  !< The ending indices

  if(associated(this%idxS)) this%idxS(:)=idxS(:)
  if(associated(this%idxE)) this%idxE(:)=idxE(:)

end subroutine setBox

!< Return a new box with expanded iteration extents
function expandBox(this,dim,n) result(new)
  class(Box_t), intent(in) :: this !< The iteration box to modify
  integer, intent(in)      :: dim  !< The dimension to expand
  integer, intent(in)      :: n    !< The extent of the expansion
  type(Box_t) :: new

  ! Local variables
  integer ::rank

  if(associated(new%idxS)) deallocate(new%idxS)
  if(associated(new%idxE)) deallocate(new%idxE)

  rank = SIZE(this%idxS)
  allocate(new%idxS(rank),new%idxE(rank))
  new%idxS(:) = this%idxS(:)
  new%idxE(:) = this%idxE(:)
  new%idxS(dim) = new%idxS(dim)-n
  new%idxE(dim) = new%idxE(dim)+n

end function expandBox

!< Return a new box with contracted iteration extents
function contractBox(this,dim,n) result(new)
  class(Box_t), intent(in) :: this   !< The iteration box to modify
  integer, intent(in)      :: dim    !< The dimension to contract
  integer, intent(in)      :: n      !< The length of the contraction
  type(Box_t) :: new

  ! Local variables
  integer ::rank

  if(associated(new%idxS)) deallocate(new%idxS)
  if(associated(new%idxE)) deallocate(new%idxE)

  rank = SIZE(this%idxS)
  allocate(new%idxS(rank),new%idxE(rank))
  new%idxS(:) = this%idxS(:)
  new%idxE(:) = this%idxE(:)
  new%idxS(dim) = new%idxS(dim)+n
  new%idxE(dim) = new%idxE(dim)-n

end function contractBox

!< Convert Fortran box to C
function to_c_Box(this) result(cdesc)
  class(Box_t), intent(in) :: this  !< The box to convert
  type(Box_C) :: cdesc              !< C compatible pointers
  cdesc%idxS = c_loc(this%idxS)
  cdesc%idxE = c_loc(this%idxE)
end function to_c_Box

end module box_mod
