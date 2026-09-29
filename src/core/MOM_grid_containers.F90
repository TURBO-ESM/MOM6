! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3.2

!> Persistent array containers holding copies of the most heavily used ocean grid fields,
!! built once at initialization beside an unchanged ocean_grid_type.
module MOM_grid_containers

use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid, only : ocean_grid_type
use array_mod, only : RealArray_t, RealArray_C

implicit none ; private

public grid_core_type, grid_core_C
public grid_core_init, grid_core_end, grid_core

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A C-interoperable mirror of grid_core_type, with its members in the same order
type, bind(C) :: grid_core_C
  type(RealArray_C) :: mask2dT       !< 0 for land and 1 for ocean points on the h-grid [nondim]
  type(RealArray_C) :: areaT         !< The area of an h-cell [L2 ~> m2]
  type(RealArray_C) :: IareaT        !< 1/areaT [L-2 ~> m-2]
  type(RealArray_C) :: dxT           !< dxT is delta x at h points [L ~> m]
  type(RealArray_C) :: dyT           !< dyT is delta y at h points [L ~> m]
  type(RealArray_C) :: IdxT          !< 1/dxT [L-1 ~> m-1]
  type(RealArray_C) :: IdyT          !< IdyT is 1/dyT [L-1 ~> m-1]
  type(RealArray_C) :: bathyT        !< Bottom depth at h points, positive below Z_ref [Z ~> m]
  type(RealArray_C) :: mask2dCu      !< 0 for boundary and 1 for ocean points on the u grid [nondim]
  type(RealArray_C) :: dy_Cu         !< The unblocked lengths of the u-faces of the h-cell [L ~> m]
  type(RealArray_C) :: IdxCu         !< 1/dxCu [L-1 ~> m-1]
  type(RealArray_C) :: dxCu          !< dxCu is delta x at u points [L ~> m]
  type(RealArray_C) :: mask2dCv      !< 0 for boundary and 1 for ocean points on the v grid [nondim]
  type(RealArray_C) :: dx_Cv         !< The unblocked lengths of the v-faces of the h-cell [L ~> m]
  type(RealArray_C) :: IdyCv         !< 1/dyCv [L-1 ~> m-1]
  type(RealArray_C) :: dyCv          !< dyCv is delta y at v points [L ~> m]
  type(RealArray_C) :: CoriolisBu    !< The Coriolis parameter at corner points [T-1 ~> s-1]
  type(RealArray_C) :: meanSL        !< Mean sea level at h points, positive above Z_ref [Z ~> m]
  type(RealArray_C) :: geoLatT       !< Geographic latitude at h points [degrees_N] or [km] or [m]
  type(RealArray_C) :: IdyCu         !< 1/dyCu [L-1 ~> m-1]
  type(RealArray_C) :: areaCu        !< The areas of the u-grid cells [L2 ~> m2]
  type(RealArray_C) :: IdxCu_OBCmask !< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1]
  type(RealArray_C) :: IdxCv         !< 1/dxCv [L-1 ~> m-1]
  type(RealArray_C) :: areaCv        !< The areas of the v-grid cells [L2 ~> m2]
  type(RealArray_C) :: IdyCv_OBCmask !< 1/dyCv or 0 at boundary or OBC points [L-1 ~> m-1]
  type(RealArray_C) :: mask2dBu      !< 0 for boundary and 1 for ocean points on the q grid [nondim]
  type(RealArray_C) :: dxBu          !< dxBu is delta x at q points [L ~> m]
  type(RealArray_C) :: dyBu          !< dyBu is delta y at q points [L ~> m]
  type(RealArray_C) :: IareaBu       !< IareaBu = 1/areaBu [L-2 ~> m-2]
end type grid_core_C

!> Copies of the ocean grid fields that are used by 2 or more of the call trees being converted
!! to array containers, with the same names as in ocean_grid_type
type :: grid_core_type
  type(RealArray_t) :: mask2dT       !< 0 for land and 1 for ocean points on the h-grid [nondim]
  type(RealArray_t) :: areaT         !< The area of an h-cell [L2 ~> m2]
  type(RealArray_t) :: IareaT        !< 1/areaT [L-2 ~> m-2]
  type(RealArray_t) :: dxT           !< dxT is delta x at h points [L ~> m]
  type(RealArray_t) :: dyT           !< dyT is delta y at h points [L ~> m]
  type(RealArray_t) :: IdxT          !< 1/dxT [L-1 ~> m-1]
  type(RealArray_t) :: IdyT          !< IdyT is 1/dyT [L-1 ~> m-1]
  type(RealArray_t) :: bathyT        !< Bottom depth at h points, positive below Z_ref [Z ~> m]
  type(RealArray_t) :: mask2dCu      !< 0 for boundary and 1 for ocean points on the u grid [nondim]
  type(RealArray_t) :: dy_Cu         !< The unblocked lengths of the u-faces of the h-cell [L ~> m]
  type(RealArray_t) :: IdxCu         !< 1/dxCu [L-1 ~> m-1]
  type(RealArray_t) :: dxCu          !< dxCu is delta x at u points [L ~> m]
  type(RealArray_t) :: mask2dCv      !< 0 for boundary and 1 for ocean points on the v grid [nondim]
  type(RealArray_t) :: dx_Cv         !< The unblocked lengths of the v-faces of the h-cell [L ~> m]
  type(RealArray_t) :: IdyCv         !< 1/dyCv [L-1 ~> m-1]
  type(RealArray_t) :: dyCv          !< dyCv is delta y at v points [L ~> m]
  type(RealArray_t) :: CoriolisBu    !< The Coriolis parameter at corner points [T-1 ~> s-1]
  type(RealArray_t) :: meanSL        !< Mean sea level at h points, positive above Z_ref [Z ~> m]
  type(RealArray_t) :: geoLatT       !< Geographic latitude at h points [degrees_N] or [km] or [m]
  type(RealArray_t) :: IdyCu         !< 1/dyCu [L-1 ~> m-1]
  type(RealArray_t) :: areaCu        !< The areas of the u-grid cells [L2 ~> m2]
  type(RealArray_t) :: IdxCu_OBCmask !< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1]
  type(RealArray_t) :: IdxCv         !< 1/dxCv [L-1 ~> m-1]
  type(RealArray_t) :: areaCv        !< The areas of the v-grid cells [L2 ~> m2]
  type(RealArray_t) :: IdyCv_OBCmask !< 1/dyCv or 0 at boundary or OBC points [L-1 ~> m-1]
  type(RealArray_t) :: mask2dBu      !< 0 for boundary and 1 for ocean points on the q grid [nondim]
  type(RealArray_t) :: dxBu          !< dxBu is delta x at q points [L ~> m]
  type(RealArray_t) :: dyBu          !< dyBu is delta y at q points [L ~> m]
  type(RealArray_t) :: IareaBu       !< IareaBu = 1/areaBu [L-2 ~> m-2]
  type(grid_core_C) :: C             !< The C-interoperable mirror of this structure, set once
end type grid_core_type

!> The one persistent copy of the grid core fields
type(grid_core_type), target, save :: the_grid_core

contains

!> Build the persistent grid core containers from the final ocean grid.
subroutine grid_core_init(G)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure

  if (the_grid_core%mask2dT%associated()) &
    call MOM_error(FATAL, "grid_core_init: the grid core containers are already built.")

  call the_grid_core%mask2dT%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call the_grid_core%areaT%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call the_grid_core%IareaT%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call the_grid_core%dxT%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call the_grid_core%dyT%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call the_grid_core%IdxT%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call the_grid_core%IdyT%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call the_grid_core%bathyT%alloc(lb=LBOUND(G%bathyT), ub=UBOUND(G%bathyT), source=G%bathyT)
  call the_grid_core%mask2dCu%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call the_grid_core%dy_Cu%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call the_grid_core%IdxCu%alloc(lb=LBOUND(G%IdxCu), ub=UBOUND(G%IdxCu), source=G%IdxCu)
  call the_grid_core%dxCu%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call the_grid_core%mask2dCv%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call the_grid_core%dx_Cv%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call the_grid_core%IdyCv%alloc(lb=LBOUND(G%IdyCv), ub=UBOUND(G%IdyCv), source=G%IdyCv)
  call the_grid_core%dyCv%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)
  call the_grid_core%CoriolisBu%alloc(lb=LBOUND(G%CoriolisBu), ub=UBOUND(G%CoriolisBu), &
                                      source=G%CoriolisBu)
  call the_grid_core%meanSL%alloc(lb=LBOUND(G%meanSL), ub=UBOUND(G%meanSL), source=G%meanSL)
  call the_grid_core%geoLatT%alloc(lb=LBOUND(G%geoLatT), ub=UBOUND(G%geoLatT), source=G%geoLatT)
  call the_grid_core%IdyCu%alloc(lb=LBOUND(G%IdyCu), ub=UBOUND(G%IdyCu), source=G%IdyCu)
  call the_grid_core%areaCu%alloc(lb=LBOUND(G%areaCu), ub=UBOUND(G%areaCu), source=G%areaCu)
  call the_grid_core%IdxCu_OBCmask%alloc(lb=LBOUND(G%IdxCu_OBCmask), ub=UBOUND(G%IdxCu_OBCmask), &
                                         source=G%IdxCu_OBCmask)
  call the_grid_core%IdxCv%alloc(lb=LBOUND(G%IdxCv), ub=UBOUND(G%IdxCv), source=G%IdxCv)
  call the_grid_core%areaCv%alloc(lb=LBOUND(G%areaCv), ub=UBOUND(G%areaCv), source=G%areaCv)
  call the_grid_core%IdyCv_OBCmask%alloc(lb=LBOUND(G%IdyCv_OBCmask), ub=UBOUND(G%IdyCv_OBCmask), &
                                         source=G%IdyCv_OBCmask)
  call the_grid_core%mask2dBu%alloc(lb=LBOUND(G%mask2dBu), ub=UBOUND(G%mask2dBu), source=G%mask2dBu)
  call the_grid_core%dxBu%alloc(lb=LBOUND(G%dxBu), ub=UBOUND(G%dxBu), source=G%dxBu)
  call the_grid_core%dyBu%alloc(lb=LBOUND(G%dyBu), ub=UBOUND(G%dyBu), source=G%dyBu)
  call the_grid_core%IareaBu%alloc(lb=LBOUND(G%IareaBu), ub=UBOUND(G%IareaBu), source=G%IareaBu)

#ifdef _TIM
  the_grid_core%C = grid_core_to_c(the_grid_core)
#endif

end subroutine grid_core_init

!> Return a pointer to the persistent grid core containers, checking that they were built from
!! a grid with the same memory extents as G.
function grid_core(G) result(Gcore)
  type(ocean_grid_type), intent(in) :: G     !< The ocean's grid structure
  type(grid_core_type),  pointer    :: Gcore !< The persistent grid core containers

  if (.not. the_grid_core%mask2dT%associated()) &
    call MOM_error(FATAL, "grid_core: grid_core_init has not been called.")
  if (any(the_grid_core%mask2dT%lb(:) /= LBOUND(G%mask2dT)) .or. &
      any(the_grid_core%mask2dT%ub(:) /= UBOUND(G%mask2dT))) &
    call MOM_error(FATAL, "grid_core: the grid core containers do not match the extents of G.")

  Gcore => the_grid_core

end function grid_core

#ifdef _TIM
!> Convert the grid core containers to their C-interoperable mirror.
function grid_core_to_c(Gcore) result(Gcore_c)
  type(grid_core_type), target, intent(in) :: Gcore   !< The grid core containers
  type(grid_core_C)                        :: Gcore_c !< The C-interoperable mirror of Gcore

  Gcore_c%mask2dT       = Gcore%mask2dT%to_c()
  Gcore_c%areaT         = Gcore%areaT%to_c()
  Gcore_c%IareaT        = Gcore%IareaT%to_c()
  Gcore_c%dxT           = Gcore%dxT%to_c()
  Gcore_c%dyT           = Gcore%dyT%to_c()
  Gcore_c%IdxT          = Gcore%IdxT%to_c()
  Gcore_c%IdyT          = Gcore%IdyT%to_c()
  Gcore_c%bathyT        = Gcore%bathyT%to_c()
  Gcore_c%mask2dCu      = Gcore%mask2dCu%to_c()
  Gcore_c%dy_Cu         = Gcore%dy_Cu%to_c()
  Gcore_c%IdxCu         = Gcore%IdxCu%to_c()
  Gcore_c%dxCu          = Gcore%dxCu%to_c()
  Gcore_c%mask2dCv      = Gcore%mask2dCv%to_c()
  Gcore_c%dx_Cv         = Gcore%dx_Cv%to_c()
  Gcore_c%IdyCv         = Gcore%IdyCv%to_c()
  Gcore_c%dyCv          = Gcore%dyCv%to_c()
  Gcore_c%CoriolisBu    = Gcore%CoriolisBu%to_c()
  Gcore_c%meanSL        = Gcore%meanSL%to_c()
  Gcore_c%geoLatT       = Gcore%geoLatT%to_c()
  Gcore_c%IdyCu         = Gcore%IdyCu%to_c()
  Gcore_c%areaCu        = Gcore%areaCu%to_c()
  Gcore_c%IdxCu_OBCmask = Gcore%IdxCu_OBCmask%to_c()
  Gcore_c%IdxCv         = Gcore%IdxCv%to_c()
  Gcore_c%areaCv        = Gcore%areaCv%to_c()
  Gcore_c%IdyCv_OBCmask = Gcore%IdyCv_OBCmask%to_c()
  Gcore_c%mask2dBu      = Gcore%mask2dBu%to_c()
  Gcore_c%dxBu          = Gcore%dxBu%to_c()
  Gcore_c%dyBu          = Gcore%dyBu%to_c()
  Gcore_c%IareaBu       = Gcore%IareaBu%to_c()

end function grid_core_to_c
#endif

!> Free the persistent grid core containers.
subroutine grid_core_end()

  call the_grid_core%mask2dT%free() ; call the_grid_core%areaT%free()
  call the_grid_core%IareaT%free() ; call the_grid_core%dxT%free()
  call the_grid_core%dyT%free() ; call the_grid_core%IdxT%free()
  call the_grid_core%IdyT%free() ; call the_grid_core%bathyT%free()
  call the_grid_core%mask2dCu%free() ; call the_grid_core%dy_Cu%free()
  call the_grid_core%IdxCu%free() ; call the_grid_core%dxCu%free()
  call the_grid_core%mask2dCv%free() ; call the_grid_core%dx_Cv%free()
  call the_grid_core%IdyCv%free() ; call the_grid_core%dyCv%free()
  call the_grid_core%CoriolisBu%free() ; call the_grid_core%meanSL%free()
  call the_grid_core%geoLatT%free() ; call the_grid_core%IdyCu%free()
  call the_grid_core%areaCu%free() ; call the_grid_core%IdxCu_OBCmask%free()
  call the_grid_core%IdxCv%free() ; call the_grid_core%areaCv%free()
  call the_grid_core%IdyCv_OBCmask%free() ; call the_grid_core%mask2dBu%free()
  call the_grid_core%dxBu%free() ; call the_grid_core%dyBu%free()
  call the_grid_core%IareaBu%free()

end subroutine grid_core_end

end module MOM_grid_containers
