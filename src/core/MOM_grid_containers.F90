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
public grid_OBC_type, grid_OBC_C
public grid_OBC_init, grid_OBC_end, grid_OBC

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A C-interoperable mirror of grid_core_type, with its members in the same order
type, bind(C) :: grid_core_C
  type(RealArray_C) :: mask2dT       !< 0 for land points and 1 for ocean points on the h-grid
                                     !! [nondim].
  type(RealArray_C) :: areaT         !< The area of an h-cell [L2 ~> m2].
  type(RealArray_C) :: IareaT        !< 1/areaT [L-2 ~> m-2].
  type(RealArray_C) :: dxT           !< dxT is delta x at h points [L ~> m].
  type(RealArray_C) :: dyT           !< dyT is delta y at h points [L ~> m].
  type(RealArray_C) :: IdxT          !< 1/dxT [L-1 ~> m-1].
  type(RealArray_C) :: IdyT          !< IdyT is 1/dyT [L-1 ~> m-1].
  type(RealArray_C) :: bathyT        !< Ocean bottom depth, referenced to Z_ref at tracer points.
                                     !! bathyT is in depth units and positive *below* Z_ref
                                     !! [Z ~> m].
  type(RealArray_C) :: mask2dCu      !< 0 for boundary points and 1 for ocean points on the u grid
                                     !! [nondim].
  type(RealArray_C) :: dy_Cu         !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
  type(RealArray_C) :: IdxCu         !< 1/dxCu [L-1 ~> m-1].
  type(RealArray_C) :: dxCu          !< dxCu is delta x at u points [L ~> m].
  type(RealArray_C) :: mask2dCv      !< 0 for boundary points and 1 for ocean points on the v grid
                                     !! [nondim].
  type(RealArray_C) :: dx_Cv         !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
  type(RealArray_C) :: IdyCv         !< 1/dyCv [L-1 ~> m-1].
  type(RealArray_C) :: dyCv          !< dyCv is delta y at v points [L ~> m].
  type(RealArray_C) :: CoriolisBu    !< The Coriolis parameter at corner points [T-1 ~> s-1].
  type(RealArray_C) :: meanSL        !< Spatially varying time mean sea level, referenced to Z_ref
                                     !! at tracer points. meanSL is in height units and positive
                                     !! *above* Z_ref. It is used a) as the height where p = p_atm
                                     !! or zero; b) to calculate time mean thickness of the water
                                     !! column, where mean thickness = max(meanSL + bathyT, 0.0).
                                     !! meanSL is 2D for the consideration of a domain with
                                     !! spatically varying mean height, e.g. the Great Lakes system
                                     !! [Z ~> m].
  type(RealArray_C) :: geoLatT       !< The geographic latitude at tracer (h) points
                                     !! [degrees_N] or [km] or [m]
  type(RealArray_C) :: IdyCu         !< 1/dyCu [L-1 ~> m-1].
  type(RealArray_C) :: areaCu        !< The areas of the u-grid cells [L2 ~> m2].
  type(RealArray_C) :: IdxCv         !< 1/dxCv [L-1 ~> m-1].
  type(RealArray_C) :: areaCv        !< The areas of the v-grid cells [L2 ~> m2].
  type(RealArray_C) :: mask2dBu      !< 0 for boundary points and 1 for ocean points on the q grid
                                     !! [nondim].
  type(RealArray_C) :: dxBu          !< dxBu is delta x at q points [L ~> m].
  type(RealArray_C) :: dyBu          !< dyBu is delta y at q points [L ~> m].
  type(RealArray_C) :: IareaBu       !< IareaBu = 1/areaBu [L-2 ~> m-2].
  type(RealArray_C) :: geoLonT       !< The geographic longitude at tracer (h) points
                                     !! [degrees_E] or [km] or [m]
  type(RealArray_C) :: df_dx         !< Derivative d/dx f (Coriolis parameter) at h-points
                                     !! [T-1 L-1 ~> s-1 m-1].
  type(RealArray_C) :: df_dy         !< Derivative d/dy f (Coriolis parameter) at h-points
                                     !! [T-1 L-1 ~> s-1 m-1].
  type(RealArray_C) :: IareaCu       !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
  type(RealArray_C) :: IareaCv       !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
  type(RealArray_C) :: IdxBu         !< 1/dxBu [L-1 ~> m-1].
  type(RealArray_C) :: IdyBu         !< 1/dyBu [L-1 ~> m-1].
  type(RealArray_C) :: Coriolis2Bu   !< The square of the Coriolis parameter at corner points
                                     !! [T-2 ~> s-2].
end type grid_core_C

!> Copies of the ocean grid fields that are used by 2 or more of the call trees being converted
!! to array containers, with the same names as in ocean_grid_type
type :: grid_core_type
  type(RealArray_t) :: mask2dT       !< 0 for land points and 1 for ocean points on the h-grid
                                     !! [nondim].
  type(RealArray_t) :: areaT         !< The area of an h-cell [L2 ~> m2].
  type(RealArray_t) :: IareaT        !< 1/areaT [L-2 ~> m-2].
  type(RealArray_t) :: dxT           !< dxT is delta x at h points [L ~> m].
  type(RealArray_t) :: dyT           !< dyT is delta y at h points [L ~> m].
  type(RealArray_t) :: IdxT          !< 1/dxT [L-1 ~> m-1].
  type(RealArray_t) :: IdyT          !< IdyT is 1/dyT [L-1 ~> m-1].
  type(RealArray_t) :: bathyT        !< Ocean bottom depth, referenced to Z_ref at tracer points.
                                     !! bathyT is in depth units and positive *below* Z_ref
                                     !! [Z ~> m].
  type(RealArray_t) :: mask2dCu      !< 0 for boundary points and 1 for ocean points on the u grid
                                     !! [nondim].
  type(RealArray_t) :: dy_Cu         !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
  type(RealArray_t) :: IdxCu         !< 1/dxCu [L-1 ~> m-1].
  type(RealArray_t) :: dxCu          !< dxCu is delta x at u points [L ~> m].
  type(RealArray_t) :: mask2dCv      !< 0 for boundary points and 1 for ocean points on the v grid
                                     !! [nondim].
  type(RealArray_t) :: dx_Cv         !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
  type(RealArray_t) :: IdyCv         !< 1/dyCv [L-1 ~> m-1].
  type(RealArray_t) :: dyCv          !< dyCv is delta y at v points [L ~> m].
  type(RealArray_t) :: CoriolisBu    !< The Coriolis parameter at corner points [T-1 ~> s-1].
  type(RealArray_t) :: meanSL        !< Spatially varying time mean sea level, referenced to Z_ref
                                     !! at tracer points. meanSL is in height units and positive
                                     !! *above* Z_ref. It is used a) as the height where p = p_atm
                                     !! or zero; b) to calculate time mean thickness of the water
                                     !! column, where mean thickness = max(meanSL + bathyT, 0.0).
                                     !! meanSL is 2D for the consideration of a domain with
                                     !! spatically varying mean height, e.g. the Great Lakes system
                                     !! [Z ~> m].
  type(RealArray_t) :: geoLatT       !< The geographic latitude at tracer (h) points
                                     !! [degrees_N] or [km] or [m]
  type(RealArray_t) :: IdyCu         !< 1/dyCu [L-1 ~> m-1].
  type(RealArray_t) :: areaCu        !< The areas of the u-grid cells [L2 ~> m2].
  type(RealArray_t) :: IdxCv         !< 1/dxCv [L-1 ~> m-1].
  type(RealArray_t) :: areaCv        !< The areas of the v-grid cells [L2 ~> m2].
  type(RealArray_t) :: mask2dBu      !< 0 for boundary points and 1 for ocean points on the q grid
                                     !! [nondim].
  type(RealArray_t) :: dxBu          !< dxBu is delta x at q points [L ~> m].
  type(RealArray_t) :: dyBu          !< dyBu is delta y at q points [L ~> m].
  type(RealArray_t) :: IareaBu       !< IareaBu = 1/areaBu [L-2 ~> m-2].
  type(RealArray_t) :: geoLonT       !< The geographic longitude at tracer (h) points
                                     !! [degrees_E] or [km] or [m]
  type(RealArray_t) :: df_dx         !< Derivative d/dx f (Coriolis parameter) at h-points
                                     !! [T-1 L-1 ~> s-1 m-1].
  type(RealArray_t) :: df_dy         !< Derivative d/dy f (Coriolis parameter) at h-points
                                     !! [T-1 L-1 ~> s-1 m-1].
  type(RealArray_t) :: IareaCu       !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
  type(RealArray_t) :: IareaCv       !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
  type(RealArray_t) :: IdxBu         !< 1/dxBu [L-1 ~> m-1].
  type(RealArray_t) :: IdyBu         !< 1/dyBu [L-1 ~> m-1].
  type(RealArray_t) :: Coriolis2Bu   !< The square of the Coriolis parameter at corner points
                                     !! [T-2 ~> s-2].
  type(grid_core_C) :: C             !< The C-interoperable mirror of this structure, set once
end type grid_core_type

!> A C-interoperable mirror of grid_OBC_type, with its members in the same order
type, bind(C) :: grid_OBC_C
  type(RealArray_C) :: OBCmaskCu     !< 0 for boundary or OBC points and 1 for ocean points on the u
                                     !! grid [nondim].
  type(RealArray_C) :: IdxCu_OBCmask !< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
  type(RealArray_C) :: OBCmaskCv     !< 0 for boundary or OBC points and 1 for ocean points on the v
                                     !! grid [nondim].
  type(RealArray_C) :: IdyCv_OBCmask !< 1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
end type grid_OBC_C

!> Copies of the ocean grid fields that mask boundary and open boundary condition points, with
!! the same names as in ocean_grid_type
type :: grid_OBC_type
  type(RealArray_t) :: OBCmaskCu     !< 0 for boundary or OBC points and 1 for ocean points on the u
                                     !! grid [nondim].
  type(RealArray_t) :: IdxCu_OBCmask !< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
  type(RealArray_t) :: OBCmaskCv     !< 0 for boundary or OBC points and 1 for ocean points on the v
                                     !! grid [nondim].
  type(RealArray_t) :: IdyCv_OBCmask !< 1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
  type(grid_OBC_C)  :: C             !< The C-interoperable mirror of this structure, set once
end type grid_OBC_type

!> The one persistent copy of the grid core fields
type(grid_core_type), target, save :: the_grid_core

!> The one persistent copy of the grid OBC mask fields
type(grid_OBC_type), target, save :: the_grid_OBC

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
  call the_grid_core%IdxCv%alloc(lb=LBOUND(G%IdxCv), ub=UBOUND(G%IdxCv), source=G%IdxCv)
  call the_grid_core%areaCv%alloc(lb=LBOUND(G%areaCv), ub=UBOUND(G%areaCv), source=G%areaCv)
  call the_grid_core%mask2dBu%alloc(lb=LBOUND(G%mask2dBu), ub=UBOUND(G%mask2dBu), source=G%mask2dBu)
  call the_grid_core%dxBu%alloc(lb=LBOUND(G%dxBu), ub=UBOUND(G%dxBu), source=G%dxBu)
  call the_grid_core%dyBu%alloc(lb=LBOUND(G%dyBu), ub=UBOUND(G%dyBu), source=G%dyBu)
  call the_grid_core%IareaBu%alloc(lb=LBOUND(G%IareaBu), ub=UBOUND(G%IareaBu), source=G%IareaBu)
  call the_grid_core%geoLonT%alloc(lb=LBOUND(G%geoLonT), ub=UBOUND(G%geoLonT), source=G%geoLonT)
  call the_grid_core%df_dx%alloc(lb=LBOUND(G%df_dx), ub=UBOUND(G%df_dx), source=G%df_dx)
  call the_grid_core%df_dy%alloc(lb=LBOUND(G%df_dy), ub=UBOUND(G%df_dy), source=G%df_dy)
  call the_grid_core%IareaCu%alloc(lb=LBOUND(G%IareaCu), ub=UBOUND(G%IareaCu), source=G%IareaCu)
  call the_grid_core%IareaCv%alloc(lb=LBOUND(G%IareaCv), ub=UBOUND(G%IareaCv), source=G%IareaCv)
  call the_grid_core%IdxBu%alloc(lb=LBOUND(G%IdxBu), ub=UBOUND(G%IdxBu), source=G%IdxBu)
  call the_grid_core%IdyBu%alloc(lb=LBOUND(G%IdyBu), ub=UBOUND(G%IdyBu), source=G%IdyBu)
  call the_grid_core%Coriolis2Bu%alloc(lb=LBOUND(G%Coriolis2Bu), ub=UBOUND(G%Coriolis2Bu), &
                                       source=G%Coriolis2Bu)

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
  Gcore_c%IdxCv         = Gcore%IdxCv%to_c()
  Gcore_c%areaCv        = Gcore%areaCv%to_c()
  Gcore_c%mask2dBu      = Gcore%mask2dBu%to_c()
  Gcore_c%dxBu          = Gcore%dxBu%to_c()
  Gcore_c%dyBu          = Gcore%dyBu%to_c()
  Gcore_c%IareaBu       = Gcore%IareaBu%to_c()
  Gcore_c%geoLonT       = Gcore%geoLonT%to_c()
  Gcore_c%df_dx         = Gcore%df_dx%to_c()
  Gcore_c%df_dy         = Gcore%df_dy%to_c()
  Gcore_c%IareaCu       = Gcore%IareaCu%to_c()
  Gcore_c%IareaCv       = Gcore%IareaCv%to_c()
  Gcore_c%IdxBu         = Gcore%IdxBu%to_c()
  Gcore_c%IdyBu         = Gcore%IdyBu%to_c()
  Gcore_c%Coriolis2Bu   = Gcore%Coriolis2Bu%to_c()

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
  call the_grid_core%areaCu%free()
  call the_grid_core%IdxCv%free() ; call the_grid_core%areaCv%free()
  call the_grid_core%mask2dBu%free()
  call the_grid_core%dxBu%free() ; call the_grid_core%dyBu%free()
  call the_grid_core%IareaBu%free() ; call the_grid_core%geoLonT%free()
  call the_grid_core%df_dx%free() ; call the_grid_core%df_dy%free()
  call the_grid_core%IareaCu%free() ; call the_grid_core%IareaCv%free()
  call the_grid_core%IdxBu%free() ; call the_grid_core%IdyBu%free()
  call the_grid_core%Coriolis2Bu%free()

end subroutine grid_core_end

!> Build the persistent grid OBC mask containers from the final ocean grid.
subroutine grid_OBC_init(G)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure

  if (the_grid_OBC%OBCmaskCu%associated()) &
    call MOM_error(FATAL, "grid_OBC_init: the grid OBC containers are already built.")

  call the_grid_OBC%OBCmaskCu%alloc(lb=LBOUND(G%OBCmaskCu), ub=UBOUND(G%OBCmaskCu), &
                                    source=G%OBCmaskCu)
  call the_grid_OBC%IdxCu_OBCmask%alloc(lb=LBOUND(G%IdxCu_OBCmask), ub=UBOUND(G%IdxCu_OBCmask), &
                                        source=G%IdxCu_OBCmask)
  call the_grid_OBC%OBCmaskCv%alloc(lb=LBOUND(G%OBCmaskCv), ub=UBOUND(G%OBCmaskCv), &
                                    source=G%OBCmaskCv)
  call the_grid_OBC%IdyCv_OBCmask%alloc(lb=LBOUND(G%IdyCv_OBCmask), ub=UBOUND(G%IdyCv_OBCmask), &
                                        source=G%IdyCv_OBCmask)

#ifdef _TIM
  the_grid_OBC%C = grid_OBC_to_c(the_grid_OBC)
#endif

end subroutine grid_OBC_init

!> Return a pointer to the persistent grid OBC mask containers, checking that they were built
!! from a grid with the same memory extents as G.
function grid_OBC(G) result(Gobc)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure
  type(grid_OBC_type),   pointer    :: Gobc !< The persistent grid OBC mask containers

  if (.not. the_grid_OBC%OBCmaskCu%associated()) &
    call MOM_error(FATAL, "grid_OBC: grid_OBC_init has not been called.")
  if (any(the_grid_OBC%OBCmaskCu%lb(:) /= LBOUND(G%OBCmaskCu)) .or. &
      any(the_grid_OBC%OBCmaskCu%ub(:) /= UBOUND(G%OBCmaskCu))) &
    call MOM_error(FATAL, "grid_OBC: the grid OBC containers do not match the extents of G.")

  Gobc => the_grid_OBC

end function grid_OBC

#ifdef _TIM
!> Convert the grid OBC mask containers to their C-interoperable mirror.
function grid_OBC_to_c(Gobc) result(Gobc_c)
  type(grid_OBC_type), target, intent(in) :: Gobc   !< The grid OBC mask containers
  type(grid_OBC_C)                        :: Gobc_c !< The C-interoperable mirror of Gobc

  Gobc_c%OBCmaskCu     = Gobc%OBCmaskCu%to_c()
  Gobc_c%IdxCu_OBCmask = Gobc%IdxCu_OBCmask%to_c()
  Gobc_c%OBCmaskCv     = Gobc%OBCmaskCv%to_c()
  Gobc_c%IdyCv_OBCmask = Gobc%IdyCv_OBCmask%to_c()

end function grid_OBC_to_c
#endif

!> Free the persistent grid OBC mask containers.
subroutine grid_OBC_end()

  call the_grid_OBC%OBCmaskCu%free() ; call the_grid_OBC%IdxCu_OBCmask%free()
  call the_grid_OBC%OBCmaskCv%free() ; call the_grid_OBC%IdyCv_OBCmask%free()

end subroutine grid_OBC_end

end module MOM_grid_containers
