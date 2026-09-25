! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3

#include "do_concurrent_compat.h"

!> Accelerations due to the Coriolis force and momentum advection
module MOM_CoriolisAdv

!> \author Robert Hallberg, April 1994 - June 2002

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary, only : OBC_VORTICITY_ZERO, OBC_VORTICITY_FREESLIP
use MOM_open_boundary, only : OBC_VORTICITY_COMPUTED, OBC_VORTICITY_SPECIFIED
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : accel_diag_ptrs, porous_barrier_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
use box_mod,           only : Box_t

implicit none ; private

public CorAdCalc, CoriolisAdv_init, CoriolisAdv_end, CoriolisAdv_stencil

#include <MOM_memory.h>

!> Control structure for mom_coriolisadv
type, public :: CoriolisAdv_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  integer :: Coriolis_Scheme !< Selects the discretization for the Coriolis terms.
                             !! Valid values are:
                             !! - SADOURNY75_ENERGY - Sadourny, 1975
                             !! - ARAKAWA_HSU90     - Arakawa & Hsu, 1990, Energy & non-div. Enstrophy
                             !! - ROBUST_ENSTRO     - Pseudo-enstrophy scheme
                             !! - SADOURNY75_ENSTRO - Sadourny, JAS 1975, Enstrophy
                             !! - ARAKAWA_LAMB81    - Arakawa & Lamb, MWR 1981, Energy & Enstrophy
                             !! - ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu and Sadourny energy.
                             !! - WENOVI3RD_PV_ENSTRO    - 3rd-order WENO scheme for PV reconstruction
                             !! - WENOVI5TH_PV_ENSTRO    - 5th-order WENO scheme for PV reconstruction
                             !! - WENOVI7TH_PV_ENSTRO    - 7th-order WENO scheme for PV reconstruction
                             !! The default, SADOURNY75_ENERGY, is the safest choice then the
                             !! deformation radius is poorly resolved.
  integer :: KE_Scheme       !< KE_SCHEME selects the discretization for
                             !! the kinetic energy. Valid values are:
                             !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
  logical :: KE_use_limiter  !< If true, use the Koren limiter for KE_UP3 scheme
  integer :: PV_Adv_Scheme   !< PV_ADV_SCHEME selects the discretization for PV advection
                             !! Valid values are:
                             !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                             !! - PV_ADV_UPWIND1  - upwind, first order
  integer :: nkblock         !< The k block size used in Coriolis/advection calculations [nondim].
  real    :: F_eff_max_blend !< The factor by which the maximum effective Coriolis
                             !! acceleration from any point can be increased when
                             !! blending different discretizations with the
                             !! ARAKAWA_LAMB_BLEND Coriolis scheme [nondim].
                             !! This must be greater than 2.0, and is 4.0 by default.
  real    :: wt_lin_blend    !< A weighting value beyond which the blending between
                             !! Sadourny and Arakawa & Hsu goes linearly to 0 [nondim].
                             !! This must be between 1 and 1e-15, often 1/8.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Coriolis  !< If true, the Coriolis terms at u points are
                             !! bounded by the four estimates of (f+rv)v from the
                             !! four neighboring v points, and similarly at v
                             !! points.  This option would have no effect on the
                             !! SADOURNY75_ENERGY scheme if it were possible to
                             !! use centered difference thickness fluxes.
  logical :: Coriolis_En_Dis !< If CORIOLIS_EN_DIS is defined, two estimates of
                             !! the thickness fluxes are used to estimate the
                             !! Coriolis term, and the one that dissipates energy
                             !! relative to the other one is used.  This is only
                             !! available at present if Coriolis scheme is
                             !! SADOURNY75_ENERGY.
  logical :: weno_velocity_smooth !< If true, use velocity to compute the smoothness indicator for WENO
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_rv = -1, id_PV = -1, id_gKEu = -1, id_gKEv = -1
  integer :: id_rvxu = -1, id_rvxv = -1
  ! integer :: id_hf_gKEu    = -1, id_hf_gKEv    = -1
  integer :: id_hf_gKEu_2d = -1, id_hf_gKEv_2d = -1
  integer :: id_intz_gKEu_2d = -1, id_intz_gKEv_2d = -1
  ! integer :: id_hf_rvxu    = -1, id_hf_rvxv    = -1
  integer :: id_hf_rvxu_2d = -1, id_hf_rvxv_2d = -1
  integer :: id_h_gKEu = -1, id_h_gKEv = -1
  integer :: id_h_rvxu = -1, id_h_rvxv = -1
  integer :: id_intz_rvxu_2d = -1, id_intz_rvxv_2d = -1
  integer :: id_CAuS = -1, id_CAvS = -1
  !>@}
end type CoriolisAdv_CS

!>@{ Enumeration values for Coriolis_Scheme
integer, parameter :: SADOURNY75_ENERGY = 1
integer, parameter :: ARAKAWA_HSU90     = 2
integer, parameter :: ROBUST_ENSTRO     = 3
integer, parameter :: SADOURNY75_ENSTRO = 4
integer, parameter :: ARAKAWA_LAMB81    = 5
integer, parameter :: AL_BLEND          = 6
integer, parameter :: wenovi7th_PV_ENSTRO = 7
integer, parameter :: wenovi5th_PV_ENSTRO = 8
integer, parameter :: wenovi3rd_PV_ENSTRO = 9
character*(20), parameter :: SADOURNY75_ENERGY_STRING = "SADOURNY75_ENERGY"
character*(20), parameter :: ARAKAWA_HSU_STRING = "ARAKAWA_HSU90"
character*(20), parameter :: ROBUST_ENSTRO_STRING = "ROBUST_ENSTRO"
character*(20), parameter :: SADOURNY75_ENSTRO_STRING = "SADOURNY75_ENSTRO"
character*(20), parameter :: ARAKAWA_LAMB_STRING = "ARAKAWA_LAMB81"
character*(20), parameter :: AL_BLEND_STRING = "ARAKAWA_LAMB_BLEND"
character*(20), parameter :: WENOVI7TH_PV_ENSTRO_STRING = "WENOVI7TH_PV_ENSTRO"
character*(20), parameter :: WENOVI5TH_PV_ENSTRO_STRING = "WENOVI5TH_PV_ENSTRO"
character*(20), parameter :: WENOVI3RD_PV_ENSTRO_STRING = "WENOVI3RD_PV_ENSTRO"
!>@}
!>@{ Enumeration values for KE_Scheme
integer, parameter :: KE_ARAKAWA        = 10
integer, parameter :: KE_SIMPLE_GUDONOV = 11
integer, parameter :: KE_GUDONOV        = 12
integer, parameter :: KE_UP3            = 13
character*(20), parameter :: KE_ARAKAWA_STRING = "KE_ARAKAWA"
character*(20), parameter :: KE_SIMPLE_GUDONOV_STRING = "KE_SIMPLE_GUDONOV"
character*(20), parameter :: KE_GUDONOV_STRING = "KE_GUDONOV"
character*(20), parameter :: KE_UP3_STRING = "KE_UP3"
!>@}
!>@{ Enumeration values for PV_Adv_Scheme
integer, parameter :: PV_ADV_CENTERED   = 21
integer, parameter :: PV_ADV_UPWIND1    = 22
character*(20), parameter :: PV_ADV_CENTERED_STRING = "PV_ADV_CENTERED"
character*(20), parameter :: PV_ADV_UPWIND1_STRING = "PV_ADV_UPWIND1"
!>@}

contains

!> Calculates the Coriolis and momentum advection contributions to the acceleration.
subroutine CorAdCalc(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves)
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uh !< Zonal transport u*h*dy
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vh !< Meridional transport v*h*dx
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out)   :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out)   :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                      intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv
  type(porous_barrier_type),                  intent(in)    :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS

  call CorAdCalc_TR(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves)

end subroutine CorAdCalc

!> The implementation of CorAdCalc (the root of its call tree): sets up the k-invariant fields,
!! then loops over iteration tiles, calling CorAdv_tile for each one.
subroutine CorAdCalc_TR(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves)
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uh !< Zonal transport u*h*dy
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vh !< Meridional transport v*h*dx
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out)   :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out)   :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                      intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv
  type(porous_barrier_type),                  intent(in)    :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS

  ! Local variables
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    Area_q      ! The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Area_h      ! The ocean area at h points [L2 ~> m2].  Area_h is used to find the
                ! average thickness in the denominator of q.  0 for land points.
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    PV, &       ! A diagnostic array of the potential vorticities [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    RV          ! A diagnostic array of the relative vorticities [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: CAuS ! Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: CAvS ! Stokes contribution to CAv [L T-2 ~> m s-2]
  real :: vol_neglect            ! A volume so small that is expected to be
                                 ! lost in roundoff [H L2 ~> m3 or kg].
  real :: area_neglect           ! An area so small that is expected to be
                                 ! lost in roundoff [L2 ~> m2].
  real :: eps_vel                ! A tiny, positive velocity [L T-1 ~> m s-1].
  real :: Fe_m2         ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: rat_lin       ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: h_tiny        ! A very small thickness [H ~> m or kg m-2].
  type(Box_t) :: bxH    ! The h-point iteration box of one tile, [isc:iec, jsc:jec, k_start:k_end]
  type(Box_t) :: bxQ    ! The B-grid-index iteration box of one tile, [IscB:IecB, JscB:JecB, k_start:k_end]
  type(Box_t) :: bxQs   ! The scheme-dependent vorticity-point box of one tile,
                        ! [Is_q:Ie_q, Js_q:Je_q, k_start:k_end]
  integer :: i, j, n, is, ie, js, je, nz, nkblock, k_start, k_end
  integer :: Is_q, Ie_q, Js_q, Je_q  ! The scheme-dependent range of values at which vorticity is set.
  logical :: Stokes_VF
  logical :: use_weno   ! True if using one of the WENO schemes
  integer :: stencil    ! Stencil size of WENO scheme

! To work, the following fields must be set outside of the usual
! is to ie range before this subroutine is called:
!   v(is-1:ie+2,js-1:je+1), u(is-1:ie+1,js-1:je+2), h(is-1:ie+2,js-1:je+2),
!   uh(is-1,ie,js:je+1) and vh(is:ie+1,js-1:je).

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_CoriolisAdv: Module must be initialized before it is used.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkblock = merge(GV%ke, CS%nkblock, CS%nkblock==0)
  vol_neglect = GV%H_subroundoff * (1e-4 * US%m_to_L)**2
  area_neglect = (1e-4 * US%m_to_L)**2
  eps_vel = 1.0e-10*US%m_s_to_L_T
  h_tiny = GV%Angstrom_H  ! Perhaps this should be set to h_neglect instead.

  stencil = CoriolisAdv_stencil(CS)

  use_weno = CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO &
      .or. CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO &
      .or. CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO

  if (use_weno) then
    Is_q = is - stencil ; Ie_q = ie + stencil - 1 ; Js_q = js - stencil ; Je_q = je + stencil - 1
  else
    Is_q = G%IscB - 1 ; Ie_q = G%IecB + 1 ; Js_q = G%JscB - 1 ; Je_q = G%JecB + 1
  endif

  !$omp target enter data map(alloc: Area_h, Area_q)

  do concurrent (j=Js_q:Je_q+1, I=Is_q:Ie_q+1)
    Area_h(i,j) = G%mask2dT(i,j) * G%areaT(i,j)
  enddo

  if (associated(OBC)) then
    !$omp target update from(Area_h)

    do n=1,OBC%number_of_segments
      if (.not. OBC%segment(n)%on_pe) cycle
      I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
      if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
        do i = max(Is_q,OBC%segment(n)%HI%isd), min(Ie_q+1,OBC%segment(n)%HI%ied)
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            Area_h(i,j+1) = Area_h(i,j)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            Area_h(i,j) = Area_h(i,j+1)
          endif
        enddo
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
        do j = max(Js_q,OBC%segment(n)%HI%jsd), min(Je_q+1,OBC%segment(n)%HI%jed)
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            Area_h(i+1,j) = Area_h(i,j)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            Area_h(i,j) = Area_h(i+1,j)
          endif
        enddo
      endif
    enddo

    !$omp target update to(Area_h)
  endif

  do concurrent (J=Js_q:Je_q, I=Is_q:Ie_q)
    Area_q(i,j) = (Area_h(i,j) + Area_h(i+1,j+1)) + &
                  (Area_h(i+1,j) + Area_h(i,j+1))
  enddo

  Stokes_VF = .false.
  if (present(Waves)) then ; if (associated(Waves)) then
    Stokes_VF = Waves%Stokes_VF
  endif ; endif

  !$omp target enter data map(alloc: CAuS, CAvS) if (Stokes_VF)

  ! Diagnostics
  !$omp target enter data map(alloc: RV) if (CS%id_RV > 0)
  !$omp target enter data map(alloc: PV) if (CS%id_PV > 0)
  !$omp target enter data map(alloc: AD%gradKEu) if (associated(AD%gradKEu))
  !$omp target enter data map(alloc: AD%gradKEv) if (associated(AD%gradKEv))
  !$omp target enter data map(alloc: AD%rv_x_u) if (associated(AD%rv_x_u))
  !$omp target enter data map(alloc: AD%rv_x_v) if (associated(AD%rv_x_v))

  ! TODO: Do this outside of the function
  !$omp target enter data map(to: pbv, pbv%por_face_areaU, pbv%por_face_areaV) &
  !$omp   if (CS%Coriolis_En_Dis)

  ! Hoist AL_BLEND k-independent scalars out of the tile loop.
  Fe_m2 = 0.0 ; rat_lin = 0.0
  if (CS%Coriolis_Scheme == AL_BLEND) then
    Fe_m2 = CS%F_eff_max_blend - 2.0
    rat_lin = 1.5 * Fe_m2 / max(CS%wt_lin_blend, 1.0e-16)
    if (CS%F_eff_max_blend <= 2.0) then ; Fe_m2 = -1. ; rat_lin = -1.0 ; endif
  endif

  ! Loop over iteration tiles.  Each tile covers the whole horizontal compute domain and
  ! nkblock layers; the tile boxes carry all of the index ranges used by CorAdv_tile.
  call bxH%safe_alloc(ndims=3) ; call bxQ%safe_alloc(ndims=3) ; call bxQs%safe_alloc(ndims=3)
  do k_start=1,nz,nkblock
    k_end = min(k_start+nkblock-1, nz)
    call bxH%set(idxS=[G%isc,G%jsc,k_start], idxE=[G%iec,G%jec,k_end])
    call bxQ%set(idxS=[G%IscB,G%JscB,k_start], idxE=[G%IecB,G%JecB,k_end])
    call bxQs%set(idxS=[Is_q,Js_q,k_start], idxE=[Ie_q,Je_q,k_end])

    call CorAdv_tile(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves, &
                     bxH, bxQ, bxQs, Area_h, Area_q, Stokes_VF, use_weno, &
                     vol_neglect, area_neglect, eps_vel, h_tiny, Fe_m2, rat_lin, &
                     RV, PV, CAuS, CAvS)
  enddo ! end of tile loop.
  call bxH%free() ; call bxQ%free() ; call bxQs%free()

  !$omp target exit data map(delete: Area_h, Area_q)

  ! TODO: Move outside function
  !$omp target exit data map(delete: pbv, pbv%por_face_areaU, pbv%por_face_areaV) &
  !$omp   if (CS%Coriolis_En_Dis)

  ! Diagnostics
  !$omp target exit data map(from: RV) if (CS%id_RV > 0)
  !$omp target exit data map(from: PV) if (CS%id_PV > 0)
  !$omp target exit data map(from: AD%gradKEu) if (associated(AD%gradKEu))
  !$omp target exit data map(from: AD%gradKEv) if (associated(AD%gradKEv))
  !$omp target exit data map(from: AD%rv_x_u) if (associated(AD%rv_x_u))
  !$omp target exit data map(from: AD%rv_x_v) if (associated(AD%rv_x_v))
  !$omp target exit data map(from: CAuS, CAvS) if (Stokes_VF)

  call CorAdv_finalize_diagnostics(RV, PV, CAuS, CAvS, Stokes_VF, AD, G, GV, CS)
end subroutine CorAdCalc_TR


!> Calculates the Coriolis and momentum advection contributions to the acceleration over one
!! iteration tile, by calling the setup, the selected scheme and the common terms in turn.
subroutine CorAdv_tile(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves, &
                       bxH, bxQ, bxQs, Area_h, Area_q, Stokes_VF, use_weno, &
                       vol_neglect, area_neglect, eps_vel, h_tiny, Fe_m2, rat_lin, &
                       RV, PV, CAuS, CAvS)
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  type(Box_t),                                intent(in)    :: bxH  !< The h-point iteration box of this
                                                                    !! tile, [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),                                intent(in)    :: bxQ  !< The B-grid-index iteration box of
                                                                    !! this tile, [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),                                intent(in)    :: bxQs !< The scheme-dependent vorticity-point
                                                                    !! box, [Is_q:Ie_q, Js_q:Je_q, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uh !< Zonal transport u*h*dy
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vh !< Meridional transport v*h*dx
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                      intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv
  type(porous_barrier_type),                  intent(in)    :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: Area_h !< The ocean area at h points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)),         intent(in)    :: Area_q !< The sum of the ocean areas at the 4
                                                                  !! adjacent thickness points [L2 ~> m2]
  logical,                                    intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  logical,                                    intent(in)    :: use_weno  !< True if using one of the WENO schemes
  real,                                       intent(in)    :: vol_neglect  !< A negligible volume [H L2 ~> m3 or kg]
  real,                                       intent(in)    :: area_neglect !< A negligible area [L2 ~> m2]
  real,                                       intent(in)    :: eps_vel !< A tiny, positive velocity [L T-1 ~> m s-1]
  real,                                       intent(in)    :: h_tiny  !< A very small thickness [H ~> m or kg m-2]
  real,                                       intent(in)    :: Fe_m2   !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  real,                                       intent(in)    :: rat_lin !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), intent(inout) :: RV !< Diagnostic relative vorticity [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), intent(inout) :: PV !< Diagnostic potential vorticity
                                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),  intent(inout) :: CAuS !< Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(inout) :: CAvS !< Stokes contribution to CAv [L T-2 ~> m s-2]

  ! Local variables
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    q, &        ! Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    qS, &       ! Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    Ih_q, &     ! The inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    h_q, &      ! The thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    abs_vort, & ! Absolute vorticity at q-points [T-1 ~> s-1].
    q2          ! Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    KEx, &      ! The zonal gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEx = d/dx KE.
    uh_center   ! Transport based on arithmetic mean h at u-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    KE          ! Kinetic energy per unit mass [L2 T-2 ~> m2 s-2], KE = (u^2 + v^2)/2.
  real, dimension(SZI_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    KEy, &      ! The meridional gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEy = d/dy KE.
    vh_center   ! Transport based on arithmetic mean h at v-points [H L2 T-1 ~> m3 s-1 or kg s-1]

  !$omp target enter data map(alloc: abs_vort, q, Ih_q)
  !$omp target enter data map(alloc: h_q) if (use_weno)
  !$omp target enter data map(alloc: KE, KEx, KEy)
  ! TODO: These Stokes_VF fields seem associated with diagnostics
  !$omp target enter data map(alloc: qS) if (Stokes_VF)
  !$omp target enter data map(alloc: uh_center, vh_center) if (CS%Coriolis_En_Dis)
  !$omp target enter data map(alloc: q2) &
  !$omp   if(associated(AD%rv_x_u) .or. associated(AD%rv_x_v))

  ! Potential vorticity and the related quantities at q points used by every scheme.
  call CorAdv_setup(u, v, h, OBC, AD, pbv, Waves, G, GV, bxH, bxQ, bxQs, Area_h, Area_q, &
                    Stokes_VF, use_weno, vol_neglect, area_neglect, q, Ih_q, abs_vort, h_q, qS, q2, &
                    uh_center, vh_center, RV, PV, CS)

  ! Calculate KE and the gradient of KE
  call gradKE(u, v, h, KE, KEx, KEy, bxH, bxQ, G, GV, US, CS)
  ! TODO: Can KE be removed from this function?

  ! Calculate the Coriolis and momentum advection accelerations, CAu and CAv, with the selected
  ! scheme.  On a Cartesian grid, CAu =  q * vh - d(KE)/dx and CAv = - q * uh - d(KE)/dy.
  select case (CS%Coriolis_Scheme)
  case (SADOURNY75_ENERGY, SADOURNY75_ENSTRO)
    call CorAdv_sadourny(u, v, uh, vh, q, uh_center, vh_center, G, GV, bxH, bxQ, CAu, CAv, CS)
  case (ARAKAWA_HSU90, ARAKAWA_LAMB81, AL_BLEND)
    call CorAdv_arakawa(uh, vh, q, Ih_q, Fe_m2, rat_lin, G, GV, bxH, bxQ, CAu, CAv, CS)
  case (ROBUST_ENSTRO)
    call CorAdv_robust_enstro(u, v, h, uh, vh, abs_vort, eps_vel, h_tiny, G, GV, bxH, bxQ, CAu, CAv, CS)
  case (wenovi7th_PV_ENSTRO, wenovi5th_PV_ENSTRO, wenovi3rd_PV_ENSTRO)
    call CorAdv_weno(u, v, uh, vh, q, abs_vort, h_q, G, GV, bxH, bxQ, CAu, CAv, CS)
  end select

  ! The Stokes-drift diagnostic, bounding, the kinetic energy gradient, and diagnostics.
  call CorAdv_common_terms(u, v, uh, vh, abs_vort, qS, q2, KEx, KEy, Stokes_VF, G, GV, bxH, bxQ, &
                           CAu, CAv, CAuS, CAvS, AD, CS)

  !$omp target exit data map(delete: abs_vort, q, Ih_q)
  !$omp target exit data map(delete: h_q) if (use_weno)
  !$omp target exit data map(delete: KE, KEx, KEy)
  !$omp target exit data map(delete: qS) if (Stokes_VF)
  !$omp target exit data map(delete: uh_center, vh_center) if (CS%Coriolis_En_Dis)
  !$omp target exit data map(delete: q2) &
  !$omp     if(associated(AD%rv_x_u) .or. associated(AD%rv_x_v))
end subroutine CorAdv_tile


!> Calculates the potential vorticity and the related quantities at q points that every
!! Coriolis scheme uses, over one iteration tile.
subroutine CorAdv_setup(u, v, h, OBC, AD, pbv, Waves, G, GV, bxH, bxQ, bxQs, Area_h, Area_q, &
                        Stokes_VF, use_weno, vol_neglect, area_neglect, q, Ih_q, abs_vort, h_q, qS, &
                        q2, uh_center, vh_center, RV, PV, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxQs !< The scheme-dependent vorticity-point box,
                                                 !! [Is_q:Ie_q, Js_q:Je_q, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),   intent(in)    :: AD  !< Storage for acceleration diagnostics
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS), optional, pointer :: Waves !< An optional pointer to Stokes drift CS
  real, dimension(SZI_(G),SZJ_(G)),   intent(in) :: Area_h !< The ocean area at h points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: Area_q !< The sum of the ocean areas at the 4
                                                           !! adjacent thickness points [L2 ~> m2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  logical,                 intent(in)    :: use_weno  !< True if using one of the WENO schemes
  real,                    intent(in)    :: vol_neglect !< A negligible volume [H L2 ~> m3 or kg]
  real,                    intent(in)    :: area_neglect !< A negligible area [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: q !< Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: Ih_q !< Inverse of thickness interpolated to q points
                                                    !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: abs_vort !< Absolute vorticity at q-points [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: h_q !< Thickness interpolated to q points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: qS !< Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: q2 !< Relative vorticity over thickness
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: uh_center !< Transport based on arithmetic mean h at u-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(out) :: vh_center !< Transport based on arithmetic mean h at v-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), intent(inout) :: RV !< Diagnostic relative vorticity [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), intent(inout) :: PV !< Diagnostic potential vorticity
                                                                   !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    dvdx, dudy, &   ! Contributions to the circulation around q-points [L2 T-1 ~> m2 s-1]
    dvSdx, duSdy, & ! idem. for Stokes drift [L2 T-1 ~> m2 s-1]
    rel_vort, &     ! Relative vorticity at q-points [T-1 ~> s-1].
    stk_vort        ! Stokes vorticity at q-points [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    hArea_u         ! The cell area weighted thickness interpolated to u points
                    ! times the effective areas [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    hArea_v         ! The cell area weighted thickness interpolated to v points
                    ! times the effective areas [H L2 ~> m3 or kg].
  real :: hArea_q   ! The sum of area times thickness of the cells
                    ! surrounding a q point [H L2 ~> m3 or kg].
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, ksc, kec, n
  integer :: Is_q, Ie_q, Js_q, Je_q  ! The scheme-dependent range of values at which vorticity is set.

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)
  Is_q = bxQs%idxS(1) ; Ie_q = bxQs%idxE(1) ; Js_q = bxQs%idxS(2) ; Je_q = bxQs%idxE(2)

  !$omp target enter data map(alloc: dvdx, dudy)
  !$omp target enter data map(alloc: hArea_u, hArea_v)
  !$omp target enter data map(alloc: rel_vort)
  !$omp target enter data map(alloc: dvSdx, duSdy, stk_vort) if (Stokes_VF)

  ! Here the second order accurate layer potential vorticities, q,
  ! are calculated.  hq is  second order accurate in space.  Relative
  ! vorticity is second order accurate everywhere with free slip b.c.s,
  ! but only first order accurate at boundaries with no slip b.c.s.
  ! First calculate the contributions to the circulation around the q-point.
  if (Stokes_VF) then
    if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
      do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
        dvSdx(I,J,k) = (-Waves%us_y(i+1,J,k)*G%dyCv(i+1,J)) - &
                         (-Waves%us_y(i,J,k)*G%dyCv(i,J))
        duSdy(I,J,k) = (-Waves%us_x(I,j+1,k)*G%dxCu(I,j+1)) - &
                         (-Waves%us_x(I,j,k)*G%dxCu(I,j))
      enddo
    endif
    if (.not. Waves%Passive_Stokes_VF) then
      do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
        dvdx(I,J,k) = ((v(i+1,J,k)-Waves%us_y(i+1,J,k))*G%dyCv(i+1,J)) - &
                        ((v(i,J,k)-Waves%us_y(i,J,k))*G%dyCv(i,J))
        dudy(I,J,k) = ((u(I,j+1,k)-Waves%us_x(I,j+1,k))*G%dxCu(I,j+1)) - &
                        ((u(I,j,k)-Waves%us_x(I,j,k))*G%dxCu(I,j))
      enddo
    else
      do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
        dvdx(I,J,k) = (v(i+1,J,k)*G%dyCv(i+1,J)) - (v(i,J,k)*G%dyCv(i,J))
        dudy(I,J,k) = (u(I,j+1,k)*G%dxCu(I,j+1)) - (u(I,j,k)*G%dxCu(I,j))
      enddo
    endif
  else
    do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
      dvdx(I,J,k) = (v(i+1,J,k)*G%dyCv(i+1,J)) - (v(i,J,k)*G%dyCv(i,J))
      dudy(I,J,k) = (u(I,j+1,k)*G%dxCu(I,j+1)) - (u(I,j,k)*G%dxCu(I,j))
    enddo
  endif

  do concurrent (k=ksc:kec, J=Js_q:Je_q, i=Is_q:Ie_q+1)
    hArea_v(i,J,k) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i,j+1) * h(i,j+1,k)))
  enddo

  do concurrent (k=ksc:kec, j=Js_q:Je_q+1, I=Is_q:Ie_q)
    hArea_u(I,j,k) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i+1,j) * h(i+1,j,k)))
  enddo

  if (CS%Coriolis_En_Dis) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq+1, I=is-1:ie)
      uh_center(I,j,k) = 0.5 * ((G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k)) * (h(i,j,k) + h(i+1,j,k))
    enddo

    do concurrent (k=ksc:kec, J=js-1:je, i=Isq:Ieq+1)
      vh_center(i,J,k) = 0.5 * ((G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k)) * (h(i,j,k) + h(i,j+1,k))
    enddo
  endif

  ! Adjust circulation components to relative vorticity and thickness projected onto
  ! velocity points on open boundaries.
  if (associated(OBC)) then
    !$omp target update from(Area_h)
    do k=ksc,kec ! TODO: port (OBC GPU path not yet implemented)
      !$omp target update from(dvdx(:,:,k), dudy(:,:,k))
      !$omp target update from(hArea_u(:,:,k), hArea_v(:,:,k))
      !$omp target update from(uh_center(:,:,k), vh_center(:,:,k))

      do n=1,OBC%number_of_segments
        if (.not. OBC%segment(n)%on_pe) cycle

        I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB

        if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
          select case (OBC%vorticity_config)
          case (OBC_VORTICITY_ZERO)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              dvdx(I,J,k) = 0. ; dudy(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_FREESLIP)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              dudy(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_COMPUTED)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J,k) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - u(I,j,k))*G%dxCu(I,j)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                dudy(I,J,k) = 2.0*(u(I,j+1,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%dxCu(I,j+1)
              endif
            enddo
          case (OBC_VORTICITY_SPECIFIED)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*G%dxCu(I,j)*G%dyBu(I,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                dudy(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*G%dxCu(I,j+1)*G%dyBu(I,J)
              endif
            enddo
          end select

          ! Project thicknesses across OBC points with a no-gradient condition.
          do i=max(Is_q,OBC%segment(n)%HI%isd), min(Ie_q+1,OBC%segment(n)%HI%ied)
            if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
              hArea_v(i,J,k) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              hArea_v(i,J,k) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
            endif
          enddo

          if (CS%Coriolis_En_Dis) then
            do i=max(Isq,OBC%segment(n)%HI%isd), min(Ieq+1,OBC%segment(n)%HI%ied)
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                vh_center(i,J,k) = (G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j,k)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                vh_center(i,J,k) = (G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j+1,k)
              endif
            enddo
          endif
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
          select case (OBC%vorticity_config)
          case (OBC_VORTICITY_ZERO)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              dvdx(I,J,k) = 0. ; dudy(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_FREESLIP)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              dvdx(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_COMPUTED)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J,k) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - v(i,J,k))*G%dyCv(i,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                dvdx(I,J,k) = 2.0*(v(i+1,J,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%dyCv(i+1,J)
              endif
            enddo
          case (OBC_VORTICITY_SPECIFIED)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*G%dyCv(i,J)*G%dxBu(I,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                dvdx(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*G%dyCv(i+1,J)*G%dxBu(I,J)
              endif
            enddo
          end select

          ! Project thicknesses across OBC points with a no-gradient condition.
          do j=max(Js_q,OBC%segment(n)%HI%jsd), min(Je_q+1,OBC%segment(n)%HI%jed)
            if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
              hArea_u(I,j,k) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              hArea_u(I,j,k) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i+1,j,k)
            endif
          enddo

          if (CS%Coriolis_En_Dis) then
            do j=max(Jsq,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                uh_center(I,j,k) = (G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k) * h(i,j,k)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                uh_center(I,j,k) = (G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k) * h(i+1,j,k)
              endif
            enddo
          endif
        endif
      enddo

      ! Now project thicknesses across cell-corner points in the OBCs.  The two
      ! projections have to occur in sequence and can not be combined easily.
      do n=1,OBC%number_of_segments
        if (.not. OBC%segment(n)%on_pe) cycle
        I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
          do I = max(Is_q,OBC%segment(n)%HI%IsdB), min(Ie_q,OBC%segment(n)%HI%IedB)
            if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
              if (Area_h(i,j) + Area_h(i+1,j) > 0.0) then
                hArea_u(I,j+1,k) = hArea_u(I,j,k) * ((Area_h(i,j+1) + Area_h(i+1,j+1)) / &
                                                         (Area_h(i,j) + Area_h(i+1,j)))
              else ; hArea_u(I,j+1,k) = 0.0 ; endif
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              if (Area_h(i,j+1) + Area_h(i+1,j+1) > 0.0) then
                hArea_u(I,j,k) = hArea_u(I,j+1,k) * ((Area_h(i,j) + Area_h(i+1,j)) / &
                                                         (Area_h(i,j+1) + Area_h(i+1,j+1)))
              else ; hArea_u(I,j,k) = 0.0 ; endif
            endif
          enddo
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
          do J = max(Js_q,OBC%segment(n)%HI%JsdB), min(Je_q,OBC%segment(n)%HI%JedB)
            if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
              if (Area_h(i,j) + Area_h(i,j+1) > 0.0) then
                hArea_v(i+1,J,k) = hArea_v(i,J,k) * ((Area_h(i+1,j) + Area_h(i+1,j+1)) / &
                                                         (Area_h(i,j) + Area_h(i,j+1)))
              else ; hArea_v(i+1,J,k) = 0.0 ; endif
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              hArea_v(i,J,k) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
              if (Area_h(i+1,j) + Area_h(i+1,j+1) > 0.0) then
                hArea_v(i,J,k) = hArea_v(i+1,J,k) * ((Area_h(i,j) + Area_h(i,j+1)) / &
                                                         (Area_h(i+1,j) + Area_h(i+1,j+1)))
              else ; hArea_v(i,J,k) = 0.0 ; endif
            endif
          enddo
        endif
      enddo

      !$omp target update to(dvdx(:,:,k), dudy(:,:,k))
      !$omp target update to(hArea_u(:,:,k), hArea_v(:,:,k))
      !$omp target update to(uh_center(:,:,k), vh_center(:,:,k))
    enddo ! k=ksc,kec OBC
  endif

  if (CS%no_slip) then
    do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
      rel_vort(I,J,k) = (2.0 - G%mask2dBu(I,J)) * (dvdx(I,J,k) - dudy(I,J,k)) * G%IareaBu(I,J)
    enddo

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        do concurrent (k=ksc:kec, J=Jsq-1:Jeq+1, I=Isq-1:Ieq+1)
          stk_vort(I,J,k) = (2.0 - G%mask2dBu(I,J)) * (dvSdx(I,J,k) - duSdy(I,J,k)) * G%IareaBu(I,J)
        enddo
      endif
    endif
  else
    do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
      rel_vort(I,J,k) = G%mask2dBu(I,J) * (dvdx(I,J,k) - dudy(I,J,k)) * G%IareaBu(I,J)
    enddo

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        do concurrent (k=ksc:kec, J=Jsq-1:Jeq+1, I=Isq-1:Ieq+1)
          stk_vort(I,J,k) = (2.0 - G%mask2dBu(I,J)) * (dvSdx(I,J,k) - duSdy(I,J,k)) * G%IareaBu(I,J)
        enddo
      endif
    endif
  endif

  do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q)
    abs_vort(I,J,k) = G%CoriolisBu(I,J) + rel_vort(I,J,k)
  enddo

  do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q) DO_LOCALITY(local(hArea_q))
    hArea_q = (hArea_u(I,j,k) + hArea_u(I,j+1,k)) + (hArea_v(i,J,k) + hArea_v(i+1,J,k))
    Ih_q(I,J,k) = Area_q(I,J) / (hArea_q + vol_neglect)
    q(I,J,k) = abs_vort(I,J,k) * Ih_q(I,J,k)
  enddo

  ! NOTE: `h_q` is only used by WENO and was pulled out of the above loop to
  !   improve GPU performance, but it may need to be moved back.
  if (use_weno) then
    do concurrent (k=ksc:kec, J=Js_q:Je_q, I=Is_q:Ie_q) DO_LOCALITY(local(hArea_q))
      hArea_q = (hArea_u(I,j,k) + hArea_u(I,j+1,k)) + (hArea_v(i,J,k) + hArea_v(i+1,J,k))
      h_q(I,J,k) = hArea_q / max(Area_q(I,J), area_neglect)
    enddo
  endif

  if (Stokes_VF) then
    if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
      do concurrent (k=ksc:kec, J=js-1:Jeq, I=is-1:Ieq)
        qS(I,J,k) = stk_vort(I,J,k) * Ih_q(I,J,k)
      enddo
    endif
  endif

  if (CS%id_rv > 0) then
    do concurrent (k=ksc:kec, J=Jsq-1:Jeq+1, I=Isq-1:Ieq+1)
      RV(I,J,k) = rel_vort(I,J,k)
    enddo
  endif

  if (CS%id_PV > 0) then
    do concurrent (k=ksc:kec, J=Jsq-1:Jeq+1, I=Isq-1:Ieq+1)
      PV(I,J,k) = q(I,J,k)
    enddo
  endif

  if (associated(AD%rv_x_v) .or. associated(AD%rv_x_u)) then
    do concurrent (k=ksc:kec, J=Jsq-1:Jeq+1, I=Isq-1:Ieq+1)
      q2(I,J,k) = rel_vort(I,J,k) * Ih_q(I,J,k)
    enddo
  endif

  !$omp target exit data map(delete: dvdx, dudy)
  !$omp target exit data map(delete: hArea_u, hArea_v)
  !$omp target exit data map(delete: rel_vort)
  !$omp target exit data map(delete: dvSdx, duSdy, stk_vort) if (Stokes_VF)
end subroutine CorAdv_setup


!> Calculates the Coriolis and momentum advection accelerations with the Sadourny (1975)
!! energy- or enstrophy-conserving schemes (SADOURNY75_ENERGY, SADOURNY75_ENSTRO), over one tile.
subroutine CorAdv_sadourny(u, v, uh, vh, q, uh_center, vh_center, G, GV, bxH, bxQ, CAu, CAv, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uh !< Zonal transport u*h*dy
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vh !< Meridional transport v*h*dx
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: q !< Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: uh_center !< Transport based on arithmetic mean h at u-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: vh_center !< Transport based on arithmetic mean h at v-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: CAu !< Zonal acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: CAv !< Meridional acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    uh_min, uh_max, &   ! The smallest and largest estimates of the zonal volume fluxes through
                        ! the faces (i.e. u*h*dy) [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh_min, vh_max      ! The smallest and largest estimates of the meridional volume fluxes through
                        ! the faces (i.e. v*h*dx) [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: temp1, temp2           ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: uhc, vhc               ! Centered estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: uhm, vhm               ! The input estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: c1, c2, c3, slope      ! Nondimensional parameters for the Coriolis limiter scheme [nondim]
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, ksc, kec

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)

  ! TODO: May also need SADOURNEY75_ENERGY
  !$omp target enter data map(alloc: uh_min, vh_min) if (CS%Coriolis_En_Dis)
  !$omp target enter data map(alloc: uh_max, vh_max) if (CS%Coriolis_En_Dis)

  ! The energy-dissipating biased limiter (Coriolis_En_Dis) is only used by these schemes.
  ! .and. SADOURNEY75_ENERGY ??
  if (CS%Coriolis_En_Dis) then
  !  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
    c1 = 1.0-1.5*0.5 ; c2 = 1.0-0.5 ; c3 = 2.0 ; slope = 0.5

    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, I=is-1:ie) DO_LOCALITY(local(uhc, uhm))
      uhc = uh_center(I,j,k)
      uhm = uh(I,j,k)
      ! This sometimes matters with some types of open boundary conditions.
      if (G%dy_Cu(I,j) == 0.0) uhc = uhm

      if (abs(uhc) < 0.1*abs(uhm)) then
        uhm = 10.0*uhc
      elseif (abs(uhc) > c1*abs(uhm)) then
        if (abs(uhc) < c2*abs(uhm)) then ; uhc = (3.0*uhc+(1.0-c2*3.0)*uhm)
        elseif (abs(uhc) <= c3*abs(uhm)) then ; uhc = uhm
        else ; uhc = slope*uhc+(1.0-c3*slope)*uhm
        endif
      endif

      if (uhc > uhm) then
        uh_min(I,j,k) = uhm ; uh_max(I,j,k) = uhc
      else
        uh_max(I,j,k) = uhm ; uh_min(I,j,k) = uhc
      endif
    enddo

    do concurrent (k=ksc:kec, J=js-1:je, i=Isq:Ieq+1) DO_LOCALITY(local(vhc, vhm))
      vhc = vh_center(i,J,k)
      vhm = vh(i,J,k)
      ! This sometimes matters with some types of open boundary conditions.
      if (G%dx_Cv(i,J) == 0.0) vhc = vhm

      if (abs(vhc) < 0.1*abs(vhm)) then
        vhm = 10.0*vhc
      elseif (abs(vhc) > c1*abs(vhm)) then
        if (abs(vhc) < c2*abs(vhm)) then ; vhc = (3.0*vhc+(1.0-c2*3.0)*vhm)
        elseif (abs(vhc) <= c3*abs(vhm)) then ; vhc = vhm
        else ; vhc = slope*vhc+(1.0-c3*slope)*vhm
        endif
      endif

      if (vhc > vhm) then
        vh_min(i,J,k) = vhm ; vh_max(i,J,k) = vhc
      else
        vh_max(i,J,k) = vhm ; vh_min(i,J,k) = vhc
      endif
    enddo
  endif

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
  if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
    if (CS%Coriolis_En_Dis) then
      ! Energy dissipating biased scheme, Hallberg 200x
      do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq) DO_LOCALITY(local(temp1, temp2))
        if (q(I,J,k)*u(I,j,k) == 0.0) then
          temp1 = q(I,J,k) * ( (vh_max(i,j,k)+vh_max(i+1,j,k)) &
                           + (vh_min(i,j,k)+vh_min(i+1,j,k)) )*0.5
        elseif (q(I,J,k)*u(I,j,k) < 0.0) then
          temp1 = q(I,J,k) * (vh_max(i,j,k)+vh_max(i+1,j,k))
        else
          temp1 = q(I,J,k) * (vh_min(i,j,k)+vh_min(i+1,j,k))
        endif
        if (q(I,J-1,k)*u(I,j,k) == 0.0) then
          temp2 = q(I,J-1,k) * ( (vh_max(i,j-1,k)+vh_max(i+1,j-1,k)) &
                             + (vh_min(i,j-1,k)+vh_min(i+1,j-1,k)) )*0.5
        elseif (q(I,J-1,k)*u(I,j,k) < 0.0) then
          temp2 = q(I,J-1,k) * (vh_max(i,j-1,k)+vh_max(i+1,j-1,k))
        else
          temp2 = q(I,J-1,k) * (vh_min(i,j-1,k)+vh_min(i+1,j-1,k))
        endif
        CAu(I,j,k) = 0.25 * G%IdxCu(I,j) * (temp1 + temp2)
      enddo
    else
      ! Energy conserving scheme, Sadourny 1975
      do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
        CAu(I,j,k) = 0.25 * &
          ((q(I,J,k) * (vh(i+1,J,k) + vh(i,J,k))) + &
           (q(I,J-1,k) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
      enddo
    endif
  elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
      CAu(I,j,k) = 0.125 * (G%IdxCu(I,j) * (q(I,J,k) + q(I,J-1,k))) * &
                   ((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
    enddo
  endif

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
  if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
    if (CS%Coriolis_En_Dis) then
      ! Energy dissipating biased scheme, Hallberg 200x
      do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie) DO_LOCALITY(local(temp1, temp2))
        if (q(I-1,J,k)*v(i,J,k) == 0.0) then
          temp1 = q(I-1,J,k) * ( (uh_max(i-1,j,k)+uh_max(i-1,j+1,k)) &
                             + (uh_min(i-1,j,k)+uh_min(i-1,j+1,k)) )*0.5
        elseif (q(I-1,J,k)*v(i,J,k) > 0.0) then
          temp1 = q(I-1,J,k) * (uh_max(i-1,j,k)+uh_max(i-1,j+1,k))
        else
          temp1 = q(I-1,J,k) * (uh_min(i-1,j,k)+uh_min(i-1,j+1,k))
        endif
        if (q(I,J,k)*v(i,J,k) == 0.0) then
          temp2 = q(I,J,k) * ( (uh_max(i,j,k)+uh_max(i,j+1,k)) &
                           + (uh_min(i,j,k)+uh_min(i,j+1,k)) )*0.5
        elseif (q(I,J,k)*v(i,J,k) > 0.0) then
          temp2 = q(I,J,k) * (uh_max(i,j,k)+uh_max(i,j+1,k))
        else
          temp2 = q(I,J,k) * (uh_min(i,j,k)+uh_min(i,j+1,k))
        endif
        CAv(i,J,k) = -0.25 * G%IdyCv(i,J) * (temp1 + temp2)
      enddo
    else
      ! Energy conserving scheme, Sadourny 1975
      do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
        CAv(i,J,k) = - 0.25* &
            ((q(I-1,J,k)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
             (q(I,J,k)*(uh(I,j,k) + uh(I,j+1,k)))) * G%IdyCv(i,J)
      enddo
    endif
  elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
      CAv(i,J,k) = -0.125 * (G%IdyCv(i,J) * (q(I-1,J,k) + q(I,J,k))) * &
                   ((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))
    enddo
  endif

  !$omp target exit data map(delete: uh_min, vh_min) if (CS%Coriolis_En_Dis)
  !$omp target exit data map(delete: uh_max, vh_max) if (CS%Coriolis_En_Dis)
end subroutine CorAdv_sadourny


!> Calculates the Coriolis and momentum advection accelerations with the Arakawa & Hsu (1990),
!! Arakawa & Lamb (1981) or blended (ARAKAWA_LAMB_BLEND) schemes, over one tile.
subroutine CorAdv_arakawa(uh, vh, q, Ih_q, Fe_m2, rat_lin, G, GV, bxH, bxQ, CAu, CAv, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uh !< Zonal transport u*h*dy
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vh !< Meridional transport v*h*dx
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: q !< Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: Ih_q !< Inverse of thickness interpolated to q points
                                                   !! [H-1 ~> m-1 or m2 kg-1]
  real,                    intent(in)    :: Fe_m2   !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  real,                    intent(in)    :: rat_lin !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: CAu !< Zonal acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: CAv !< Meridional acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    a, b, c, d    ! a, b, c, & d are combinations of the potential vorticities
                  ! surrounding an h grid point.  At small scales, a = q/4,
                  ! b = q/4, etc.  All are in [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1],
                  ! and use the indexing of the corresponding u point.
  real, dimension(SZI_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)) :: &
    ep_u, ep_v    ! Additional pseudo-Coriolis terms in the Arakawa and Lamb
                  ! discretization [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, parameter :: C1_12 = 1.0 / 12.0 ! C1_12 = 1/12 [nondim]
  real, parameter :: C1_24 = 1.0 / 24.0 ! C1_24 = 1/24 [nondim]
  real :: max_Ihq, min_Ihq       ! The maximum and minimum of the nearby Ihq [H-1 ~> m-1 or m2 kg-1].
  real :: rat_m1        ! The ratio of the maximum neighboring inverse thickness
                        ! to the minimum inverse thickness minus 1 [nondim]. rat_m1 >= 0.
  real :: AL_wt         ! The relative weight of the Arakawa & Lamb scheme to the
                        ! Arakawa & Hsu scheme [nondim], between 0 and 1.
  real :: Sad_wt        ! The relative weight of the Sadourny energy scheme to
                        ! the other two with the ARAKAWA_LAMB_BLEND scheme [nondim],
                        ! between 0 and 1.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, ksc, kec

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)

  !$omp target enter data map(alloc: a, b, c, d, ep_u, ep_v)

  !   a, b, c, and d are combinations of neighboring potential
  ! vorticities which form the Arakawa and Hsu vorticity advection
  ! scheme.  All are defined at u grid points.

  if (CS%Coriolis_Scheme == ARAKAWA_HSU90) then
    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, I=is-1:Ieq)
      a(I,j,k) = (q(I,J,k) + (q(I+1,J,k) + q(I,J-1,k))) * C1_12
      d(I,j,k) = ((q(I,J,k) + q(I+1,J-1,k)) + q(I,J-1,k)) * C1_12
    enddo

    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, I=Isq:Ieq)
      b(I,j,k) = (q(I,J,k) + (q(I-1,J,k) + q(I,J-1,k))) * C1_12
      c(I,j,k) = ((q(I,J,k) + q(I-1,J-1,k)) + q(I,J-1,k)) * C1_12
    enddo
  elseif (CS%Coriolis_Scheme == ARAKAWA_LAMB81) then
    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, I=Isq:Ieq+1)
      a(I-1,j,k) = (2.0*(q(I,J,k) + q(I-1,J-1,k)) + (q(I-1,J,k) + q(I,J-1,k))) * C1_24
      d(I-1,j,k) = ((q(I,j,k) + q(I-1,J-1,k)) + 2.0*(q(I-1,J,k) + q(I,J-1,k))) * C1_24
      b(I,j,k) =   ((q(I,J,k) + q(I-1,J-1,k)) + 2.0*(q(I-1,J,k) + q(I,J-1,k))) * C1_24
      c(I,j,k) =   (2.0*(q(I,J,k) + q(I-1,J-1,k)) + (q(I-1,J,k) + q(I,J-1,k))) * C1_24
      ep_u(i,j,k) = ((q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
      ep_v(i,j,k) = (-(q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
    enddo
  elseif (CS%Coriolis_Scheme == AL_BLEND) then
    ! Fe_m2 and rat_lin are k-independent; computed in CorAdCalc_TR before the tile loop.
    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, I=Isq:Ieq+1) &
        DO_LOCALITY(local(min_Ihq, max_Ihq, rat_m1, AL_wt, Sad_wt))
      min_Ihq = MIN(Ih_q(I-1,J-1,k), Ih_q(I,J-1,k), Ih_q(I-1,J,k), Ih_q(I,J,k))
      max_Ihq = MAX(Ih_q(I-1,J-1,k), Ih_q(I,J-1,k), Ih_q(I-1,J,k), Ih_q(I,J,k))
      rat_m1 = 1.0e15
      if (max_Ihq < 1.0e15*min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0
      ! The weights used here are designed to keep the effective Coriolis
      ! acceleration from any one point on its neighbors within a factor
      ! of F_eff_max.  The minimum permitted value is 2 (the factor for
      ! Sadourny's energy conserving scheme).

      ! Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
      if (rat_m1 <= Fe_m2) then ; AL_wt = 1.0
      elseif (rat_m1 < 1.5*Fe_m2) then ; AL_wt = 3.0*Fe_m2 / rat_m1 - 2.0
      else ; AL_wt = 0.0 ; endif

      ! Determine the relative weights of Sadourny Energy vs. the other two.
      if (rat_m1 <= 1.5*Fe_m2) then ; Sad_wt = 0.0
      elseif (rat_m1 <= rat_lin) then
        Sad_wt = 1.0 - (1.5*Fe_m2) / rat_m1
      elseif (rat_m1 < 2.0*rat_lin) then
        Sad_wt = 1.0 - (CS%wt_lin_blend / rat_lin) * (rat_m1 - 2.0*rat_lin)
      else ; Sad_wt = 1.0 ; endif

      a(I-1,j,k) = Sad_wt * 0.25 * q(I-1,J,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I-1,J,k) + AL_wt*q(I,J-1,k)) + &
                       2.0 * (q(I,J,k) + q(I-1,J-1,k)) ) * C1_24
      d(I-1,j,k) = Sad_wt * 0.25 * q(I-1,J-1,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I-1,J-1,k) + AL_wt*q(I,J,k)) + &
                       2.0 * (q(I-1,J,k) + q(I,J-1,k)) ) * C1_24
      b(I,j,k) =   Sad_wt * 0.25 * q(I,J,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I,J,k) + AL_wt*q(I-1,J-1,k)) + &
                       2.0 * (q(I-1,J,k) + q(I,J-1,k)) ) * C1_24
      c(I,j,k) =   Sad_wt * 0.25 * q(I,J-1,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I,J-1,k) + AL_wt*q(I-1,J,k)) + &
                       2.0 * (q(I,J,k) + q(I-1,J-1,k)) ) * C1_24
      ep_u(i,j,k) = AL_wt  * ((q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
      ep_v(i,j,k) = AL_wt * (-(q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
    enddo
  endif

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
    ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
      CAu(I,j,k) = (((a(I,j,k) * vh(i+1,J,k)) +  (c(I,j,k) * vh(i,J-1,k)))  + &
                    ((b(I,j,k) * vh(i,J,k)) +  (d(I,j,k) * vh(i+1,J-1,k)))) * G%IdxCu(I,j)
    enddo

  ! Add in the additional terms with Arakawa & Lamb.
  if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
      (CS%Coriolis_Scheme == AL_BLEND)) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
      CAu(I,j,k) = CAu(I,j,k) + &
            ((ep_u(i,j,k)*uh(I-1,j,k)) - (ep_u(i+1,j,k)*uh(I+1,j,k))) * G%IdxCu(I,j)
    enddo
  endif

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
    ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
      CAv(i,J,k) = - (((a(I-1,j,k)   * uh(I-1,j,k)) + &
                       (c(I,j+1,k)   * uh(I,j+1,k)))  &
                    + ((b(I,j,k)     * uh(I,j,k)) +   &
                       (d(I-1,j+1,k) * uh(I-1,j+1,k)))) * G%IdyCv(i,J)
    enddo
  ! Add in the additonal terms with Arakawa & Lamb.
  if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
      (CS%Coriolis_Scheme == AL_BLEND)) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
      CAv(i,J,k) = CAv(i,J,k) + &
            ((ep_v(i,j,k)*vh(i,J-1,k)) - (ep_v(i,j+1,k)*vh(i,J+1,k))) * G%IdyCv(i,J)
    enddo
  endif

  !$omp target exit data map(delete: a, b, c, d, ep_u, ep_v)
end subroutine CorAdv_arakawa


!> Calculates the Coriolis and momentum advection accelerations with the pseudo-enstrophy
!! scheme that is robust to vanishing layers (ROBUST_ENSTRO), over one tile.
subroutine CorAdv_robust_enstro(u, v, h, uh, vh, abs_vort, eps_vel, h_tiny, G, GV, bxH, bxQ, CAu, &
                                CAv, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uh !< Zonal transport u*h*dy
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vh !< Meridional transport v*h*dx
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: abs_vort !< Absolute vorticity at q-points [T-1 ~> s-1]
  real,                    intent(in)    :: eps_vel !< A tiny, positive velocity [L T-1 ~> m s-1]
  real,                    intent(in)    :: h_tiny  !< A very small thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: CAu !< Zonal acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: CAv !< Meridional acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real :: Heff1, Heff2  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: Heff3, Heff4  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: UHeff, VHeff  ! More temporary variables [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: QUHeff,QVHeff ! More temporary variables [H L2 T-2 ~> m3 s-2 or kg s-2].
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, ksc, kec

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
    ! An enstrophy conserving scheme robust to vanishing layers
    ! Note: Heffs are in lieu of h_at_v that should be returned by the
    !       continuity solver. AJA
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq) &
        DO_LOCALITY(local(Heff1, Heff2, Heff3, Heff4, VHeff, QVHeff))
      Heff1 = abs(vh(i,J,k) * G%IdxCv(i,J)) / (eps_vel+abs(v(i,J,k)))
      Heff1 = max(Heff1, min(h(i,j,k),h(i,j+1,k)))
      Heff1 = min(Heff1, max(h(i,j,k),h(i,j+1,k)))
      Heff2 = abs(vh(i,J-1,k) * G%IdxCv(i,J-1)) / (eps_vel+abs(v(i,J-1,k)))
      Heff2 = max(Heff2, min(h(i,j-1,k),h(i,j,k)))
      Heff2 = min(Heff2, max(h(i,j-1,k),h(i,j,k)))
      Heff3 = abs(vh(i+1,J,k) * G%IdxCv(i+1,J)) / (eps_vel+abs(v(i+1,J,k)))
      Heff3 = max(Heff3, min(h(i+1,j,k),h(i+1,j+1,k)))
      Heff3 = min(Heff3, max(h(i+1,j,k),h(i+1,j+1,k)))
      Heff4 = abs(vh(i+1,J-1,k) * G%IdxCv(i+1,J-1)) / (eps_vel+abs(v(i+1,J-1,k)))
      Heff4 = max(Heff4, min(h(i+1,j-1,k),h(i+1,j,k)))
      Heff4 = min(Heff4, max(h(i+1,j-1,k),h(i+1,j,k)))
      if (CS%PV_Adv_Scheme == PV_ADV_CENTERED) then
        CAu(I,j,k) = 0.5*(abs_vort(I,J,k)+abs_vort(I,J-1,k)) * &
                     ((vh(i,J,k) + vh(i+1,J-1,k)) + (vh(i,J-1,k) + vh(i+1,J,k)) ) /  &
                     (h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) * G%IdxCu(I,j)
      elseif (CS%PV_Adv_Scheme == PV_ADV_UPWIND1) then
        VHeff = ((vh(i,J,k) + vh(i+1,J-1,k)) + (vh(i,J-1,k) + vh(i+1,J,k)) )
        QVHeff = 0.5*( ((abs_vort(I,J,k)+abs_vort(I,J-1,k))*VHeff) &
                     - ((abs_vort(I,J,k)-abs_vort(I,J-1,k))*abs(VHeff)) )
        CAu(I,j,k) = (QVHeff / ( h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) ) * G%IdxCu(I,j)
      endif
    enddo

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
    ! An enstrophy conserving scheme robust to vanishing layers
    ! Note: Heffs are in lieu of h_at_u that should be returned by the
    !       continuity solver. AJA
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie) &
        DO_LOCALITY(local(Heff1, Heff2, Heff3, Heff4, UHeff, QUHeff))
      Heff1 = abs(uh(I,j,k) * G%IdyCu(I,j)) / (eps_vel+abs(u(I,j,k)))
      Heff1 = max(Heff1, min(h(i,j,k),h(i+1,j,k)))
      Heff1 = min(Heff1, max(h(i,j,k),h(i+1,j,k)))
      Heff2 = abs(uh(I-1,j,k) * G%IdyCu(I-1,j)) / (eps_vel+abs(u(I-1,j,k)))
      Heff2 = max(Heff2, min(h(i-1,j,k),h(i,j,k)))
      Heff2 = min(Heff2, max(h(i-1,j,k),h(i,j,k)))
      Heff3 = abs(uh(I,j+1,k) * G%IdyCu(I,j+1)) / (eps_vel+abs(u(I,j+1,k)))
      Heff3 = max(Heff3, min(h(i,j+1,k),h(i+1,j+1,k)))
      Heff3 = min(Heff3, max(h(i,j+1,k),h(i+1,j+1,k)))
      Heff4 = abs(uh(I-1,j+1,k) * G%IdyCu(I-1,j+1)) / (eps_vel+abs(u(I-1,j+1,k)))
      Heff4 = max(Heff4, min(h(i-1,j+1,k),h(i,j+1,k)))
      Heff4 = min(Heff4, max(h(i-1,j+1,k),h(i,j+1,k)))
      if (CS%PV_Adv_Scheme == PV_ADV_CENTERED) then
        CAv(i,J,k) = - 0.5*(abs_vort(I,J,k)+abs_vort(I-1,J,k)) * &
                       ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                        (uh(I-1,j  ,k)+uh(I  ,j+1,k)) ) /    &
                    (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdyCv(i,J)
      elseif (CS%PV_Adv_Scheme == PV_ADV_UPWIND1) then
        UHeff = ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                 (uh(I-1,j  ,k)+uh(I  ,j+1,k)) )
        QUHeff = 0.5*( ((abs_vort(I,J,k)+abs_vort(I-1,J,k))*UHeff) &
                     - ((abs_vort(I,J,k)-abs_vort(I-1,J,k))*abs(UHeff)) )
        CAv(i,J,k) = - QUHeff / &
                     (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdyCv(i,J)
      endif
    enddo
end subroutine CorAdv_robust_enstro


!> Calculates the Coriolis and momentum advection accelerations with the WENO PV-reconstruction
!! schemes (WENOVI7TH/5TH/3RD_PV_ENSTRO), over one tile.  Near land the reconstruction falls back to
!! narrower stencils, down to first-order upwind.
subroutine CorAdv_weno(u, v, uh, vh, q, abs_vort, h_q, G, GV, bxH, bxQ, CAu, CAv, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uh !< Zonal transport u*h*dy
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vh !< Meridional transport v*h*dx
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: q !< Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: abs_vort !< Absolute vorticity at q-points [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: h_q !< Thickness interpolated to q points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: CAu !< Zonal acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: CAv !< Meridional acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real :: u_v, v_u      ! u_v is the u velocity at v point, v_u is the v velocity at u point [L T-1 ~> m s-1]
  real :: q_v, q_u      ! PV at the u and v points [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  integer :: seventh_order, fifth_order, third_order ! Order of accuracy for the WENO calculations
  real :: u_q8(8) ! Eight-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: u_q6(6) ! Six-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: u_q4(4) ! Four-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q8(8) ! Eight-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q6(6) ! Six-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q4(4) ! Four-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, ksc, kec

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq) &
        DO_LOCALITY(local(v_u, q_u, third_order, fifth_order, seventh_order, u_q8, u_q6, u_q4))
      v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      ! check whether there is masked land points in the stencil
      third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * &
                     G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))

      fifth_order   = third_order * G%mask2dCu(I,j-3) * G%mask2dCu(I,j+3)
      seventh_order = fifth_order * G%mask2dCu(I,j-4) * G%mask2dCu(I,j+4)


      ! compute the masking to make sure that inland values are not used
      if (seventh_order == 1) then
        ! all values are valid, we use seventh order reconstruction
        u_q8(:) = (u(I,j-4:j+3,k) + u(I,j-3:j+4,k)) * 0.5
        call weno_seven_h_weight_reconstruction(abs_vort(I,J-4:J+3,k), &
                                       h_q(I,J-4:J+3,k), &
                                       u_q8, &
                                       GV%H_subroundoff, v_u, q_u, cs%weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      elseif (fifth_order == 1) then
        ! all values are valid, we use fifth order reconstruction
        u_q6(:) = (u(I,j-3:j+2,k) + u(I,j-2:j+3,k)) * 0.5
        call weno_five_h_weight_reconstruction(abs_vort(I,J-3:J+2,k), &
                                      h_q(I,J-3:J+2,k), &
                                      u_q6, &
                                      GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      elseif (third_order == 1) then
        ! only the middle values are valid, we use third order reconstruction
        u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
        call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1,k), &
                                       h_q(I,J-2:J+1,k), &
                                       u_q4, &
                                       GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)
      else ! Upwind first order
        if (v_u>0.) then
            q_u = q(I,J-1,k)
        else
            q_u = q(I,J,k)
        endif
        CAu(I,j,k) = (q_u * v_u)

      endif
    enddo
  elseif (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq) &
        DO_LOCALITY(local(v_u, q_u, third_order, fifth_order, u_q6, u_q4))
      v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * &
                     G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))

      fifth_order   = third_order * G%mask2dCu(I,j-3) * G%mask2dCu(I,j+3)

      if (fifth_order == 1) then
        ! all values are valid, we use fifth order reconstruction
        u_q6(:) = (u(I,j-3:j+2,k) + u(I,j-2:j+3,k)) * 0.5
        call weno_five_h_weight_reconstruction(abs_vort(I,J-3:J+2,k), &
                                      h_q(I,J-3:J+2,k), &
                                      u_q6, &
                                      GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      elseif (third_order == 1) then
        ! only the middle values are valid, we use third order reconstruction
        u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
        call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1,k), &
                                       h_q(I,J-2:J+1,k), &
                                       u_q4, &
                                       GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      else ! Upwind first order
        if (v_u>0.) then
            q_u = q(I,J-1,k)
        else
            q_u = q(I,J,k)
        endif
        CAu(I,j,k) = (q_u * v_u)
      endif
    enddo
  elseif (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq) &
        DO_LOCALITY(local(v_u, q_u, third_order, u_q4))
      v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * &
                     G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))


      if (third_order == 1) then
        ! only the middle values are valid, we use third order reconstruction
        u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
        call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1,k), &
                                       h_q(I,J-2:J+1,k), &
                                       u_q4, &
                                       GV%H_subroundoff, v_u, q_u, CS%weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      else ! Upwind first order
        if (v_u>0.) then
            q_u = q(I,J-1,k)
        else
            q_u = q(I,J,k)
        endif
        CAu(I,j,k) = (q_u * v_u)
      endif
    enddo
  endif

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie) &
        DO_LOCALITY(local(u_v, q_v, third_order, fifth_order, seventh_order, v_q8, v_q6, v_q4))
      u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

      ! check whether there is any masked land values within the stencils
      third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * &
                     G%mask2dCv(i+2,J))
      fifth_order   = third_order * G%mask2dCv(i-3,J) * G%mask2dCv(i+3,J)
      seventh_order = fifth_order * G%mask2dCv(i-4,J) * G%mask2dCv(i+4,J)



      ! compute the masking to make sure that inland values are not used
      if (seventh_order == 1) then
        v_q8(:) = (v(i-4:i+3,J,k) + v(i-3:i+4,J,k)) * 0.5
        ! all values are valid, we use seventh order reconstruction
        call weno_seven_h_weight_reconstruction(abs_vort(I-4:I+3,J,k), &
                                       h_q(I-4:I+3,J,k), &
                                       v_q8, &
                                       GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      elseif (fifth_order == 1) then
        v_q6(:) = (v(i-3:i+2,J,k) + v(i-2:i+3,J,k)) * 0.5
        ! all values are valid, we use fifth order reconstruction
        call weno_five_h_weight_reconstruction(abs_vort(I-3:I+2,J,k), &
                                      h_q(I-3:I+2,J,k), &
                                      v_q6, &
                                      GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      elseif (third_order == 1) then
        v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
        call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J,k), &
                                               h_q(I-2:I+1,J,k), &
                                               v_q4, &
                                               GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)
      else ! Upwind first order!
        if (u_v>0.) then
            q_v = q(I-1,J,k)
        else
            q_v = q(I,J,k)
        endif
        CAv(i,J,k) = - (q_v * u_v)
      endif

    enddo
  elseif (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie) &
        DO_LOCALITY(local(u_v, q_v, third_order, fifth_order, v_q6, v_q4))
      u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

      third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * &
                     G%mask2dCv(i+2,J))
      fifth_order   = third_order * G%mask2dCv(i-3,J) * G%mask2dCv(i+3,J)


      ! compute the masking to make sure that inland values are not used
      if (fifth_order == 1) then
        v_q6(:) = (v(i-3:i+2,J,k) + v(i-2:i+3,J,k)) * 0.5
        ! all values are valid, we use fifth order reconstruction
        call weno_five_h_weight_reconstruction(abs_vort(I-3:I+2,J,k), &
                                      h_q(I-3:I+2,J,k), &
                                      v_q6, &
                                      GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      elseif (third_order == 1) then
        v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
        call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J,k), &
                                               h_q(I-2:I+1,J,k), &
                                               v_q4, &
                                               GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      else
        if (u_v>0.) then
            q_v = q(I-1,J,k)
        else
            q_v = q(I,J,k)
        endif
        CAv(i,J,k) = - (q_v * u_v)
      endif

    enddo
  elseif (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie) &
        DO_LOCALITY(local(u_v, q_v, third_order, v_q4))
      u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

      third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * &
                     G%mask2dCv(i+2,J))


      ! compute the masking to make sure that inland values are not used
      if (third_order == 1) then
        v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
        call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J,k), &
                                               h_q(I-2:I+1,J,k), &
                                               v_q4, &
                                               GV%H_subroundoff, u_v, q_v, CS%weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      else
        if (u_v>0.) then
            q_v = q(I-1,J,k)
        else
            q_v = q(I,J,k)
        endif
        CAv(i,J,k) = - (q_v * u_v)
      endif

    enddo
  endif
end subroutine CorAdv_weno


!> Adds the terms that are common to every Coriolis scheme to the accelerations over one tile:
!! the Stokes-drift diagnostic, the optional bounding of the Coriolis terms and the kinetic energy
!! gradient, and the diagnostics of the kinetic energy gradient and relative-vorticity terms.
subroutine CorAdv_common_terms(u, v, uh, vh, abs_vort, qS, q2, KEx, KEy, Stokes_VF, G, GV, bxH, &
                               bxQ, CAu, CAv, CAuS, CAvS, AD, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uh !< Zonal transport u*h*dy
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vh !< Meridional transport v*h*dx
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: abs_vort !< Absolute vorticity at q-points [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: qS !< Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: q2 !< Relative vorticity over thickness
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: KEx !< Zonal gradient of kinetic energy per unit mass [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                intent(in) :: KEy !< Meridional gradient of kinetic energy per unit mass
                                                  !! [L T-2 ~> m s-2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: CAu !< Zonal acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: CAv !< Meridional acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: CAuS !< Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: CAvS !< Stokes contribution to CAv [L T-2 ~> m s-2]
  type(accel_diag_ptrs),   intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real :: fv1, fv2, fv3, fv4   ! (f+rv)*v at the 4 points surrounding a u points[L T-2 ~> m s-2]
  real :: fu1, fu2, fu3, fu4   ! -(f+rv)*u at the 4 points surrounding a v point [L T-2 ~> m s-2]
  real :: max_fv, max_fu       ! The maximum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
  real :: min_fv, min_fu       ! The minimum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
  real, parameter :: C1_12 = 1.0 / 12.0 ! C1_12 = 1/12 [nondim]
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, ksc, kec

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)

  if (Stokes_VF) then
    if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
      ! Computing the diagnostic Stokes contribution to CAu
      do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
        CAuS(I,j,k) = 0.25 * &
              ((qS(I,J,k) * (vh(i+1,J,k) + vh(i,J,k))) + &
               (qS(I,J-1,k) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
      enddo
    endif
  endif

  if (CS%bound_Coriolis) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq) DO_LOCALITY(local(fv1, fv2, fv3, fv4, max_fv, min_fv))
      fv1 = abs_vort(I,J,k) * v(i+1,J,k)
      fv2 = abs_vort(I,J,k) * v(i,J,k)
      fv3 = abs_vort(I,J-1,k) * v(i+1,J-1,k)
      fv4 = abs_vort(I,J-1,k) * v(i,J-1,k)

      max_fv = max(fv1, fv2, fv3, fv4)
      min_fv = min(fv1, fv2, fv3, fv4)

      CAu(I,j,k) = min(CAu(I,j,k), max_fv)
      CAu(I,j,k) = max(CAu(I,j,k), min_fv)
    enddo
  endif

  ! Term - d(KE)/dx.
  do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
    CAu(I,j,k) = CAu(I,j,k) - KEx(I,j,k)
  enddo

  if (associated(AD%gradKEu)) then
    do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
      AD%gradKEu(I,j,k) = -KEx(I,j,k)
    enddo
  endif

  if (Stokes_VF) then
    if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
      ! Computing the diagnostic Stokes contribution to CAv
      do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
        CAvS(i,J,k) = 0.25 * &
              ((qS(I,J,k) * (uh(I,j+1,k) + uh(I,j,k))) + &
               (qS(I-1,J,k) * (uh(I-1,j,k) + uh(I-1,j+1,k)))) * G%IdyCv(i,J)
      enddo
    endif
  endif

  if (CS%bound_Coriolis) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie) DO_LOCALITY(local(fu1, fu2, fu3, fu4, max_fu, min_fu))
      fu1 = -abs_vort(I,J,k) * u(I,j+1,k)
      fu2 = -abs_vort(I,J,k) * u(I,j,k)
      fu3 = -abs_vort(I-1,J,k) * u(I-1,j+1,k)
      fu4 = -abs_vort(I-1,J,k) * u(I-1,j,k)

      max_fu = max(fu1, fu2, fu3, fu4)
      min_fu = min(fu1, fu2, fu3, fu4)

      CAv(I,j,k) = min(CAv(I,j,k), max_fu)
      CAv(I,j,k) = max(CAv(I,j,k), min_fu)
    enddo
  endif

  ! Term - d(KE)/dy.
  do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
    CAv(i,J,k) = CAv(i,J,k) - KEy(i,J,k)
  enddo
  if (associated(AD%gradKEv)) then
    do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
      AD%gradKEv(i,J,k) = -KEy(i,J,k)
    enddo
  endif

  if (associated(AD%rv_x_u) .or. associated(AD%rv_x_v)) then
    ! Calculate the Coriolis-like acceleration due to relative vorticity.
    if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
      if (associated(AD%rv_x_u)) then
        do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
          AD%rv_x_u(i,J,k) = - 0.25* &
            ((q2(I-1,j,k)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
             (q2(I,j,k)*(uh(I,j,k) + uh(I,j+1,k)))) * G%IdyCv(i,J)
        enddo
      endif

      if (associated(AD%rv_x_v)) then
        do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
          AD%rv_x_v(I,j,k) = 0.25 * &
            ((q2(I,j,k) * (vh(i+1,J,k) + vh(i,J,k))) + &
             (q2(I,j-1,k) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
        enddo
      endif
    else
      if (associated(AD%rv_x_u)) then
        do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
          AD%rv_x_u(i,J,k) = -G%IdyCv(i,J) * C1_12 * &
            (((((q2(I,J,k) + q2(I-1,J-1,k)) + q2(I-1,J,k)) * uh(I-1,j,k)) + &
              (((q2(I-1,J,k) + q2(I,J+1,k)) + q2(I,J,k)) * uh(I,j+1,k))) + &
             ((((q2(I-1,J,k) + q2(I,J-1,k)) + q2(I,J,k)) * uh(I,j,k))+ &
              (((q2(I,J,k) + q2(I-1,J+1,k)) + q2(I-1,J,k)) * uh(I-1,j+1,k))))
        enddo
      endif

      if (associated(AD%rv_x_v)) then
        do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
          AD%rv_x_v(I,j,k) = G%IdxCu(I,j) * C1_12 * &
            (((((q2(I+1,J,k) + q2(I,J-1,k)) + q2(I,J,k)) * vh(i+1,J,k)) + &
              (((q2(I-1,J-1,k) + q2(I,J,k)) + q2(I,J-1,k)) * vh(i,J-1,k))) + &
             ((((q2(I-1,J,k) + q2(I,J-1,k)) + q2(I,J,k)) * vh(i,J,k)) + &
              (((q2(I+1,J-1,k) + q2(I,J,k)) + q2(I,J-1,k)) * vh(i+1,J-1,k))))
        enddo
      endif
    endif
  endif
end subroutine CorAdv_common_terms


!> Offers the Coriolis-related derived quantities for averaging, once the accelerations over
!! every tile have been calculated.
subroutine CorAdv_finalize_diagnostics(RV, PV, CAuS, CAvS, Stokes_VF, AD, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), intent(in) :: RV !< Diagnostic relative vorticity [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), intent(in) :: PV !< Diagnostic potential vorticity
                                                                !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: CAuS !< Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: CAvS !< Stokes contribution to CAv [L T-2 ~> m s-2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  type(accel_diag_ptrs),   intent(in)    :: AD  !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  integer :: nz

  nz = GV%ke

  ! Here the various Coriolis-related derived quantities are offered for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_rv > 0) call post_data(CS%id_rv, RV, CS%diag)
    if (CS%id_PV > 0) call post_data(CS%id_PV, PV, CS%diag)
    if (CS%id_gKEu>0) call post_data(CS%id_gKEu, AD%gradKEu, CS%diag)
    if (CS%id_gKEv>0) call post_data(CS%id_gKEv, AD%gradKEv, CS%diag)
    if (CS%id_rvxu > 0) call post_data(CS%id_rvxu, AD%rv_x_u, CS%diag)
    if (CS%id_rvxv > 0) call post_data(CS%id_rvxv, AD%rv_x_v, CS%diag)
    if (Stokes_VF) then
      if (CS%id_CAuS > 0) call post_data(CS%id_CAuS, CAuS, CS%diag)
      if (CS%id_CAvS > 0) call post_data(CS%id_CAvS, CAvS, CS%diag)
    endif

    ! Diagnostics for terms multiplied by fractional thicknesses

    ! 3D diagnostics hf_gKEu etc. are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    ! if (CS%id_hf_gKEu > 0) call post_product_u(CS%id_hf_gKEu, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_gKEv > 0) call post_product_v(CS%id_hf_gKEv, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    ! if (CS%id_hf_rvxv > 0) call post_product_u(CS%id_hf_rvxv, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_rvxu > 0) call post_product_v(CS%id_hf_rvxu, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_hf_gKEu_2d > 0) call post_product_sum_u(CS%id_hf_gKEu_2d, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_gKEv_2d > 0) call post_product_sum_v(CS%id_hf_gKEv_2d, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    if (CS%id_intz_gKEu_2d > 0) call post_product_sum_u(CS%id_intz_gKEu_2d, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_gKEv_2d > 0) call post_product_sum_v(CS%id_intz_gKEv_2d, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_hf_rvxv_2d > 0) call post_product_sum_u(CS%id_hf_rvxv_2d, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_rvxu_2d > 0) call post_product_sum_v(CS%id_hf_rvxu_2d, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_h_gKEu > 0) call post_product_u(CS%id_h_gKEu, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_gKEv > 0) call post_product_v(CS%id_h_gKEv, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)
    if (CS%id_h_rvxv > 0) call post_product_u(CS%id_h_rvxv, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_rvxu > 0) call post_product_v(CS%id_h_rvxu, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_intz_rvxv_2d > 0) call post_product_sum_u(CS%id_intz_rvxv_2d, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_rvxu_2d > 0) call post_product_sum_v(CS%id_intz_rvxu_2d, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)
  endif
end subroutine CorAdv_finalize_diagnostics


!> Calculates the acceleration due to the gradient of kinetic energy over one iteration tile.
subroutine gradKE(u, v, h, KE, KEx, KEy, bxH, bxQ, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(Box_t),                                intent(in)  :: bxH !< The h-point iteration box of this
                                                                 !! tile, [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),                                intent(in)  :: bxQ !< The B-grid-index iteration box of
                                                                 !! this tile, [IscB:IecB, JscB:JecB, ksc:kec]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: u   !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: v   !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)), &
                                              intent(out) :: KE  !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJ_(G),bxH%idxS(3):bxH%idxE(3)), &
                                              intent(out) :: KEx !< Zonal acceleration due to kinetic
                                                                 !! energy gradient [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)), &
                                              intent(out) :: KEy !< Meridional acceleration due to kinetic
                                                                 !! energy gradient [L T-2 ~> m s-2]
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)  :: CS  !< Control structure for MOM_CoriolisAdv
  ! Local variables
  real :: um, up, vm, vp         ! Temporary variables [L T-1 ~> m s-1].
  real :: um2, up2, vm2, vp2     ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: um2a, up2a, vm2a, vp2a ! Temporary variables [L4 T-2 ~> m4 s-2].
  real :: third_order_u, third_order_v  ! Product of mask values to determine the boundary
  integer :: i, j, k, ksc, kec, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  real, parameter     :: C1_12 = 1.0/12.0   ! The ratio of 1/12 [nondim]

  is = bxH%idxS(1) ; ie = bxH%idxE(1) ; js = bxH%idxS(2) ; je = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)

  ! Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
  if (CS%KE_Scheme == KE_ARAKAWA) then
    ! The following calculation of Kinetic energy includes the metric terms
    ! identified in Arakawa & Lamb 1982 as important for KE conservation.  It
    ! also includes the possibility of partially-blocked tracer cell faces.
    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, i=Isq:Ieq+1)
      KE(i,j,k) = ( ( (G%areaCu( I ,j)*(u( I ,j,k)*u( I ,j,k))) + &
                       (G%areaCu(I-1,j)*(u(I-1,j,k)*u(I-1,j,k))) ) + &
                     ( (G%areaCv(i, J )*(v(i, J ,k)*v(i, J ,k))) + &
                       (G%areaCv(i,J-1)*(v(i,J-1,k)*v(i,J-1,k))) ) )*0.25*G%IareaT(i,j)
    enddo
  elseif (CS%KE_Scheme == KE_SIMPLE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme which does not take into account any geometric factors
    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, i=Isq:Ieq+1) &
        DO_LOCALITY(local(up, um, vp, vm, up2, um2, vp2, vm2))
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2 = up*up
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2 = um*um
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2 = vp*vp
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2 = vm*vm
      KE(i,j,k) = ( max(up2,um2) + max(vp2,vm2) ) *0.5
    enddo
  elseif (CS%KE_Scheme == KE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme but has been adapted to take horizontal grid factors into account
    do concurrent (k=ksc:kec, j=Jsq:Jeq+1, i=Isq:Ieq+1) &
        DO_LOCALITY(local(up, um, vp, vm, up2a, um2a, vp2a, vm2a))
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2a = up*up*G%areaCu(I-1,j)
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2a = um*um*G%areaCu( I ,j)
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2a = vp*vp*G%areaCv(i,J-1)
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2a = vm*vm*G%areaCv(i, J )
      KE(i,j,k) = ( max(um2a,up2a) + max(vm2a,vp2a) )*0.5*G%IareaT(i,j)
    enddo
  elseif (CS%KE_Scheme == KE_UP3) then
    ! The following discretization of KE is based on the one-dimensional third-order
    ! upwind scheme which does not take horizontal grid factors into account

    if (CS%KE_use_limiter) then
      do concurrent (k=ksc:kec, j=Jsq:Jeq+1, i=Isq:Ieq+1) &
          DO_LOCALITY(local(up, um, vp, vm, third_order_u, third_order_v))
        ! compute the masking to make sure that inland values are not used
        third_order_u = (G%mask2dCu(I-2,j) * G%mask2dCu(I-1,j)* &
                       G%mask2dCu(I,j) * G%mask2dCu(I+1,j))

        if (third_order_u == 1) then
          up = (7.0 * (u(I-1,j,k) + u(I,j,k)) - (u(I-2,j,k) + u(I+1,j,k))) * C1_12
          call UP3_Koren_limiter_reconstruction(u(I-2:I+1,j,k), up, um)
        else
          up = (u(I-1,j,k) + u(I,j,k))*0.5
          if (up>0.) then
            um = u(I-1,j,k)
          elseif (up<0.) then
            um = u(I,j,k)
          else
            um = up
          endif
        endif

        third_order_v = (G%mask2dCv(i,J-2) * G%mask2dCv(i,J-1)* &
                       G%mask2dCv(i,J) * G%mask2dCv(i,J+1))
        if (third_order_v ==1) then
          vp = (7.0 * (v(i,J-1,k) + v(i,J,k)) - (v(i,J-2,k) + v(i,J+1,k))) * C1_12
          call UP3_Koren_limiter_reconstruction(v(i,J-2:J+1,k), vp, vm)
        else
          vp = (v(i,J-1,k) + v(i,J,k))*0.5
          if (vp>0.) then
            vm = v(i,J-1,k)
          elseif (vp<0.) then
            vm = v(i,J,k)
          else
            vm = vp
          endif
        endif

        KE(i,j,k) = ( (um*um) + (vm*vm) )*0.5
      enddo
    else
      do concurrent (k=ksc:kec, j=Jsq:Jeq+1, i=Isq:Ieq+1) &
          DO_LOCALITY(local(up, um, vp, vm, third_order_u, third_order_v))
        ! compute the masking to make sure that inland values are not used
        third_order_u = (G%mask2dCu(I-2,j) * G%mask2dCu(I-1,j)* &
                       G%mask2dCu(I,j) * G%mask2dCu(I+1,j))

        if (third_order_u == 1) then
          up = (7.0 * (u(I-1,j,k) + u(I,j,k)) - (u(I-2,j,k) + u(I+1,j,k))) * C1_12
          call UP3_reconstruction(u(I-2:I+1,j,k), up, um)
        else
          up = (u(I-1,j,k) + u(I,j,k))*0.5
          if (up>0.) then
            um = u(I-1,j,k)
          elseif (up<0.) then
            um = u(I,j,k)
          else
            um = up
          endif
        endif

        third_order_v = (G%mask2dCv(i,J-2) * G%mask2dCv(i,J-1)* &
                       G%mask2dCv(i,J) * G%mask2dCv(i,J+1))
        if (third_order_v ==1) then
          vp = (7.0 * (v(i,J-1,k) + v(i,J,k)) - (v(i,J-2,k) + v(i,J+1,k))) * C1_12
          call UP3_reconstruction(v(i,J-2:J+1,k), vp, vm)
        else
          vp = (v(i,J-1,k) + v(i,J,k))*0.5
          if (vp>0.) then
            vm = v(i,J-1,k)
          elseif (vp<0.) then
            vm = v(i,J,k)
          else
            vm = vp
          endif
        endif

        KE(i,j,k) = ( (um*um) + (vm*vm) )*0.5
      enddo
    endif
  endif

  ! Term - d(KE)/dx.
  do concurrent (k=ksc:kec, j=js:je, I=Isq:Ieq)
    KEx(I,j,k) = (KE(i+1,j,k) - KE(i,j,k)) * G%IdxCu_OBCmask(I,j)
  enddo

  ! Term - d(KE)/dy.
  do concurrent (k=ksc:kec, J=Jsq:Jeq, i=is:ie)
    KEy(i,J,k) = (KE(i,j+1,k) - KE(i,j,k)) * G%IdyCv_OBCmask(i,J)
  enddo
end subroutine gradKE

!> Reconstruct the scalar (e.g., pv, vorticity) onto point i-1/2 using a third-order upwind scheme
pure subroutine UP3_reconstruction(q4,u,qr)
  real, intent(in)    :: q4(4)            !< Tracer values on points i-2, i-1, i, i+1 [A ~> a]
  real, intent(in)    :: u                !< Velocity or thickness flux on point i-1/2
                                          !! [l t-1 ~> m s-1] or [l2 t-1 ~> m2 s-1]
  real, intent(inout) :: qr               !< Reconstruction of tracer q at point i-1/2 [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0       ! The ratio of 1/6 [nondim]

  if (u>0.) then
    qr = ((2.*q4(3) + 5.*q4(2)) - q4(1)) * C1_6
  else
    qr = ((2.*q4(2) + 5.*q4(3)) - q4(4)) * C1_6
  endif

end subroutine UP3_reconstruction


!> Reconstruct the scalar (e.g., PV, vorticity) onto point i-1/2
!! using a third-order upwind scheme with the Koren flux limiter
pure subroutine UP3_Koren_limiter_reconstruction(q4,u,qr)
  real, intent(in)    :: q4(4)            !< Tracer values on points i-2, i-1, i, i+1 [A ~> a]
  real, intent(in)    :: u                !< Velocity or thickness flux on point i-1/2
                                          !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
  real, intent(inout) :: qr               !< Reconstruction of tracer q on point i-1/2 [A ~> a]
  real                :: theta            ! Ratio of gradient [nondim]
  real                :: psi              ! Limiter function [nondim]
  real, parameter     :: C1_3 = 1.0/3.0   ! The ratio of 1/3 [nondim]
  real, parameter     :: C1_6 = 1.0/6.0   ! The ratio of 1/6 [nondim]

  if (u>0.) then
    if (q4(3) == q4(2)) then
      qr = q4(2)
    else
      theta = (q4(2) - q4(1))/(q4(3) - q4(2))
      psi = max(0., min(1., C1_3 + C1_6*theta, theta)) ! limiter introduced by Koren (1993)
      qr = q4(2) + psi*(q4(3) - q4(2))
    endif
  else
    if (q4(3) == q4(2)) then
      qr = q4(3)
    else
      theta = (q4(4) - q4(3))/(q4(3) - q4(2))
      psi = max(0., min(1., C1_3 + C1_6*theta, theta))
      qr = q4(3) + psi*(q4(2) - q4(3))
    endif
  endif

end subroutine UP3_Koren_limiter_reconstruction

!> Compute the factor for the WENO weights
pure function fac_fn(tau, b) result(fac)
  real, intent(in)  :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in)  :: b    !< The smoothness indicator [A ~> a]
  real :: fac               !< The factor for the weight [nondim]

  fac = 1.0e40 ; if (abs(b) > 1.0e-20*tau) fac = (1 + tau / b)**2

end function fac_fn


!> Reconstruct the tracer (e.g., PV, vorticity) onto the point i-1/2 using a third-order WENO scheme
!! This reconstruction is thickness-weighted
pure subroutine weno_three_h_weight_reconstruction(q4, h4, u4, &
                                              h_tiny, u, qr, velocity_smoothing)
    real, intent(in)    :: q4(4)   !< Tracer value times thickness on points i-2, i-1, i, i+1 [A ~> a]
    real, intent(in)    :: h4(4)   !< Thickness values on points i-2, i-1, i, i+1 [L ~> m]
    real, optional, intent(in)    :: u4(4) !< Velocity values on points i-2, i-1, i, i+1
                                                    !![L T-1 ~> m s-1]
    real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
    real, intent(in)    :: u              !< Velocity or thickness flux on point i-1/2
                                          !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
    real, intent(inout) :: qr             !< Reconstruction of tracer q on point i-1/2 [A ~> a]
    logical, intent(in) :: velocity_smoothing !< If true, use velocity to compute smoothness indicator
    real :: vr                            ! Reconstruction of hq [A ~> a]
    real :: hr                            ! Reconstruction of h [L ~> m]
    real :: c0, c1                        ! Intermediate reconstruction of q [A ~> a]
    real :: d0, d1                        ! Intermediate reconstruction of h [L ~> m]
    real :: b0, b1                        ! Smoothness indicator [A ~> a]
    real :: tau                           ! Difference of smoothness indicator [A ~> a]
    real :: w0, w1                        ! Weights [nondim]
    real :: s                             ! Temporary variables [nondim]
    real, parameter :: C2_3 = 2.0/3.0     ! The ratio of 2/3 [nondim]
    real, parameter :: C1_3 = 1.0/3.0     ! The ratio of 1/3 [nondim]

    if (u>0.) then
      call weno_three_reconstruction_0(q4(2:3), c0) ! Reconstruction in the second upwind stencil
      call weno_three_reconstruction_1(q4(1:2), c1) ! Reconstruction in the first upwind stencil

      call weno_three_reconstruction_0(h4(2:3), d0)
      call weno_three_reconstruction_1(h4(1:2), d1)
      if (velocity_smoothing) then
        call weno_three_weight(u4(2:3), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(u4(1:2), b1)  ! Smoothness indicator the first upwind stencil
      else
        call weno_three_weight(q4(2:3), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(q4(1:2), b1)  ! Smoothness indicator the first upwind stencil
      endif
    else
      call weno_three_reconstruction_0(q4(3:2:-1), c0) ! Reconstruction in the second upwind stencil
      call weno_three_reconstruction_1(q4(4:3:-1), c1) ! Reconstruction in the first upwind stencil

      call weno_three_reconstruction_0(h4(3:2:-1), d0)
      call weno_three_reconstruction_1(h4(4:3:-1), d1)
      if (velocity_smoothing) then
        call weno_three_weight(u4(3:2:-1), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(u4(4:3:-1), b1)  ! Smoothness indicator the first upwind stencil
      else
        call weno_three_weight(q4(3:2:-1), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(q4(4:3:-1), b1)  ! Smoothness indicator the first upwind stencil
      endif
    endif

    tau = abs(b0-b1)
    w0  = C2_3 * fac_fn(tau, b0)
    w1  = C1_3 * fac_fn(tau, b1)

    s = 1. / (w0 + w1)
    w0 = w0 * s   ! Weights of stencils
    w1 = w1 * s

    vr = (w0 * c0) + (w1 * c1)
    hr = (w0 * d0) + (w1 * d1)
!    vr = min(max(q4(3), q4(2)), vr) ; vr = max(min(q4(3), q4(2)), vr) !Impose a monotonicity limiter
    hr = min(max(h4(3), h4(2)), hr) ; hr = max(min(h4(3), h4(2)), hr) ! A monotonicity limiter

    qr = vr / max(hr, h_tiny)

end subroutine weno_three_h_weight_reconstruction

!> Compute the smoothness indicator for the two-point stencil of the third-order WENO scheme
pure subroutine weno_three_weight(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Smoothness indicator for this stencil [A2 ~> a2]

    w0 = (q2(1) - q2(2))**2

end subroutine weno_three_weight

!> Reconstruction in the second upwind stencil of the third-order WENO scheme
pure subroutine weno_three_reconstruction_0(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Reconstruction of the quantity [A2 ~> a2]

    w0 = (q2(1) + q2(2)) * 0.5

end subroutine weno_three_reconstruction_0

!> Reconstruction in the first upwind stencil for third-order WENO scheme
pure subroutine weno_three_reconstruction_1(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Reconstruction of the quantity [A ~> a]

    w0 = (- q2(1) + 3 * q2(2)) * 0.5

end subroutine weno_three_reconstruction_1


!> Reconstruct the tracer (e.g., PV, vorticity) onto point i-1/2 using a fifth-order WENO scheme
!! The reconstruction is weighted by the thickness
pure subroutine weno_five_h_weight_reconstruction(q6, h6, u6, &
                                             h_tiny, u, qr, velocity_smoothing)
    real, intent(in)    :: q6(6)
    !< Tracer values on points i-3, i-2, i-1, i, i+1, i+2 [A ~> a]
    real, intent(in)    :: h6(6)
    !< Thickness values on points i-3, i-2, i-1, i, i+1, i+2 [L ~> m]
    real, optional, intent(in)    :: u6(6)
    !< Velocity values on points i-3, i-2, i-1, i, i+1, i+2 [L T-1 ~> m s-1]
    real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
    real, intent(in)    :: u                      !< Velocity or thickness flux on point i-1/2
                                                  !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
    logical, intent(in) :: velocity_smoothing     !< If ture, use velocity to compute the smoothness indicator
    real, intent(inout) :: qr                     !< Reconstruction of tracer q on point i-1/2 [A ~> a]
    real :: vr                                    ! Reconstruction of hq [A ~> a]
    real :: hr                                    ! Reconstruction of h [L ~> m]
    real :: c0, c1, c2                            ! Intermediate reconstruction of hq[A ~> a]
    real :: d0, d1, d2                            ! Intermediate reconstruction of h [L ~> m]
    real :: b0, b1, b2                            ! Smoothness indicator [A ~> a]
    real :: tau                                   ! Difference of smoothness indicators [A ~> a]
    real :: w0, w1, w2                            ! Weights [nondim]
    real :: s                                     ! Temporary variables [nondim]
    real, parameter :: C3_10 = 3.0/10.0           ! The ratio of 3/10 [nondim]
    real, parameter :: C3_5 = 3.0/5.0             ! The ratio of 3/5 [nondim]
    real, parameter :: C1_10 = 1.0/10.0           ! The ratio of 1/10 [nondim]

    if (u>0.) then
      call weno_five_reconstruction_0(q6(3:5), c0) ! Reconstruction in the third upwind stencil
      call weno_five_reconstruction_1(q6(2:4), c1) ! Reconstruction in the second upwind stencil
      call weno_five_reconstruction_2(q6(1:3), c2) ! Reconstruction in the first upwind stencil

      call weno_five_reconstruction_0(h6(3:5), d0)
      call weno_five_reconstruction_1(h6(2:4), d1)
      call weno_five_reconstruction_2(h6(1:3), d2)
      if (velocity_smoothing) then
        call weno_five_weight_0(u6(3:5), b0)   ! Smoothness indicator of the third upwind stencil
        call weno_five_weight_1(u6(2:4), b1)   ! Smoothness indicator of the second upwind stencil
        call weno_five_weight_2(u6(1:3), b2)   ! Smoothness indicator of the first upwind stencil
      else
        call weno_five_weight_0(q6(3:5), b0)
        call weno_five_weight_1(q6(2:4), b1)
        call weno_five_weight_2(q6(1:3), b2)
      endif
    else
      call weno_five_reconstruction_0(q6(4:2:-1), c0) ! Reconstruction in the third upwind stencil
      call weno_five_reconstruction_1(q6(5:3:-1), c1) ! Reconstruction in the second upwind stencil
      call weno_five_reconstruction_2(q6(6:4:-1), c2) ! Reconstruction in the first upwind stencil

      call weno_five_reconstruction_0(h6(4:2:-1), d0)
      call weno_five_reconstruction_1(h6(5:3:-1), d1)
      call weno_five_reconstruction_2(h6(6:4:-1), d2)
      if (velocity_smoothing) then
        call weno_five_weight_0(u6(4:2:-1), b0) ! Smoothness indicator of the third upwind stencil
        call weno_five_weight_1(u6(5:3:-1), b1) ! Smoothness indicator of the second upwind stencil
        call weno_five_weight_2(u6(6:4:-1), b2) ! Smoothness indicator of the first upwind stencil
      else
        call weno_five_weight_0(q6(4:2:-1), b0)
        call weno_five_weight_1(q6(5:3:-1), b1)
        call weno_five_weight_2(q6(6:4:-1), b2)
      endif
    endif

    tau = abs(b0 - b2)
    w0  = C3_10 * fac_fn(tau, b0)
    w1  = C3_5  * fac_fn(tau, b1)
    w2  = C1_10 * fac_fn(tau, b2)

    s = 1. / ((w0 + w1) + w2)
    w0 = w0 * s   ! Weights of stencils
    w1 = w1 * s
    w2 = w2 * s

    vr = ((w0 * c0) + (w1 * c1)) + (w2 * c2)
    hr = ((w0 * d0) + (w1 * d1)) + (w2 * d2)
!    vr = min(max(q6(3), q6(4)), vr) ; vr = max(min(q6(3), q6(4)), vr) !Impose a monotonicity limiter
    hr = min(max(h6(3), h6(4)), hr) ; hr = max(min(h6(3), h6(4)), hr) !Impose a monotonicity limiter

    qr = vr / max(hr, h_tiny)

end subroutine weno_five_h_weight_reconstruction

!> Compute the smoothness indicator for the third upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_weight_0(q3, w0)
  real, intent(in) :: q3(3)       !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w0       !< Smoothness indicator for this stencil [A2 ~> a2]

  w0 = (q3(1) * ((10 * q3(1) - 31 * q3(2)) + 11 * q3(3))) + &
       ((q3(2) * (25 * q3(2) - 19 * q3(3))) + 4 * (q3(3) * q3(3)))

end subroutine weno_five_weight_0

!> Compute the smoothness indicator for the second upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_weight_1(q3, w1)
  real, intent(in) :: q3(3)        !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w1        !< Smoothness indicator for this stencil [A2 ~> a2]

  w1 = (q3(1) * ((4 * q3(1) - 13 * q3(2)) + 5 * q3(3))) + &
       ((q3(2) * (13 * q3(2) - 13 * q3(3))) + 4 * (q3(3) * q3(3)))

end subroutine weno_five_weight_1

!> Compute the smoothness indicator for the first upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_weight_2(q3, w2)
  real, intent(in) :: q3(3)        !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w2        !< Smoothness indicator for this stencil [A2 ~> a2]

  w2 = (q3(1) * ((4 * q3(1) - 19 * q3(2)) + 11 * q3(3))) + &
       ((q3(2) * (25 * q3(2) - 31 * q3(3))) + 10 * (q3(3) * q3(3)))

end subroutine weno_five_weight_2

!> Reconstruction in the third upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_reconstruction_0(q3, p0)
  real, intent(in) :: q3(3)        !< Tracer values on three points [A ~> a]
  real, intent(inout) :: p0        !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0 ! One sixth [nondim]

  p0 = ((2*q3(1) + 5*q3(2)) - q3(3)) * C1_6

end subroutine weno_five_reconstruction_0

!> Reconstruction in the second upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_reconstruction_1(q3, p1)
  real, intent(in) :: q3(3)         !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: p1         !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0 ! One sixth [nondim]

  p1 = ((-q3(1) + 5*q3(2)) + 2*q3(3)) * C1_6

end subroutine weno_five_reconstruction_1

!> Reconstruction in the first upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_reconstruction_2(q3, p2)
  real, intent(in) :: q3(3)          !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: p2          !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0  ! One sixth [nondim]

  p2 = ((2*q3(1) - 7*q3(2)) + 11*q3(3)) * C1_6

end subroutine weno_five_reconstruction_2


!> Reconstruct the tracer (e.g., PV, vorticity) onto point i-1/2 using a seventh-order WENO scheme
!! This reconstruction computes a thickness weighted average of PV
pure subroutine weno_seven_h_weight_reconstruction(q8, h8, u8, &
                                            h_tiny, u, qr, velocity_smoothing)
  real, intent(in)    :: q8(8)
  !< Tracer values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3
  real, intent(in)    :: h8(8)
  !< Thickness on the same tracer points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L ~> m]
  real, optional, intent(in)    :: u8(8)
  !< Velocity values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L T-1 ~> m s-1]
  real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
  real, intent(in)    :: u    !< Velocity or thickness flux on point i-1/2
                              !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
  logical, intent(in) :: velocity_smoothing !< If true, use velocity to compute the smoothness indicator
  real, intent(inout) :: qr   !< Reconstruction of tracer q on point i-1/2 [A ~> a]
  real :: vr                  ! Reconstruction of hq [A ~> a]
  real :: hr                  ! Reconstruction of h [L ~> m]
  real :: c0, c1, c2, c3      ! Intermediate reconstruction of hq [A ~> a]
  real :: d0, d1, d2, d3      ! Intermediate reconstruction of h [L ~> m]
  real :: b0, b1, b2, b3      ! Smoothness indicator [A ~> a]
  real :: tau                 ! Difference of smoothness indicators [A ~> a]
  real :: w0, w1, w2, w3      ! Weights [nondim]
  real :: s                   ! Temporary variables [nondim]
  real, parameter :: C4_35 = 4.0/35.0 ! The ratio of 4/35 [nondim]
  real, parameter :: C18_35 = 18.0/35.0 ! The ratio of 18/35 [nondim]
  real, parameter :: C12_35 = 12.0/35.0 ! The ratio of 12/35 [nondim]
  real, parameter :: C1_35 = 1.0/35.0   ! The ratio of 1/35 [nondim]

  if (u>0.) then
    call weno_seven_reconstruction_0(q8(4:7), c0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(q8(3:6), c1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(q8(2:5), c2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(q8(1:4), c3) ! Reconstruction in the first upwind stencil

    call weno_seven_reconstruction_0(h8(4:7), d0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(h8(3:6), d1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(h8(2:5), d2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(h8(1:4), d3) ! Reconstruction in the first upwind stencil
    if (velocity_smoothing) then
      call weno_seven_weight_0(u8(4:7), b0)       ! Smoothness indicator of the fourth upwind stencil
      call weno_seven_weight_1(u8(3:6), b1)       ! Smoothness indicator of the third upwind stencil
      call weno_seven_weight_2(u8(2:5), b2)       ! Smoothness indicator of the second upwind stencil
      call weno_seven_weight_3(u8(1:4), b3)       ! Smoothness indicator of the first upwind stencil
    else
      call weno_seven_weight_0(q8(4:7), b0)
      call weno_seven_weight_1(q8(3:6), b1)
      call weno_seven_weight_2(q8(2:5), b2)
      call weno_seven_weight_3(q8(1:4), b3)
    endif
  else
    call weno_seven_reconstruction_0(q8(5:2:-1), c0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(q8(6:3:-1), c1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(q8(7:4:-1), c2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(q8(8:5:-1), c3) ! Reconstruction in the first upwind stencil

    call weno_seven_reconstruction_0(h8(5:2:-1), d0)
    call weno_seven_reconstruction_1(h8(6:3:-1), d1)
    call weno_seven_reconstruction_2(h8(7:4:-1), d2)
    call weno_seven_reconstruction_3(h8(8:5:-1), d3)
    if (velocity_smoothing) then
      call weno_seven_weight_0(u8(5:2:-1), b0)    ! Smoothness indicator of the fourth upwind stencil
      call weno_seven_weight_1(u8(6:3:-1), b1)    ! Smoothness indicator of the third upwind stencil
      call weno_seven_weight_2(u8(7:4:-1), b2)    ! Smoothness indicator of the second upwind stencil
      call weno_seven_weight_3(u8(8:5:-1), b3)    ! Smoothness indicator of the first upwind stencil
    else
      call weno_seven_weight_0(q8(5:2:-1), b0)
      call weno_seven_weight_1(q8(6:3:-1), b1)
      call weno_seven_weight_2(q8(7:4:-1), b2)
      call weno_seven_weight_3(q8(8:5:-1), b3)
    endif
  endif

  tau = abs((b0 - b3) + 3 * (b1 - b2))
  w0  = C4_35  * fac_fn(tau, b0)
  w1  = C18_35 * fac_fn(tau, b1)
  w2  = C12_35 * fac_fn(tau, b2)
  w3  = C1_35  * fac_fn(tau, b3)

  s = 1. / ((w0 + w1) + (w2 + w3))
  w0 = w0 * s   ! Weights of the stencils
  w1 = w1 * s
  w2 = w2 * s
  w3 = w3 * s

  vr = ((w0 * c0) + (w1 * c1)) + ((w2 * c2) + (w3 * c3))
  hr = ((w0 * d0) + (w1 * d1)) + ((w2 * d2) + (w3 * d3))

!  vr = min(max(q4, q5), vr) ; vr = max(min(q4, q5), vr)
  hr = min(max(h8(4), h8(5)), hr) ; hr = max(min(h8(4), h8(5)), hr) ! Impose a monotonicity limiter

  qr = vr / max(hr, h_tiny)

end subroutine weno_seven_h_weight_reconstruction

!> Compute the smoothness indicator for the fourth upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_0(q4, w0)
  real, intent(in) :: q4(4)          !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w0          !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w0 = ((q4(1) * ((2.107 * q4(1) - 9.402 * q4(2)) + (7.042 * q4(3) - 1.854 * q4(4)))) + &
     (q4(2) * ((11.003 * q4(2) - 17.246 * q4(3)) + 4.642 * q4(4)))) + &
     ((q4(3) * (7.043 * q4(3) - 3.882 * q4(4))) + 0.547 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_0

!> Compute the smoothness indicator for the third upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_1(q4, w1)
  real, intent(in) :: q4(4)          !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w1          !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w1 = ((q4(1) * ((0.547 * q4(1) - 2.522 * q4(2)) + (1.922 * q4(3) - 0.494 * q4(4)))) + &
     (q4(2) * ((3.443 * q4(2) - 5.966 * q4(3)) + 1.602 * q4(4)))) + &
     ((q4(3) * (2.843 * q4(3) - 1.642 * q4(4))) + 0.267 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_1

!> Compute the smoothness indicator for the second upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_2(q4, w2)
  real, intent(in) :: q4(4)           !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w2           !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w2 = ((q4(1) * ((0.267 * q4(1) - 1.642 * q4(2)) + (1.602 * q4(3) - 0.494 * q4(4)))) + &
     (q4(2) * ((2.843 * q4(2) - 5.966 * q4(3)) + 1.922 * q4(4)))) + &
     ((q4(3) * (3.443 * q4(3) - 2.522 * q4(4))) + 0.547 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_2

!> Compute smoothness indicator for the first upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_3(q4, w3)
  real, intent(in) :: q4(4)           !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w3           !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w3 = ((q4(1) * ((0.547  * q4(1) - 3.882 * q4(2)) + (4.642 * q4(3) - 1.854 * q4(4)))) + &
     (q4(2) * ((7.043 * q4(2) - 17.246 * q4(3)) + 7.042 * q4(4)))) + &
     ((q4(3) * (11.003 * q4(3) - 9.402 * q4(4))) + 2.107 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_3

!> Reconstruction in the fourth upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_0(q4, p0)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p0            !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]

  p0 = (((6 * q4(1) + 26 * q4(2)) - 10 * q4(3)) + 2 * q4(4)) * C1_24

end subroutine weno_seven_reconstruction_0

!> Reconstruction in the third upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_1(q4, p1)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p1            !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]

  p1 = (14 * (q4(2) + q4(3)) - 2 * (q4(1) + q4(4))) * C1_24

end subroutine weno_seven_reconstruction_1

!> Reconstruction in the second upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_2(q4, p2)
  real, intent(in) :: q4(4)             !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p2             !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0   ! One twenty fourth [nondim]

  p2 = (((2 * q4(1) - 10 * q4(2)) + 26 * q4(3)) + 6 * q4(4)) * C1_24

end subroutine weno_seven_reconstruction_2

!> Reconstruction in the first upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_3(q4, p3)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p3            !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]

  p3 = (((-6 * q4(1) + 26 * q4(2)) - 46 * q4(3)) + 50 * q4(4)) * C1_24

end subroutine weno_seven_reconstruction_3

function CoriolisAdv_stencil(CS) result(stencil)
  type(CoriolisAdv_CS), intent(in)  :: CS  !< Control structure for MOM_CoriolisAdv
  integer :: stencil  !< The halo stencil size for the Coriolis advection scheme

  stencil = 2
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) stencil = 4
  if (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) stencil = 3

end function CoriolisAdv_stencil


!> Initializes the control structure for MOM_CoriolisAdv
subroutine CoriolisAdv_init(Time, G, GV, US, param_file, diag, AD, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(accel_diag_ptrs),   target, intent(inout) :: AD !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(inout) :: CS   !< Control structure for MOM_CoriolisAdv
  ! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_CoriolisAdv" ! This module's name.
  character(len=20)  :: tmpstr
  character(len=400) :: mesg
  logical :: use_weno
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
#ifdef __NVCOMPILER_OPENMP_GPU
  integer, parameter :: default_nkblock = 0
#else
  integer, parameter :: default_nkblock = 1
#endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CORIOLIS_ADV_NKBLOCK", CS%nkblock, &
                 "The k-direction block size used in Coriolis and momentum advection "//&
                 "calculations. The default 0 setting dynamically uses the full vertical column.", &
                 default=default_nkblock, layoutParam=.true.)
  if (CS%nkblock < 0) call MOM_error(FATAL, "CORIOLIS_ADV_NKBLOCK must be >= 0.")
  call get_param(param_file, mdl, "NOSLIP", CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise "//&
                 "free slip boundary conditions are assumed. The "//&
                 "implementation of the free slip BCs on a C-grid is much "//&
                 "cleaner than the no slip BCs. The use of free slip BCs "//&
                 "is strongly encouraged, and no slip BCs are not used with "//&
                 "the biharmonic viscosity.", default=.false.)

  call get_param(param_file, mdl, "CORIOLIS_EN_DIS", CS%Coriolis_En_Dis, &
                 "If true, two estimates of the thickness fluxes are used "//&
                 "to estimate the Coriolis term, and the one that "//&
                 "dissipates energy relative to the other one is used.", &
                 default=.false.)

  ! Set %Coriolis_Scheme
  ! (Select the baseline discretization for the Coriolis term)
  call get_param(param_file, mdl, "CORIOLIS_SCHEME", tmpstr, &
                 "CORIOLIS_SCHEME selects the discretization for the "//&
                 "Coriolis terms. Valid values are: \n"//&
                 "\t SADOURNY75_ENERGY - Sadourny, 1975; energy cons. \n"//&
                 "\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990 \n"//&
                 "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy cons. \n"//&
                 "\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; En. + Enst.\n"//&
                 "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with \n"//&
                 "\t                      Arakawa & Hsu and Sadourny energy \n"//&
                 "\t WENOVI5TH_PV_ENSTRO   - 5th-order WENO PV enstrophy \n"//&
                 "\t WENOVI3RD_PV_ENSTRO   - 3rd-order WENO PV enstrophy \n"//&
                 "\t WENOVI7TH_PV_ENSTRO  - 7th-order WENO PV enstrophy \n", &
                 default=SADOURNY75_ENERGY_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (SADOURNY75_ENERGY_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENERGY
    case (ARAKAWA_HSU_STRING)
      CS%Coriolis_Scheme = ARAKAWA_HSU90
    case (SADOURNY75_ENSTRO_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENSTRO
    case (ARAKAWA_LAMB_STRING)
      CS%Coriolis_Scheme = ARAKAWA_LAMB81
    case (AL_BLEND_STRING)
      CS%Coriolis_Scheme = AL_BLEND
    case (ROBUST_ENSTRO_STRING)
      CS%Coriolis_Scheme = ROBUST_ENSTRO
      CS%Coriolis_En_Dis = .false.
    case (WENOVI7TH_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi7th_PV_ENSTRO
    case (WENOVI5TH_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi5th_PV_ENSTRO
    case (WENOVI3RD_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi3rd_PV_ENSTRO
    case default
      call MOM_mesg('CoriolisAdv_init: Coriolis_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: Unrecognized setting "// &
            "#define CORIOLIS_SCHEME "//trim(tmpstr)//" found in input file.")
  end select

  use_weno = CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO &
      .or. CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO &
      .or. CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO

  if (use_weno) then
    call get_param(param_file, mdl, "WENO_VELOCITY_SMOOTH", CS%weno_velocity_smooth, &
            "If true, use velocity to compute weighting for WENO. ", &
                  default=.false.)
  endif

  if (CS%Coriolis_Scheme == AL_BLEND) then
    call get_param(param_file, mdl, "CORIOLIS_BLEND_WT_LIN", CS%wt_lin_blend, &
                 "A weighting value for the ratio of inverse thicknesses, "//&
                 "beyond which the blending between Sadourny Energy and "//&
                 "Arakawa & Hsu goes linearly to 0 when CORIOLIS_SCHEME "//&
                 "is ARAWAKA_LAMB_BLEND. This must be between 1 and 1e-16.", &
                 units="nondim", default=0.125)
    call get_param(param_file, mdl, "CORIOLIS_BLEND_F_EFF_MAX", CS%F_eff_max_blend, &
                 "The factor by which the maximum effective Coriolis "//&
                 "acceleration from any point can be increased when "//&
                 "blending different discretizations with the "//&
                 "ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be "//&
                 "greater than 2.0 (the max value for Sadourny energy).", &
                 units="nondim", default=4.0)
    CS%wt_lin_blend = min(1.0, max(CS%wt_lin_blend,1e-16))
    if (CS%F_eff_max_blend < 2.0) call MOM_error(WARNING, "CoriolisAdv_init: "//&
           "CORIOLIS_BLEND_F_EFF_MAX should be at least 2.")
  endif

  mesg = "If true, the Coriolis terms at u-points are bounded by "//&
         "the four estimates of (f+rv)v from the four neighboring "//&
         "v-points, and similarly at v-points."
  if (CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) then
    mesg = trim(mesg)//"  This option is "//&
                 "always effectively false with CORIOLIS_EN_DIS defined and "//&
                 "CORIOLIS_SCHEME set to "//trim(SADOURNY75_ENERGY_STRING)//"."
  else
    mesg = trim(mesg)//"  This option would "//&
                 "have no effect on the SADOURNY Coriolis scheme if it "//&
                 "were possible to use centered difference thickness fluxes."
  endif
  call get_param(param_file, mdl, "BOUND_CORIOLIS", CS%bound_Coriolis, mesg, &
                 default=.false.)
  if ((CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) .or. &
      (CS%Coriolis_Scheme == ROBUST_ENSTRO)) CS%bound_Coriolis = .false.

  ! Set KE_Scheme (selects discretization of KE)
  call get_param(param_file, mdl, "KE_SCHEME", tmpstr, &
                 "KE_SCHEME selects the discretization for acceleration "//&
                 "due to the kinetic energy gradient. Valid values are: \n"//&
                 "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV, KE_UP3", &
                 default=KE_ARAKAWA_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (KE_ARAKAWA_STRING); CS%KE_Scheme = KE_ARAKAWA
    case (KE_SIMPLE_GUDONOV_STRING); CS%KE_Scheme = KE_SIMPLE_GUDONOV
    case (KE_GUDONOV_STRING); CS%KE_Scheme = KE_GUDONOV
    case (KE_UP3_STRING); CS%KE_Scheme = KE_UP3
    case default
      call MOM_mesg('CoriolisAdv_init: KE_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
               "#define KE_SCHEME "//trim(tmpstr)//" in input file is invalid.")
  end select

  if (CS%KE_Scheme == KE_UP3) then
    call get_param(param_file, mdl, "KE_USE_LIMITER", CS%KE_use_limiter, &
            "If true, use Koren limiter for KE_UP3 scheme", &
                  default=.True.)
  endif

  ! Set PV_Adv_Scheme (selects discretization of PV advection)
  call get_param(param_file, mdl, "PV_ADV_SCHEME", tmpstr, &
                 "PV_ADV_SCHEME selects the discretization for PV "//&
                 "advection. Valid values are: \n"//&
                 "\t PV_ADV_CENTERED - centered (aka Sadourny, 75) \n"//&
                 "\t PV_ADV_UPWIND1  - upwind, first order", &
                 default=PV_ADV_CENTERED_STRING)
  select case (uppercase(tmpstr))
    case (PV_ADV_CENTERED_STRING)
      CS%PV_Adv_Scheme = PV_ADV_CENTERED
    case (PV_ADV_UPWIND1_STRING)
      CS%PV_Adv_Scheme = PV_ADV_UPWIND1
    case default
      call MOM_mesg('CoriolisAdv_init: PV_Adv_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
                     "#DEFINE PV_ADV_SCHEME in input file is invalid.")
  end select

  CS%id_rv = register_diag_field('ocean_model', 'RV', diag%axesBL, Time, &
     'Relative Vorticity', 's-1', conversion=US%s_to_T)

  CS%id_PV = register_diag_field('ocean_model', 'PV', diag%axesBL, Time, &
     'Potential Vorticity', 'm-1 s-1', conversion=GV%m_to_H*US%s_to_T)

  CS%id_gKEu = register_diag_field('ocean_model', 'gKEu', diag%axesCuL, Time, &
     'Zonal Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_gKEv = register_diag_field('ocean_model', 'gKEv', diag%axesCvL, Time, &
     'Meridional Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxu = register_diag_field('ocean_model', 'rvxu', diag%axesCvL, Time, &
     'Meridional Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxv = register_diag_field('ocean_model', 'rvxv', diag%axesCuL, Time, &
     'Zonal Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_CAuS = register_diag_field('ocean_model', 'CAu_Stokes', diag%axesCuL, Time, &
     'Zonal Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  CS%id_CAvS = register_diag_field('ocean_model', 'CAv_Stokes', diag%axesCvL, Time, &
     'Meridional Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  !CS%id_hf_gKEu = register_diag_field('ocean_model', 'hf_gKEu', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEu_2d = register_diag_field('ocean_model', 'hf_gKEu_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_gKEv = register_diag_field('ocean_model', 'hf_gKEv', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEv_2d = register_diag_field('ocean_model', 'hf_gKEv_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_gKEu = register_diag_field('ocean_model', 'h_gKEu', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEu_2d = register_diag_field('ocean_model', 'intz_gKEu_2d', diag%axesCu1, Time, &
     'Depth-integral of Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_gKEv = register_diag_field('ocean_model', 'h_gKEv', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEv_2d = register_diag_field('ocean_model', 'intz_gKEv_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  !CS%id_hf_rvxu = register_diag_field('ocean_model', 'hf_rvxu', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxu_2d = register_diag_field('ocean_model', 'hf_rvxu_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_rvxv = register_diag_field('ocean_model', 'hf_rvxv', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxv_2d = register_diag_field('ocean_model', 'hf_rvxv_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_rvxu = register_diag_field('ocean_model', 'h_rvxu', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxu_2d = register_diag_field('ocean_model', 'intz_rvxu_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_rvxv = register_diag_field('ocean_model', 'h_rvxv', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxv_2d = register_diag_field('ocean_model', 'intz_rvxv_2d', diag%axesCu1, Time, &
     'Depth-integral of Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  ! Allocate memory needed for the diagnostics that have been enabled.
  if ((CS%id_gKEu > 0) .or. (CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. &
      (CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_gKEv > 0) .or. (CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. &
      (CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEv, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxu > 0) .or. (CS%id_hf_rvxu_2d > 0) .or. &
    ! (CS%id_hf_rvxu > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_u, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxv > 0) .or. (CS%id_hf_rvxv_2d > 0) .or. &
    ! (CS%id_hf_rvxv > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_v, IsdB, IedB, jsd, jed, nz)
  endif

  if ((CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. (CS%id_hf_rvxu > 0) .or. &
      (CS%id_hf_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_v, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. (CS%id_hf_rvxv > 0) .or. &
      (CS%id_hf_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_u, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hv, isd, ied, JsdB, JedB, nz)
  endif

end subroutine CoriolisAdv_init

!> Destructor for coriolisadv_cs
subroutine CoriolisAdv_end(CS)
  type(CoriolisAdv_CS), intent(inout) :: CS !< Control structure for MOM_CoriolisAdv
end subroutine CoriolisAdv_end

!> \namespace mom_coriolisadv
!!
!! This file contains the subroutine that calculates the time
!! derivatives of the velocities due to Coriolis acceleration and
!! momentum advection.  This subroutine uses either a vorticity
!! advection scheme from Arakawa and Hsu, Mon. Wea. Rev. 1990, or
!! Sadourny's (JAS 1975) energy conserving scheme.  Both have been
!! modified to use general orthogonal coordinates as described in
!! Arakawa and Lamb, Mon. Wea. Rev. 1981.  Both schemes are second
!! order accurate, and allow for vanishingly small layer thicknesses.
!! The Arakawa and Hsu scheme globally conserves both total energy
!! and potential enstrophy in the limit of nondivergent flow.
!! Sadourny's energy conserving scheme conserves energy if the flow
!! is nondivergent or centered difference thickness fluxes are used.
!!
!! A small fragment of the grid is shown below:
!! \verbatim
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, CAv, vh
!!    j    x ^ x ^ x   At >:  u, CAu, uh, a, b, c, d
!!    j    > o > o >   At o:  h, KE
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!! The boundaries always run through q grid points (x).

end module MOM_CoriolisAdv
