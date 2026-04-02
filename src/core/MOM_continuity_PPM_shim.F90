module MOM_continuity_PPM_shim
  use iso_c_binding, only : c_double, c_int
  use MOM_grid, only : ocean_grid_type
  use TIM_helperF, only : TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN, &
          TIMH_CAPTURE_INPUT, TIMH_CAPTURE_OUTPUT, TIMH_RUN, get_mode_env
  implicit none ; private

#include <MOM_memory.h>

  logical, parameter :: use_AMREX = .TRUE. 
  !----------------------------------------
  ! C interface (bridge to C++)
  !----------------------------------------
  interface
    subroutine ppm_limit_pos_c(h_in, h_L, h_R, h_min,  &
                               lo_i, hi_i, lo_j, hi_j, &
			       i_min, i_max, j_min, j_max, mode) bind(C)
      use iso_c_binding
      implicit none

      real(c_double), intent(in)    :: h_in(*)
      real(c_double), intent(inout) :: h_L(*)
      real(c_double), intent(inout) :: h_R(*)

      real(c_double), intent(in) :: h_min

      integer(c_int), intent(in) :: lo_i, hi_i
      integer(c_int), intent(in) :: lo_j, hi_j
      integer(c_int), intent(in) :: i_min, i_max
      integer(c_int), intent(in) :: j_min, j_max
      integer(c_int), intent(in) :: mode
	
    end subroutine ppm_limit_pos_c
  end interface
  interface
    subroutine ppm_limit_cw84_c(h_in, h_L, h_R,  &
                               lo_i, hi_i, lo_j, hi_j, &
			       i_min, i_max, j_min, j_max, mode) bind(C)
      use iso_c_binding
      implicit none

      real(c_double), intent(in)    :: h_in(*)
      real(c_double), intent(inout) :: h_L(*)
      real(c_double), intent(inout) :: h_R(*)

      integer(c_int), intent(in) :: lo_i, hi_i
      integer(c_int), intent(in) :: lo_j, hi_j
      integer(c_int), intent(in) :: i_min, i_max
      integer(c_int), intent(in) :: j_min, j_max
      integer(c_int), intent(in) :: mode
	
    end subroutine ppm_limit_cw84_c
  end interface

  public :: PPM_limit_pos_shim
  public :: PPM_limit_cw84_shim

contains

  !----------------------------------------
  ! Drop-in replacement for PPM_limit_pos
  !----------------------------------------
  subroutine PPM_limit_pos_shim(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
    implicit none

    type(ocean_grid_type), intent(in) :: G
    real(c_double), dimension(SZI_(G), SZJ_(G)), intent(in)    :: h_in
    real(c_double), dimension(SZI_(G), SZJ_(G)), intent(inout) :: h_L
    real(c_double), dimension(SZI_(G), SZJ_(G)), intent(inout) :: h_R
    real(c_double), intent(in)    :: h_min
    integer, intent(in) :: iis, iie, jis, jie

    ! local variables
    integer :: imin, imax, jmin, jmax
    integer :: ppm_limit_pos_mode

    imin = LBOUND(h_in,dim=1)
    imax = UBOUND(h_in,dim=1)
    jmin = LBOUND(h_in,dim=2)
    jmax = UBOUND(h_in,dim=2)

    ppm_limit_pos_mode = get_mode_env("PPM_LIMIT_POS_MODE",default=TIMH_runFORTRAN)

    ! Call C++ bridge
    select case (ppm_limit_pos_mode)
       case (TIMH_runFORTRAN)

          ! Run Fortran code
          call ppm_limit_pos(h_in, h_L, h_R, h_min,  &
             G, iis, iie, jis, jie) 

       case (TIMH_capture)

           ! capture the input state
           call ppm_limit_pos_c(h_in, h_L, h_R, h_min,  &
              iis, iie, jis, jie, imin, imax, jmin, jmax, TIMH_CAPTURE_INPUT)

          ! Run Fortran truth 
          call ppm_limit_pos(h_in, h_L, h_R, h_min,  &
             G, iis, iie, jis, jie) 

          ! capture the output state
           call ppm_limit_pos_c(h_in, h_L, h_R, h_min,  &
              iis, iie, jis, jie, imin, imax, jmin, jmax, TIMH_CAPTURE_OUTPUT)

       case (TIMH_runAMREX)

          ! Run AMReX code
          call ppm_limit_pos_c(h_in, h_L, h_R, h_min,  &
               iis, iie, jis, jie, imin, imax, jmin, jmax, TIMH_RUN)

    end select

end subroutine PPM_limit_pos_shim

  !----------------------------------------

  ! Drop-in replacement for PPM_limit_cw84
  !----------------------------------------
  subroutine PPM_limit_cw84_shim(h_in, h_L, h_R, G, iis, iie, jis, jie)
    implicit none

    type(ocean_grid_type), intent(in) :: G
    real(c_double), dimension(SZI_(G), SZJ_(G)), intent(in)    :: h_in
    real(c_double), dimension(SZI_(G), SZJ_(G)), intent(inout) :: h_L
    real(c_double), dimension(SZI_(G), SZJ_(G)), intent(inout) :: h_R
    integer, intent(in) :: iis, iie, jis, jie

    ! local variables
    integer :: imin, imax, jmin, jmax
    integer :: mode

    imin = LBOUND(h_in,dim=1)
    imax = UBOUND(h_in,dim=1)
    jmin = LBOUND(h_in,dim=2)
    jmax = UBOUND(h_in,dim=2)

    mode = get_mode_env("PPM_LIMIT_CW84_MODE", default=TIMH_runFORTRAN)
    ! Call C++ bridge
    select case (mode)
       case (TIMH_runFORTRAN)

          ! Run Fortran code
          call ppm_limit_cw84(h_in, h_L, h_R, &
              G, iis, iie, jis, jie) 

       case (TIMH_capture)

          ! capture the input state
          call ppm_limit_cw84_c(h_in, h_L, h_R,  &
             iis, iie, jis, jie, imin, imax, jmin, jmax, TIMH_CAPTURE_INPUT)

          ! Run Fortran truth
          call ppm_limit_cw84(h_in, h_L, h_R, &
              G, iis, iie, jis, jie) 

          ! capture the input state
          call ppm_limit_cw84_c(h_in, h_L, h_R,  &
             iis, iie, jis, jie, imin, imax, jmin, jmax, TIMH_CAPTURE_OUTPUT)

       case (TIMH_runAMREX)

          ! Run AMReX code
          call ppm_limit_cw84_c(h_in, h_L, h_R,  &
             iis, iie, jis, jie, imin, imax, jmin, jmax, TIMH_RUN)

     end select

end subroutine PPM_limit_cw84_shim

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant thickness if the mean thickness is less
!! than h_min, with a minimum of h_min otherwise.
subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_L !< Left thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_R !< Right thickness in the reconstruction [H ~> m or kg m-2].
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  integer,                           intent(in)  :: iis      !< Start of i index range.
  integer,                           intent(in)  :: iie      !< End of i index range.
  integer,                           intent(in)  :: jis      !< Start of j index range.
  integer,                           intent(in)  :: jie      !< End of j index range.

! Local variables
  real    :: curv  ! The grid-normalized curvature of the three thicknesses  [H ~> m or kg m-2]
  real    :: dh    ! The difference between the edge thicknesses             [H ~> m or kg m-2]
  real    :: scale ! A scaling factor to reduce the curvature of the fit               [nondim]
  integer :: i,j

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*((h_L(i,j) + h_R(i,j)) - 2.0*h_in(i,j))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j) - h_L(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j) <= h_min) then
          h_L(i,j) = h_in(i,j) ; h_R(i,j) = h_in(i,j)
        elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j) = h_in(i,j) + scale*(h_L(i,j) - h_in(i,j))
          h_R(i,j) = h_in(i,j) + scale*(h_R(i,j) - h_in(i,j))
        endif
      endif
    endif
  enddo ; enddo

end subroutine PPM_limit_pos

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! according to the monotonic prescription of Colella and Woodward, 1984.
subroutine PPM_limit_CW84(h_in, h_L, h_R, G, iis, iie, jis, jie)
  type(ocean_grid_type),             intent(in)  :: G     !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in  !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_L !< Left thickness in the reconstruction,
                                                          !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h_R !< Right thickness in the reconstruction,
                                                          !! [H ~> m or kg m-2].
  integer,                           intent(in)  :: iis   !< Start of i index range.
  integer,                           intent(in)  :: iie   !< End of i index range.
  integer,                           intent(in)  :: jis   !< Start of j index range.
  integer,                           intent(in)  :: jie   !< End of j index range.

  ! Local variables
  real    :: h_i      ! A copy of the cell-average layer thickness                [H ~> m or kg m-2]
  real    :: RLdiff   ! The difference between the input edge values              [H ~> m or kg m-2]
  real    :: RLdiff2  ! The squared difference between the input edge values   [H2 ~> m2 or kg2 m-4]
  real    :: RLmean   ! The average of the input edge thicknesses                 [H ~> m or kg m-2]
  real    :: FunFac   ! A curious product of the thickness slope and curvature [H2 ~> m2 or kg2 m-4]
  integer :: i, j

  do j=jis,jie ; do i=iis,iie
    ! This limiter monotonizes the parabola following
    ! Colella and Woodward, 1984, Eq. 1.10
    h_i = h_in(i,j)
    if ( ( h_R(i,j) - h_i ) * ( h_i - h_L(i,j) ) <= 0. ) then
      h_L(i,j) = h_i ; h_R(i,j) = h_i
    else
      RLdiff = h_R(i,j) - h_L(i,j)            ! Difference of edge values
      RLmean = 0.5 * ( h_R(i,j) + h_L(i,j) )  ! Mean of edge values
      FunFac = 6. * RLdiff * ( h_i - RLmean ) ! Some funny factor
      RLdiff2 = RLdiff * RLdiff               ! Square of difference
      if ( FunFac >  RLdiff2 ) h_L(i,j) = 3. * h_i - 2. * h_R(i,j)
      if ( FunFac < -RLdiff2 ) h_R(i,j) = 3. * h_i - 2. * h_L(i,j)
    endif
  enddo ; enddo

  return
end subroutine PPM_limit_CW84


end module MOM_continuity_PPM_shim
