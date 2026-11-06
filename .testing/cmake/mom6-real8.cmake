# This file is part of MOM6, the Modular Ocean Model version 6.
# See the LICENSE file for licensing information.
# SPDX-License-Identifier: Apache-2.0

# Included by mom6-variant-flags.cmake. Finds the compiler flag that makes
# default reals 8 bytes, as ac/m4/ax_fc_real8.m4 does for the autoconf build:
# the first candidate under which both a default real and a double precision
# variable have kind selected_real_kind(15, 307) wins. The flag is added to
# CMAKE_Fortran_FLAGS, so that every target, including the bundled packages, is
# compiled with it, as autoconf adds it to FCFLAGS.

include(CheckFortranSourceCompiles)

if(NOT DEFINED MOM6_REAL8_FCFLAGS)
  # A zero stride (m or n = 0) is a compile-time error
  set(_mom6_real8_source [=[
      program real8
      real :: x(4)
      double precision :: y(4)
      integer, parameter :: &
        m = merge(1, 0, kind(x(1)) == selected_real_kind(15, 307)), &
        n = merge(1, 0, kind(y(1)) == selected_real_kind(15, 307))
      print *, x(::m)
      print *, y(::n)
      end program real8
  ]=])

  set(_mom6_real8_flag "unsupported")
  set(_i 0)
  foreach(_flag IN ITEMS
      none
      "-fdefault-real-8"
      "-fdefault-real-8 -fdefault-double-8"
      "-r8"
      "-real-kind 64"
      "-Mr8"
      "-s real64"
      "-double")
    if(_flag STREQUAL "none")
      set(CMAKE_REQUIRED_FLAGS "")
    else()
      set(CMAKE_REQUIRED_FLAGS "${_flag}")
    endif()
    check_fortran_source_compiles("${_mom6_real8_source}" MOM6_REAL8_FLAG_${_i}
      SRC_EXT F90)
    unset(CMAKE_REQUIRED_FLAGS)
    if(MOM6_REAL8_FLAG_${_i})
      set(_mom6_real8_flag "${_flag}")
      break()
    endif()
    math(EXPR _i "${_i} + 1")
  endforeach()

  if(_mom6_real8_flag STREQUAL "unsupported")
    message(WARNING "No Fortran flag found to force 8-byte reals")
    set(_mom6_real8_flag "")
  elseif(_mom6_real8_flag STREQUAL "none")
    set(_mom6_real8_flag "")
  endif()
  set(MOM6_REAL8_FCFLAGS "${_mom6_real8_flag}" CACHE STRING
    "Fortran flags that force 8-byte default reals")
  unset(_mom6_real8_flag)
  unset(_mom6_real8_source)
  unset(_i)
endif()

message(STATUS "Fortran flags for 8-byte reals: '${MOM6_REAL8_FCFLAGS}'")
string(APPEND CMAKE_Fortran_FLAGS " ${MOM6_REAL8_FCFLAGS}")
