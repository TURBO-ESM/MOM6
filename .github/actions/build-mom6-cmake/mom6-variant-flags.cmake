# This file is part of MOM6, the Modular Ocean Model version 6.
# See the LICENSE file for licensing information.
# SPDX-License-Identifier: Apache-2.0

# Loaded by build-mom6-cmake as CMAKE_PROJECT_MOM6_INCLUDE, at the end of
# project(MOM6). Adds MOM6_VARIANT_FCFLAGS (the .testing debug or repro flags)
# to MOM6's own targets only. The bundled packages keep the global flags, as in
# the .testing build, which compiles them with FCFLAGS_FMS; the debug flags
# include -std=f2018, which they do not conform to.

set(_mom6_bundled_targets mom6_cvmix mom6_gsw mom6_marbl)

# Collect the build-system targets defined in a directory and its subdirectories
function(_mom6_collect_targets dir out)
  get_property(_targets DIRECTORY "${dir}" PROPERTY BUILDSYSTEM_TARGETS)
  get_property(_subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
  foreach(_subdir IN LISTS _subdirs)
    _mom6_collect_targets("${_subdir}" _sub_targets)
    list(APPEND _targets ${_sub_targets})
  endforeach()
  set(${out} ${_targets} PARENT_SCOPE)
endfunction()

function(_mom6_add_variant_flags)
  separate_arguments(_flags UNIX_COMMAND "${MOM6_VARIANT_FCFLAGS}")
  _mom6_collect_targets("${MOM6_SOURCE_DIR}" _targets)
  foreach(_target IN LISTS _targets)
    get_target_property(_type ${_target} TYPE)
    if(_target IN_LIST _mom6_bundled_targets OR _type STREQUAL "INTERFACE_LIBRARY")
      continue()
    endif()
    target_compile_options(${_target} PRIVATE ${_flags})
    message(STATUS "MOM6 variant flags on ${_target}: ${MOM6_VARIANT_FCFLAGS}")
  endforeach()
endfunction()

# Run once every subdirectory has defined its targets
cmake_language(DEFER DIRECTORY "${MOM6_SOURCE_DIR}" CALL _mom6_add_variant_flags)
