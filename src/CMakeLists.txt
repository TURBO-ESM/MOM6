set(_EXT    "${MOM6_SOURCE_DIR}/config_src/external")

# mom6_ocean: the MOM6 model — every src/ file except the model-agnostic
# framework (mom6_framework) and the infra wrapper (mom6_infra), plus the
# config_src/external/ stubs (MARBL excluded — provided by MARBL::marbl).
# One library; CMake's Fortran scanner orders compilation from the module graph.
add_library(mom6_ocean STATIC)
add_library(MOM6::ocean ALIAS mom6_ocean)

set_target_properties(mom6_ocean PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(mom6_ocean
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
    PRIVATE   ${MOM6_MEMORY_DIR}
)
target_link_libraries(mom6_ocean
    PUBLIC MOM6::framework MOM6::CVMix MOM6::GSW MARBL::marbl
)

# _TIM gates the AMReX bridge branches in core/MOM_continuity_PPM.F90 -- the only
# file in the tree with `#ifdef _TIM` guards. PRIVATE: it changes how this target
# compiles and is not part of its interface.
#
# Keyed on the backend rather than on an opt-in flag. The guarded branches exist
# only on this branch -- dev/turbo, dev/ncar and the rest carry none -- so a build
# option would only ever be turned on, by everyone, everywhere this file is used.
# The backend is the one thing that genuinely varies: FMS2 cannot provide
# turbotmp_*_bridge, so under FMS2 the branches must stay compiled out.
#
# The bridge is protoMOMxx::mom_bridge, from an installed protoMOMxx (built
# against the same TIM). It is found and linked here, on the target whose sources
# call it, rather than in the top-level dependency block or onto the infra
# interface library: when turbo-stack builds MOM6 with add_subdirectory() it
# defines TURBO::infra_r8 itself, and MOM6's whole standalone dependency block is
# skipped. Doing it here works either way. The kernels (protoMOMxx::mom_kernels)
# come in through the bridge's link interface and are not linked directly.
#
# MOM_BRIDGE is required, not optional: this branch's
# src/core/MOM_continuity_PPM.F90 calls the bridge unconditionally under TIM.
# Asking for it as a component means a missing bridge fails here, with the
# component named, instead of at link time as undefined turbotmp_*_bridge
# references.
if(MOM6_INFRA STREQUAL "TIM")
    if(NOT TARGET protoMOMxx::mom_bridge)
        find_package(protoMOMxx REQUIRED COMPONENTS MOM_BRIDGE)
    endif()
    if(NOT TARGET protoMOMxx::mom_bridge)
        message(FATAL_ERROR
            "MOM6_INFRA=TIM, but the protoMOMxx in use does not provide "
            "protoMOMxx::mom_bridge. Build and install protoMOMxx with "
            "-DPROTOMOM_INSTALL=ON (and PROTOMOM_MOM_BRIDGE=ON) and put its "
            "install prefix on CMAKE_PREFIX_PATH.")
    endif()
    target_compile_definitions(mom6_ocean PRIVATE _TIM)
    target_link_libraries(mom6_ocean PRIVATE protoMOMxx::mom_bridge)
endif()

# I/O modules from src/framework/: they depend on grid (which also lives in
# mom6_ocean), so they compile here rather than in mom6_framework.
target_sources(mom6_ocean PRIVATE
    framework/MOM_io.F90
    framework/MOM_get_input.F90
    framework/MOM_write_cputime.F90
    framework/MOM_interpolate.F90
    framework/MOM_restart.F90
)

# Other src/framework/ files that depend on model modules (EOS, regridding,
# diagnostics, core) and so compile here rather than in mom6_framework.
target_sources(mom6_ocean PRIVATE
    framework/MOM_ANN.F90
    framework/MOM_diag_buffers.F90
    framework/MOM_diag_mediator.F90
    framework/MOM_diag_remap.F90
    framework/MOM_horizontal_regridding.F90
)

# External stubs (MARBL excluded — provided by MARBL::marbl).
target_sources(mom6_ocean PRIVATE
    ${_EXT}/GFDL_ocean_BGC/MOM_generic_tracer.F90
    ${_EXT}/ODA_hooks/kdtree.f90
    ${_EXT}/ODA_hooks/ocean_da_core.F90
    ${_EXT}/ODA_hooks/ocean_da_types.F90
    ${_EXT}/ODA_hooks/write_ocean_obs.F90
    ${_EXT}/database_comms/MOM_database_comms.F90
    ${_EXT}/database_comms/database_client_interface.F90
    ${_EXT}/drifters/MOM_particles.F90
    ${_EXT}/drifters/MOM_particles_types.F90
    ${_EXT}/stochastic_physics/get_stochy_pattern.F90
    ${_EXT}/stochastic_physics/stochastic_physics.F90
)

# Compilation order between these subdirectories no longer matters: they all
# add sources to the single mom6_ocean target and CMake orders by module graph.
add_subdirectory(core)
add_subdirectory(framework)
add_subdirectory(ALE)
add_subdirectory(diagnostics)
add_subdirectory(equation_of_state)
add_subdirectory(ice_shelf)
add_subdirectory(initialization)
add_subdirectory(ocean_data_assim)
add_subdirectory(parameterizations)
add_subdirectory(tracer)
add_subdirectory(user)
