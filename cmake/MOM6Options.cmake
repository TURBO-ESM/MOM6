set(MOM6_INFRA "FMS2" CACHE STRING "Infrastructure backend: FMS2 or TIM")
set_property(CACHE MOM6_INFRA PROPERTY STRINGS FMS2 TIM)

set(MOM6_MEMORY_MODE "dynamic_symmetric" CACHE STRING
    "MOM6 memory layout: dynamic_symmetric or dynamic_nonsymmetric")
set_property(CACHE MOM6_MEMORY_MODE PROPERTY STRINGS
    dynamic_symmetric dynamic_nonsymmetric)

# Compiles the `#ifdef _TIM` branches in src/core/MOM_continuity_PPM.F90, which
# call TIM's AMReX bridge kernels. Off by default: the bridge is an opt-in build
# of TIM (TIM_ENABLE_MOM_BRIDGE), and defining _TIM against a TIM built without it
# produces six undefined references at link time.
#
# Compiling the branches in does not execute them -- each is one arm of a runtime
# dispatch that defaults to TIMH_runFORTRAN, selected per kernel by environment
# variables (ZONAL_EDGE_THICKNESS_MODE and friends). What this buys is that the
# bridge is type-checked against the C++ signatures and linked on every build.
option(MOM6_ENABLE_TIM_BRIDGE
    "Compile the TIM AMReX bridge call sites (-D_TIM); requires MOM6_INFRA=TIM" OFF)

if(NOT MOM6_INFRA MATCHES "^(FMS2|TIM)$")
    message(FATAL_ERROR
        "Unknown MOM6_INFRA='${MOM6_INFRA}'. Valid values: FMS2, TIM")
endif()

# Rejected rather than ignored: the guarded branches call turbotmp_*_bridge, which
# only TIM provides, so with FMS2 this would compile and then fail at link. A
# caller that sets it on the wrong backend has made a mistake worth reporting.
if(MOM6_ENABLE_TIM_BRIDGE AND NOT MOM6_INFRA STREQUAL "TIM")
    message(FATAL_ERROR
        "MOM6_ENABLE_TIM_BRIDGE=ON requires MOM6_INFRA=TIM (got '${MOM6_INFRA}')")
endif()

if(NOT MOM6_MEMORY_MODE MATCHES "^(dynamic_symmetric|dynamic_nonsymmetric)$")
    message(FATAL_ERROR
        "Unknown MOM6_MEMORY_MODE='${MOM6_MEMORY_MODE}'. "
        "Valid values: dynamic_symmetric, dynamic_nonsymmetric")
endif()

# Resolved path to the selected memory-mode sources. Interpreting the memory
# mode once here lets every target just reference ${MOM6_MEMORY_DIR}.
set(MOM6_MEMORY_DIR "${MOM6_SOURCE_DIR}/config_src/memory/${MOM6_MEMORY_MODE}")

