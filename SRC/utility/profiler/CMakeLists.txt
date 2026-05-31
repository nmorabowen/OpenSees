#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#
# Ladruno Stack Profiler — Phase P1 build wiring.
#
# SINGLE-BINARY model (06_profiler.md, revised 2026-05-30):
#
#   * The WHOLE profiler (coarse + deep + memory) is compiled into the one shipped
#     binary and gated by RUNTIME flags (enabled/deep/mem). No special build is
#     needed for deep/memory — the user toggles them per run. Cost when off is one
#     predicted branch. PerfClock.cpp + Profiler.cpp are HDF5-free with no OpenSees
#     dependency, so they link unconditionally into OPS_Utilities.
#
#   * The HDF5 writer (ProfilerHDF5Writer.cpp) is the ONLY TU that #include
#     <hdf5.h>. It is compiled only when HDF5 >= 1.12.0 was found — mirroring the
#     exact HDF5_FOUND / HDF5_VERSION idiom SRC/recorder/CMakeLists.txt uses for
#     MPCORecorder.cpp. The top-level CMakeLists.txt already does
#     find_package(HDF5 REQUIRED), include_directories(${HDF5_INCLUDE_DIR}), and
#     links HDF5::HDF5 into the final OpenSees target, so this TU needs no extra
#     link wiring here (OPS_Utilities is an OBJECT library — its objects inherit
#     the HDF5 include path globally and the symbols resolve at the final link).
#
#   * OPS_PROFILE_DISABLE is an OPTIONAL opt-OUT (default OFF = profiler present).
#     When ON it adds -DOPS_PROFILE_DISABLE and every instrumentation macro becomes
#     literally nothing — a benchmark/purist build with not even the runtime branch.
#     This is the ONLY build switch; there is no opt-in compile flag anymore.
#
# All sources register onto the existing OPS_Utilities OBJECT library (defined at
# the top-level CMakeLists.txt), exactly like SRC/utility/CMakeLists.txt registers
# Timer.cpp / SimulationInformation.cpp.
#==============================================================================

# ---- coarse timing core: ALWAYS compiled, HDF5-free ----
target_sources(OPS_Utilities
    PRIVATE
        PerfClock.cpp
        Profiler.cpp
    PUBLIC
        PerfClock.h
        Profiler.h
        ProfilerHDF5Writer.h
        ProfilerMacros.h
)

# ---- HDF5 writer: compiled only when HDF5 >= 1.12.0 is available ----
# Same guard idiom as SRC/recorder/CMakeLists.txt (MPCORecorder.cpp).
if(HDF5_FOUND)
    if(HDF5_VERSION VERSION_GREATER_EQUAL 1.12.0)
        target_sources(OPS_Utilities
            PRIVATE
                ProfilerHDF5Writer.cpp
        )
        # HDF5::HDF5 is linked into the final OpenSees target at the top level;
        # OPS_Utilities is an OBJECT library so we only need the include path,
        # which include_directories(${HDF5_INCLUDE_DIR}) provides globally. The
        # explicit link below is harmless and future-proofs a standalone build.
        target_link_libraries(OPS_Utilities PUBLIC HDF5::HDF5)
        message(STATUS "OPS >>> Profiler: HDF5 writer enabled (HDF5 ${HDF5_VERSION})")
    else()
        message(STATUS "OPS >>> Profiler: HDF5 < 1.12.0 — coarse timing only (no writer)")
    endif()
else()
    message(STATUS "OPS >>> Profiler: no HDF5 — coarse timing only (no writer)")
endif()

# ---- OPS_PROFILE_DISABLE: optional opt-OUT, compile the profiler fully out ----
# Default OFF = the profiler is present and runtime-gated (the normal shipped
# binary). When ON, every OPS_PROFILE_* macro becomes ((void)0) everywhere it is
# used — so the define must be GLOBAL: any engine TU that includes ProfilerMacros.h
# (P2+ seam sites) must see it too. add_compile_definitions here applies to this
# directory and below; for a whole-tree compile-out promote it to the top-level
# CMakeLists. This is a rare benchmark/purist switch, not the normal path.
option(OPS_PROFILE_DISABLE "Compile the stack profiler entirely out (no runtime branch)" OFF)
if(OPS_PROFILE_DISABLE)
    add_compile_definitions(OPS_PROFILE_DISABLE)
    message(STATUS "OPS >>> Profiler: DISABLED at compile time (OPS_PROFILE_DISABLE)")
else()
    message(STATUS "OPS >>> Profiler: present, runtime-gated (coarse/deep/memory)")
endif()

target_include_directories(OPS_Utilities PUBLIC ${CMAKE_CURRENT_LIST_DIR})
