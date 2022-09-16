# find libvdwxc if in non-standard location set environment variabled
# `VDWXC_ROOT` to the root directory

include(FindPackageHandleStandardArgs)
include(CheckSymbolExists)
include(cp2k_utils)

find_package(PkgConfig REQUIRED)

cp2k_set_default_paths(LIBVDWXC)

pkg_search_module(CP2K_LIBVDWXC libvdwxc>=${LibVDWXC_FIND_VERSION})

if (NOT CP2K_LIBVDWXC_FOUND)
  cp2k_find_libraries(LIBVDWXC "vdwxc")
endif()

if (NOT DEFINED CP2K_LIBVDWXC_INCLUDE_DIRS)
  cp2k_include_dirs(LIBVDWXC "vdwxc.h;vdwxc_mpi.h")
endif()

# try linking in C (C++ fails because vdwxc_mpi.h includes mpi.h inside extern
# "C"{...})
set(CMAKE_REQUIRED_LIBRARIES "${CP2K_LIBVDWXC_LINK_LIBRARIES}")
check_symbol_exists(vdwxc_init_mpi "${CP2K_LIBVDWXC_INCLUDE_DIRS}/vdwxc_mpi.h"
                    HAVE_LIBVDW_WITH_MPI)

find_package_handle_standard_args(LibVDWXC DEFAULT_MSG CP2K_LIBVDWXC_LINK_LIBRARIES
                                  CP2K_LIBVDWXC_INCLUDE_DIRS)

if(LibVDWXC_FOUND AND NOT TARGET CP2K_libvdwxc::libvdwxc)
  add_library(CP2K_libvdwxc::libvdwxc INTERFACE IMPORTED)
  set_target_properties(
    CP2K_libvdwxc::libvdwxc
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBVDWXC_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${CP2K_LIBVDWXC_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_LIBVDWXC_INCLUDE_DIRS CP2K_LIBVDWXC_LINK_LIBRARIES CP2K_LIBVDWXC_FOUND)
