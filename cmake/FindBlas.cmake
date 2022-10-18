# Copyright (c) 2022- ETH Zurich
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if(NOT
   (CMAKE_C_COMPILER_LOADED
    OR CMAKE_CXX_COMPILER_LOADED
    OR CMAKE_Fortran_COMPILER_LOADED))
  message(FATAL_ERROR "FindBLAS requires Fortran, C, or C++ to be enabled.")
endif()

set(CP2K_BLAS_VENDOR_LIST
    "auto"
    "MKL"
    "OpenBLAS"
    "SCI"
    "GenericBLAS"
    "Armpl"
    "FlexiBLAS"
    "Atlas")

set(__BLAS_VENDOR_LIST
    "MKL"
    "OpenBLAS"
    "SCI"
    "GenericBLAS"
    "Armpl"
    "FlexiBLAS"
    "Atlas")

set(CP2K_BLAS_VENDOR
    "auto"
    CACHE STRING "Blas library for computations on host")
set_property(CACHE CP2K_BLAS_VENDOR PROPERTY STRINGS ${CP2K_BLAS_VENDOR_LIST})

if(NOT ${CP2K_BLAS_VENDOR} IN_LIST CP2K_BLAS_VENDOR_LIST)
  message(FATAL_ERROR "Invalid Host BLAS backend")
endif()

set(CP2K_BLAS_THREAD_LIST "sequential" "thread" "gnu-thread" "intel-thread"
                          "tbb-thread" "openmp")

set(CP2K_BLAS_THREADING
    "sequential"
    CACHE STRING "threaded blas library")
set_property(CACHE CP2K_BLAS_THREADING PROPERTY STRINGS
                                                ${CP2K_BLAS_THREAD_LIST})

if(NOT ${CP2K_BLAS_THREADING} IN_LIST CP2K_BLAS_THREAD_LIST)
  message(FATAL_ERROR "Invalid threaded BLAS backend")
endif()

set(CP2K_BLAS_INTERFACE_BITS_LIST "32bits" "64bits")
set(CP2K_BLAS_INTERFACE
    "32bits"
    CACHE STRING
          "32 bits integers are used for indices, matrices and vectors sizes")
set_property(CACHE CP2K_BLAS_INTERFACE
             PROPERTY STRINGS ${CP2K_BLAS_INTERFACE_BITS_LIST})

if(NOT ${CP2K_BLAS_INTERFACE} IN_LIST CP2K_BLAS_INTERFACE_BITS_LIST)
  message(
    FATAL_ERROR
      "Invalid parameters. Blas and lapack can exist in two flavors 32 or 64 bits interfaces (relevant mostly for mkl)"
  )
endif()

set(CP2K_BLAS_FOUND FALSE)

# first check for a specific implementation if requested

if(NOT CP2K_BLAS_VENDOR MATCHES "auto")
  find_package(${CP2K_BLAS_VENDOR} REQUIRED)
  if(TARGET CP2K_${CP2K_BLAS_VENDOR}::blas)
    get_target_property(CP2K_BLAS_INCLUDE_DIRS CP2K_${CP2K_BLAS_VENDOR}::blas
                        INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(CP2K_BLAS_LINK_LIBRARIES CP2K_${CP2K_BLAS_VENDOR}::blas
                        INTERFACE_LINK_LIBRARIES)
    set(CP2K_BLAS_FOUND TRUE)
  endif()
else()
  # search for any blas implementation
  foreach(_libs ${__BLAS_VENDOR_LIST})
    # i exclude the first item of the list
    find_package(${_libs})
    if(TARGET CP2K_${_libs}::blas)
      get_target_property(CP2K_BLAS_INCLUDE_DIRS CP2K_${_libs}::blas
                          INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(CP2K_BLAS_LINK_LIBRARIES CP2K_${_libs}::blas
                          INTERFACE_LINK_LIBRARIES)
      set(CP2K_BLAS_VENDOR "${_libs}")
      set(CP2K_BLAS_FOUND TRUE)
    endif()
  endforeach()
endif()

if(CP2K_BLAS_INCLUDE_DIRS)
  find_package_handle_standard_args(
    Blas REQUIRED_VARS CP2K_BLAS_LINK_LIBRARIES CP2K_BLAS_INCLUDE_DIRS
                       CP2K_BLAS_VENDOR)
else()
  find_package_handle_standard_args(Blas REQUIRED_VARS CP2K_BLAS_LINK_LIBRARIES
                                                       CP2K_BLAS_VENDOR)
endif()

if(NOT TARGET CP2K_BLAS::blas)
  add_library(CP2K_BLAS::blas INTERFACE IMPORTED)
endif()

set_target_properties(CP2K_BLAS::blas PROPERTIES INTERFACE_LINK_LIBRARIES
                                                 "${CP2K_BLAS_LINK_LIBRARIES}")

if(CP2K_BLAS_INCLUDE_DIRS)
  set_target_properties(CP2K_BLAS::blas PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                   "${CP2K_BLAS_INCLUDE_DIRS}")
endif()

mark_as_advanced(CP2K_BLAS_INCLUDE_DIRS)
mark_as_advanced(CP2K_BLAS_LINK_LIBRARIES)
mark_as_advanced(CP2K_BLAS_VENDOR)
mark_as_advanced(CP2K_BLAS_FOUND)
