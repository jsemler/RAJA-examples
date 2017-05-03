###############################################################################
#
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-689114
#
# All rights reserved.
#
# This file is part of RAJA.
#
# For additional details, please also read raja/README-license.txt.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

add_definitions(-DUSE_MPI=0 -DUSE_OMP=1)

# Build options
option(LULESH_ENABLE_MPI "Build Lulesh with MPI" Off)

if (LULESH_ENABLE_MPI)
  set(MPI_CXX_SKIP_FLAGS
    -DMPI_NO_CPPBIND -DOMPI_SKIP_MPICXX -D_MPICC_H -DMPICH_SKIP_MPICXX -DMPIBULL_SKIP_MPICXX
    CACHE STRING "Default skip-flags used (for SGI, OpenMPI, HP-MPI, MPICH, BULL-MPI)")
  find_package(MPI REQUIRED)
  add_definitions(-DUSE_MPI=1 -DUSE_OMP=0)
endif()

if (LULESH_ENABLE_MPI)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_INCLUDE_PATH}
    ${MPI_CXXFLAGS})
else()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

if (RAJA_ENABLE_CUDA)
  add_definitions(-DUSE_CASE=9)
endif ()

set (lulesh_sources
  lulesh.cc
  lulesh-comm.cc
  lulesh-init.cc
  lulesh-util.cc
  lulesh-viz.cc)

raja_add_executable(
  NAME lulesh2.0_RAJA-ISet.exe
  SOURCES ${lulesh_sources}
  DEPENDS_ON ${RT_LIBRARIES})

if (LULESH_ENABLE_MPI)
  target_link_libraries(lulesh2.0_RAJA-ISet.exe ${MPI_LIBRARIES})
endif()
