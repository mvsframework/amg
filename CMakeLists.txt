# 
# Copyright (c) 2011 Claus Christmann.
# All rights reserved.       
#   
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  
# 1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#  
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# 


project(amg)
cmake_minimum_required(VERSION 2.8)

# load my CMake Stuff
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../CMake/Modules/")
# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/src/3rdParty")
include(cc_general) 

# install into the source dir of this project
set( CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/.. )

# enable C++11
add_definitions(-std=gnu++11)

# show all warnings during the compilation
add_definitions(-Wall)


# Deal with Eigen
# ===============
# Eigen is a matrix library that is required for the current Vroni interface.
# Eigen is a headers-only library, hence no libraries are used. 
#
find_package( Eigen3 REQUIRED QUIET)
if( EIGEN3_FOUND )
#   message(STATUS "Eigen ${EIGEN3_VERSION} has been found at ${EIGEN3_INCLUDE_DIR}")
  include_directories( ${EIGEN3_INCLUDE_DIR} )
  message_variable(EIGEN3_INCLUDE_DIR)
endif( EIGEN3_FOUND )


#
# Deal with Boost
# ===============
# openSuse currently doesn't have a RPM for boost_1_47 hence a small manual 
# "hack". As boost is headers only, no libraries are used.
# NOTE: this is set to the trunk version of a (partially checked out) boost repo
# so that the geometries extensions are available.
#
set( Boost_INCLUDE_DIR "~/3rdParty/boost/1_53_0/" )
include_directories( ${Boost_INCLUDE_DIR} )
message(STATUS "Boost has been found at ${Boost_INCLUDE_DIR}")


set( amg_SRCS  
  frameofreference.cpp
  vector.cpp
  gis.cpp
  rigidbody.cpp
  commoncoordinatesystems.cpp
#  point.cpp
)

set( amg_HDRS
  amg.hpp
  typedefs.hpp
  units.hpp
  frameofreference.h
  vector.h
  gis.h
  commoncoordinatesystems.h
  rigidbody.h
)

add_library_shared_and_static(amg amg_SRCS)

add_executable(test main.cpp )
target_link_libraries(test ${LIB_AMG} )

install(TARGETS ${LIB_AMG_TARGETS}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${archive_DIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${library_DIR}
)

install(FILES ${amg_HDRS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${headers_DIR}/amg 
)

