# 
# Copyright (c) 2015 Claus Christmann <hcc |ä| gatech.edu>.  
#   
# Licensed under the Apache Lice*nse, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
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

