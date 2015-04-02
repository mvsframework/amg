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

#
# Find the include directory
# =======================
#
find_path(AMG_INCLUDE_DIR amg.hpp PATH_SUFFIXES amg)

#debug output
if( NOT Amg_FIND_QUIETLY )
  if( NOT AMG_INCLUDE_DIR )
    message(SEND_ERROR ${AMG_INCLUDE_DIR})
  else( NOT Amg_FIND_QUIETLY )
    message(STATUS "AMG include directory is " ${AMG_INCLUDE_DIR} )
  endif( NOT AMG_INCLUDE_DIR )
endif( NOT Amg_FIND_QUIETLY )


#
# Find the actual library
# =======================
#

# Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
if( Amg_USE_STATIC_LIBS )
  set( _Amg_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a )
  endif()
endif( Amg_USE_STATIC_LIBS )

# do the finding
find_library(AMG_Amg_LIBRARY amg)

# Restore the original find library ordering
if( Amg_USE_STATIC_LIBS )
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_Amg_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif(Amg_USE_STATIC_LIBS)




# debug output
if( NOT Amg_FIND_QUIETLY )
  if( NOT AMG_Amg_LIBRARY )
    message(SEND_ERROR ${AMG_Amg_LIBRARY})
  else( NOT AMG_Amg_LIBRARY )
    message(STATUS "AMG has been found as " ${AMG_Amg_LIBRARY} )
  endif( NOT AMG_Amg_LIBRARY )
endif( NOT Amg_FIND_QUIETLY )

# find_path(AMG_INCLUDE_DIR xxx.h)
# find_library(AMG_xxx_LIBRARY xxx)
# find_library(AMG_yyy_LIBRARY yyy)


#
# Prepare the standart return
# ===========================
#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMG DEFAULT_MSG
  AMG_INCLUDE_DIR AMG_Amg_LIBRARY)
if(AMG_FOUND)
  set(AMG_INCLUDE_DIRS ${AMG_INCLUDE_DIR})
  set(AMG_LIBRARIES ${AMG_Amg_LIBRARY} )
endif() 


# don't show the internal variables outside the advanced view
mark_as_advanced( 
  AMG_INCLUDE_DIR
  AMG_Amg_LIBRARY
)