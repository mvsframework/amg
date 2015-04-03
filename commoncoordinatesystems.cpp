/*
 * Copyright (c) 2015 Claus Christmann <hcc |ä| gatech.edu>.  
 *   
 * Licensed under the Apache Lice*nse, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */


#include "commoncoordinatesystems.h"
#include "gis.h"
#include "units.hpp"

#include <assert.h>

namespace {
  AMG::FrameOfReference* ecef_frame = NULL;
  AMG::FrameOfReference* datum_frame = NULL;
}

namespace AMG {
namespace CoSy {

void setECEF(FrameOfReference* frame)
{
  ecef_frame = frame;
}

FrameOfReference* getECEF()
{
  return ecef_frame;
}


void setDatumFrame(FrameOfReference* frame)
{
  datum_frame = frame;
}


FrameOfReference*const getDatumFrame()
{
  assert( datum_frame );
  return datum_frame;
}


template<>
FrameOfReference createNorthEastDown<AMG::Units::radian>(const double& longitude,
                                                         const double& latitude)
{
  // ensure that the ECEF frame has alrady been set.
  assert(ecef_frame);
  assert( -M_PIl   <= longitude && longitude <= M_PIl   );
  assert( -M_PI_2l <= latitude  && latitude <= M_PI_2l );
  
  // the NED frame will be situated on the refernce ellipsoid
  double altitude(0.0); 
  
  // crate a copy of the ECEF which will become the NED frame
  FrameOfReference ned_frame(ecef_frame);
  
  // some temporary coordinates holding the origin of the NED frame.
  double x,y,z;
  
  // compute the desired orign in ecef coordinates
  AMG::GIS::geodetic2ecef(longitude,latitude,altitude,x,y,z);
  
  // translate the frame to the given point
  ned_frame.translate(x,y,z);
  
  if( 0.0L <= latitude )
  { // Northern half, latitude >= 0
    ned_frame.rotate<Units::radian>(longitude+M_PIl,-(M_PI_2l-latitude),M_PIl);
  }
  else
  { // Southern half, latitude is < 0
    ned_frame.rotate<Units::radian>(longitude,-(M_PI_2l+latitude),0.0L);
  }
    
  return ned_frame; // copy by value
}


template<>
FrameOfReference createNorthEastDown<Units::degree>( const double& longitude,
                                                         const double& latitude)
{
 return createNorthEastDown<Units::radian>(Units::degree2radian(longitude),
                                                Units::degree2radian(latitude));
}
  
template<>
FrameOfReference* newNorthEastDown<Units::radian>(const double& longitude,
                                                         const double& latitude)
{
  // ensure that the ECEF frame has alrady been set.
  assert(ecef_frame);
  assert( -M_PIl   <= longitude && longitude <= M_PIl   );
  assert( -M_PI_2l <= latitude  && latitude <= M_PI_2l );
  
  // the NED frame will be situated on the refernce ellipsoid
  double altitude(0.0); 
  
  // crate a copy of the ECEF which will become the NED frame
  FrameOfReference* ned_frame = new FrameOfReference(ecef_frame);
  
  // some temporary coordinates holding the origin of the NED frame.
  double x,y,z;
  
  // compute the desired orign in ecef coordinates
  AMG::GIS::geodetic2ecef(longitude,latitude,altitude,x,y,z);
  
  // translate the frame to the given point
  ned_frame->translate(x,y,z);
  
  // rotate the ned frame
  // 1) pos. latitude is North, i.e. a negative rotation aroung the Y-axis
  ned_frame->rotate<Units::radian>(longitude,-latitude,0.0);
  // 2) Add -90.0deg/-M_PI_2 to turn the "Up-East-North" into an NED frame.
  ned_frame->rotate<Units::radian>(0.0,-M_PI_2l,0.0);
  
  return ned_frame; // copy by value
}

template<>
FrameOfReference* newNorthEastDown<Units::degree>( const double& longitude,
                                                     const double& latitude)
{
  return newNorthEastDown<Units::radian>(Units::degree2radian(longitude),
                                            Units::degree2radian(latitude));
}


  
} // end namespace CoSy  
} // end namespace AMG