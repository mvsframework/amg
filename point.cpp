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

#include "point.h"

namespace AMG {
  
Point::Point()
{

}


Point::Point(const Vector& positionVector)
  : position(positionVector)
{
  
}


Point::Point(const CoordinateTupel& coords, const FrameOfReference* referenceFrame)
  : position(coords,referenceFrame)
{

}


Point::Point(const double& c1, const double& c2, const double& c3, const FrameOfReference* referenceFrame)
  : position(c1,c2,c3,referenceFrame)
{

}


CoordinateTupel Point::absoluteCoordsIn(const FrameOfReference* targetFrame) const
{
  return position.absoluteCoordsIn(targetFrame);
}

  
} // end namespace AMG
