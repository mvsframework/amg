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


#ifndef POINT_H
#define POINT_H

#include "typedefs.hpp"
#include "vector.h"

namespace AMG {

/** \deprecated Does not provide any benefit over the proper use of AMG::Vector. 
 * 
 * \brief A point in the (3-dimensional) space.
 * The major difference between a \c Point and a \c Vector is the how 
 * coordinates are expressed in different frames:
 * Points incorporate the origin offset of different frames (i.e.
 * FrameOfReference::translate ), whereas 
 * Vectors do not.
 * 
 * \sa Point::coordsIn
 */
class Point
{
public:
  /** \brief Default (empty) constructor */
  Point();
  

  /** \brief Constructing a point at the given coordinates (in the given frame)
   * The given position information is internally stored in a position vector, 
   * i.e. the combination of coordinates and a related frame of reference.
   * 
   * @param[in] positionVector A position vector for this point.
   *  The position vector determines the position of the point by indcating the 
   *  displacement from the the origin of the frame of reference the vector is 
   *  represented in to the actual position of \c this point in space.
   */
  Point(const Vector& positionVector);
  
  
  /** \overload
   * The given parameters are used to generate a position vector.
   * 
   * @param[in] coords The tupel representing the point's coordinates.
   * @param[in] referenceFrame The frame the tupel is in referecen to.
   */
  Point(const CoordinateTupel& coords, const FrameOfReference * referenceFrame);

  
  /** \overload 
   * The given parameters are used to generate a position vector.
   * 
   * @param[in] c1 Value of the 1st coordinate of the point.
   * @param[in] c2 Value of the 2nd coordinate of the point.
   * @param[in] c3 Value of the 3rd coordinate of the point.
   * @param[in] referenceFrame Reference frame the coordinates are relative to.
   */
  Point(const double& c1,
        const double& c2,
        const double& c3,
        const FrameOfReference * referenceFrame);
  
  /** \brief Destructor */
  virtual ~Point();
  
  
  /** \brief Get the position coordinates of this point in a different frame 
   */
  CoordinateTupel absoluteCoordsIn(const FrameOfReference* targetFrame) const;
  
private:
  Vector position;
 
};


} // end namespace AMG

#endif // POINT_H
