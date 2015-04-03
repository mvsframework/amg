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


#include "vector.h"

namespace AMG {


Vector::Vector()
  : coordinateTupel(0.0,0.0,0.0)
  , frameOfReference(NULL)
{

}


Vector::Vector(const double& c1,
               const double& c2,
               const double& c3,
               const FrameOfReference * referenceFrame)
 : coordinateTupel(c1,c2,c3)
 , frameOfReference(referenceFrame)
{

}


Vector::Vector(const CoordinateTupel& coords, const FrameOfReference* referenceFrame)
  : coordinateTupel(coords)
  , frameOfReference(referenceFrame)
{

}


Vector::Vector(const Vector& other)
  : coordinateTupel(other.coordinateTupel)
  , frameOfReference(other.frameOfReference)
{

}



Vector& Vector::operator=(const Vector& other)
{
  if( this != &other )
  {
    this->frameOfReference = other.frameOfReference;
    this->coordinateTupel = other.coordinateTupel;
  }
  return *this;
}

Vector& Vector::operator+=(const Vector& rhs)
{
  coordinateTupel += rhs.coordsIn(frameOfReference);
  return *this;
}

Vector& Vector::operator-=(const Vector& rhs)
{
  coordinateTupel -= rhs.coordsIn(frameOfReference);
  return *this;
}

// scalar multiplication
Vector& Vector::operator*=(const double& rhs)
{
  coordinateTupel *= rhs;
  return *this;
}

Vector& Vector::operator/=(const double& rhs)
{
  assert( 0.0 != rhs );
  return *this *= 1.0/rhs;
}


Vector Vector::operator+(const AMG::Vector& rhs) const
{
  Vector result = *this;  // create the resulting Vector as a copy of *this
  result += rhs;          // use the compound operator to do the work
  return result;          // return by value 
}

Vector Vector::operator-(const AMG::Vector& rhs) const
{
  Vector result = *this;  // create the resulting Vector
  result -= rhs;          // use the compound operator to do the work
  return result;          // return by value 
}

Vector Vector::operator*(const double& rhs) const
{
  Vector result = *this;
  return result *= rhs;
}

Vector Vector::operator/(const double& rhs) const
{
  Vector result = *this;
  return result /= rhs;
}



/** \todo Deal with numerical precision issues when comparing transformed 
 *    coordinates.
 */
bool Vector::operator==(const Vector& other) const
{
  if( this == &other ) return true;
  if( this->coordinateTupel == other.coordsIn(this->frame()) ) return true; //TODO: numerical epsilon
  
  return false;
}

bool Vector::operator!=(const Vector& other) const
{
  return !(*this == other);
}

std::ostream& operator<<(std::ostream& stream, const AMG::Vector& vector)
{
  stream << vector.coordinateTupel << "\n";
  return stream;
} 


//
// Accessors
//

CoordinateTupel& Vector::coords(void )
{
  return coordinateTupel;
}

CoordinateTupel const& Vector::coords(void ) const
{
  return coordinateTupel;
}



const double& Vector::coords(const int& dimension) const
{
  assert(0<= dimension && dimension <=3);
  return coordinateTupel[dimension];
}



CoordinateTupel Vector::coordsIn(const FrameOfReference* targetFrame) const
{
  if( targetFrame == frameOfReference)
  {
    return coordinateTupel;
  }
  else
  {

    Eigen::Affine3d transformationToRoot =  
      frameOfReference->transformationToRoot();
      
    Eigen::Affine3d transformationFromRootToTargetFrame =  
      (targetFrame->transformationToRoot()).inverse();
    
    //NOTE: the ().linear() is necessary to take translations out of the picture      
    return (transformationFromRootToTargetFrame*transformationToRoot).linear()*coordinateTupel;
  }
}

CoordinateTupel Vector::absoluteCoordsIn(const FrameOfReference* targetFrame) const
{
  if( targetFrame == frameOfReference)
  {
    return coordinateTupel;
  }
  else
  {
    
    Eigen::Affine3d transformationToRoot =  
    frameOfReference->transformationToRoot();
    
    Eigen::Affine3d transformationFromRootToTargetFrame =  
    (targetFrame->transformationToRoot()).inverse();
    
    return transformationFromRootToTargetFrame*transformationToRoot*coordinateTupel;
  }
}


const FrameOfReference* Vector::frame(void ) const
{
  return frameOfReference;
}

void Vector::setFrameOfReference(const FrameOfReference* frame, bool transformCoords)
{
  if( transformCoords )
    coordinateTupel = absoluteCoordsIn(frame);
  frameOfReference = frame;
}


double Vector::norm() const
{
  return coordinateTupel.norm();
}




} // namespace AMG
