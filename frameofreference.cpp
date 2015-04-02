/*
 * Copyright (c) 2011 Claus Christmann.
 * All rights reserved.       
 *   
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  
 * 1) Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *  
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 */


#include "frameofreference.h"
#include "vector.h"
#include "assert.h"
#include <math.h>
#include <iostream>


namespace AMG { 

//
// Some declarations/definitions for statics...
//
  
const FrameOfReference* FrameOfReference::root = nullptr;



//
// ... and the usual stuff
//

FrameOfReference::FrameOfReference()
  : inertial(true)
  , parent(NULL)
  , translationFromParent(Eigen::Vector3d(0.0,0.0,0.0))
  , rotationFromParent(Eigen::Matrix3d::Identity())
{
  assert( NULL==FrameOfReference::root );
    
  if( NULL == FrameOfReference::root )
  { // this is the first time the empty constructor is called
    FrameOfReference::root = this;
  }
  else
  {
    ///\todo Create an error message even in the #define NDEBUG scenario.
  }
}

FrameOfReference::FrameOfReference(const FrameOfReference* other)
  : inertial(other->isInertial())
  , parent(other)
  , translationFromParent(Eigen::Vector3d(0.0,0.0,0.0))
  , rotationFromParent(Eigen::Matrix3d::Identity())
{

}

FrameOfReference::~FrameOfReference()
{

}



bool FrameOfReference::isInertial(void ) const
{
  return inertial;
}

FrameOfReference const * FrameOfReference::inertialRootFrame(void ) 
{
  assert( NULL != root );
  
  return root;
}

const Vector FrameOfReference::positionOfOrigin() const
{
  return Vector(expressInParentFrame(CoordinateTupel::Zero()),parent);
}







const FrameOfReference* FrameOfReference::parentFrame() const
{
  return parent;
}


QuaternionTupel FrameOfReference::attitude() const
{
  //WARNING: Note that Eigen-internally the coefficients are stored in the 
  // following order: [x, y, z, w]  for q=w+xi+yj+zk
  Eigen::Quaterniond::Coefficients qEigen = Eigen::Quaterniond(rotationFromParent).coeffs();
  
  return QuaternionTupel(qEigen[3],qEigen[0],qEigen[1],qEigen[2]);
}

template<>
EulerAngleTupel FrameOfReference::attitude<Units::radian>() const
{
  
  //NOTE: the -1 premultiplication is necessary due to the choice of the 
  //  rotation matrix in Eigen. See the note in 
  //  FrameOfReference::EulerRotationMatrix for more details.
  return -1*rotationFromParent.eulerAngles(0,1,2);
}

template<>
EulerAngleTupel FrameOfReference::attitude<Units::degree>() const
{
 return Units::radian2degree(attitude<Units::radian>()); 
}

/** \todo It should not be necessary to split the handed in coordinate tupel 
 *    when calling the Translation3d constructor
 */
void FrameOfReference::translate(const double& xTranslation, const double& yTranslation, const double& zTranslation)
{
  assert( this != root );

//   std::cout << "Translating a frame" << std::endl;
  
  // convert the translation to parent's coordinates
  CoordinateTupel parentCoords=expressInParentFrame(CoordinateTupel(xTranslation,yTranslation,zTranslation));
  
  // the translations are given individually, use them to create an affine 
  // transformation construct
  Eigen::Translation3d translation( Eigen::Translation3d(parentCoords(0),
                                                    parentCoords(1),
                                                    parentCoords(2)) );
  
  // The converstion of the locally defined translation vector to the parent's
  // coordinates already incorporates the previous offset.
  // Hence the obtained coordinates represent the new total offset between this
  // coordinate system's origin and the parent's origin and can as such be 
  // stored as the total translation from the parent frame.
  translationFromParent = translation;
}

void FrameOfReference::translate(const CoordinateTupel& xyzTupel)
{
  this->translate(xyzTupel(0),xyzTupel(1),xyzTupel(2));
}

void FrameOfReference::translate(const Vector& vector)
{
  this->translate(vector.coordsIn(this));
}

void FrameOfReference::setPositionOfOrigin(const Vector absoluteLocation)
{
  this->translate(absoluteLocation.absoluteCoordsIn(this));
}


// A template specialization for AMG::radian
template<>
void FrameOfReference::rotate<Units::radian>(const double& psi, 
                                             const double& theta, 
                                             const double& phi)
{
  assert( this != root );
  assert( -M_PIl   <= psi   && psi   <  M_PIl   );
  assert( -M_PI_2l <= theta && theta <= M_PI_2l );
  assert( -M_PIl   <= phi   && phi   <  M_PIl   );
  
  
  // "add" the created rotation to the overall affinte transformations to keep
  // track of the overall trasformation
  rotationFromParent = EulerRotationMatrix(psi,theta,phi) * rotationFromParent;
}


// A template specialization for AMG::degree
template<>
void FrameOfReference::rotate<Units::degree>(const double& psi, 
                                             const double& theta, 
                                             const double& phi)
{
  // do the conversion and call the radian function
  rotate<AMG::Units::radian>(Units::degree2radian(psi)   ,
                             Units::degree2radian(theta) ,
                             Units::degree2radian(phi)   );
}

/**
 * \todo Currently only relative attitudes wrt. to the (default) parent
 *  frame are supported.
 */
template<>
void FrameOfReference::setAttitude<Units::radian>(const double psi, const double theta, const double phi, const FrameOfReference* referenceFrame)
{
  if( NULL == referenceFrame )
    referenceFrame = parent;
 
  assert( referenceFrame == parent );
    
  rotationFromParent = EulerRotationMatrix(psi, theta, phi);  
}

/**
 * \todo Currently only relative attitudes wrt. to the (default) parent
 *  frame are supported.
 */
template<>
void FrameOfReference::setAttitude<Units::degree>(const double psi, const double theta, const double phi, const FrameOfReference* referenceFrame)
{
  if( NULL == referenceFrame )
    referenceFrame = parent;
 
  assert( referenceFrame == parent );
    
  rotationFromParent = EulerRotationMatrix(Units::degree2radian(psi)  ,
                                           Units::degree2radian(theta),
                                           Units::degree2radian(phi)  );  
}







void FrameOfReference::setAttitude(const QuaternionTupel relativeAttitude, 
                                   const FrameOfReference* referenceFrame)
{
  double psi_rad, theta_rad, phi_rad;
  quaternions2euler<Units::radian>(relativeAttitude[0],
                                   relativeAttitude[1],
                                   relativeAttitude[2],
                                   relativeAttitude[3],
                                   psi_rad,theta_rad, phi_rad  );
  
  setAttitude<Units::radian>(psi_rad, theta_rad, phi_rad, referenceFrame);
}


CoordinateTupel FrameOfReference::expressInParentFrame(const AMG::CoordinateTupel& tupel) const
{
//   std::cout 
//     << "Original tupel      : "  << tupel.transpose() << "\n"
//     << "Undo the rotations  : " << (rotationFromParent.inverse()*tupel).transpose() << "\n"
//     << "Undo the translation: " << (translationFromParent*rotationFromParent*tupel).transpose() << "\n"
//     << std::endl;
    
   //NOTE:
   // the rotation needs to be inverted, the translation MUST NOT be inverted.
   // The translation essentially captures a vector _from_ the origin of the 
   // parent _to_ the origin of the child. Hence this vector (and not it's 
   // inverse, i.e. -1*vector) simply needs to be added to any vector in the 
   // local frame to get the total vector in the parent frame.
   
//    return translationFromParent*rotationFromParent.inverse()*tupel;
   return transformationToParent() * tupel;
}


// CoordinateTupel FrameOfReference::expressInParentFrame(const double& c1,
//                                                        const double& c2,
//                                                        const double& c3) const
// {
//   return transformationToParent() * CoordinateTupel(c1,c2,c3);
// }


void FrameOfReference::outputRelativePositionToParent(void ) const
{
  std::cout << "All coordinates in the parent frame!\n" 
  << "\tOffset  : " << expressInParentFrame(CoordinateTupel(0.0,0.0,0.0)).transpose() << "\n"
  << "\t(1,0,0) : " << expressInParentFrame(CoordinateTupel(1.0,0.0,0.0)).transpose() << "\n"
  << "\t(0,1,0) : " << expressInParentFrame(CoordinateTupel(0.0,1.0,0.0)).transpose() << "\n"
  << "\t(0,0,1) : " << expressInParentFrame(CoordinateTupel(0.0,0.0,1.0)).transpose() << "\n"
  << std::endl;
}


Eigen::Matrix3d FrameOfReference::EulerRotationMatrix(const double& psi, const double& theta, const double& phi)
{
  using namespace Eigen;
  
  assert( -M_PIl   <= psi   && psi <  M_PIl   );
  assert( -M_PI_2l <= theta && theta <= M_PI_2l );
  assert( -M_PIl   <= phi   && phi <  M_PIl   );
  
  //NOTE: all the angles are multiplied by -1 in order to matcht usage in flight
  // dynamics. The rotation matrix created by eigen represents the matrx for 
  // rotating a vector inside one frame of reference.
  // I.e. yawing a (1,0,0) vector by 90.0deg results a (0,1,0) vector -- which
  // is correct. However, the rotation matrix that we want would have to give us
  // (0,-1,0), i.e. the answer to "how is a north pointing vector, (1,0,0) in
  // in NED, represented in the 90.0deg yawed body frame?"...
  // Right, the answer would be as (0,-1,0)...
 
  return (  Eigen::AngleAxis<double>(-1*phi, Eigen::Vector3d::UnitX())
  * Eigen::AngleAxis<double>(-1*theta, Eigen::Vector3d::UnitY())
  * Eigen::AngleAxis<double>(-1*psi, Eigen::Vector3d::UnitZ()) ).toRotationMatrix();
}


Eigen::Affine3d FrameOfReference::transformationToParent(void ) const
{
  
  //NOTE:
  // the rotation needs to be inverted, the translation MUST NOT be inverted.
  // The translation essentially captures a vector _from_ the origin of the 
  // parent _to_ the origin of the child. Hence this vector (and not it's 
  // inverse, i.e. -1*vector) simply needs to be added to any vector in the 
  // local frame to get the total vector in the parent frame.
  
  return translationFromParent*rotationFromParent.inverse();
}


Eigen::Affine3d FrameOfReference::transformationToRoot(void ) const{
  if( this != root )
  {
    return (parent->transformationToRoot()) * (this->transformationToParent()); 
  }
  else
  {
    return Eigen::Affine3d::Identity();
  }
}


std::ostream& operator<<(std::ostream& stream, const AMG::FrameOfReference& f)
{
  stream
    << "All coordinates in the parent frame!\n" 
    << "\tOffset  : " << f.expressInParentFrame(CoordinateTupel(0.0,0.0,0.0)).transpose() << "\n"
    << "\t(1,0,0) : " << f.expressInParentFrame(CoordinateTupel(1.0,0.0,0.0)).transpose() << "\n"
    << "\t(0,1,0) : " << f.expressInParentFrame(CoordinateTupel(0.0,1.0,0.0)).transpose() << "\n"
    << "\t(0,0,1) : " << f.expressInParentFrame(CoordinateTupel(0.0,0.0,1.0)).transpose() << "\n" ;

    return stream;
}





} // namespace AMG
