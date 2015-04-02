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

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <Eigen/Dense>

namespace AMG {
  
  /** \brief A tupel of 3 values, representing a set of coordinates.
   * 
   * The purpose of this type is to allow faster computations by accessing the
   * underlying Eigen::Vector3d type inside functions and algorithms that have
   * sorted out the FrameOfReference issues, i.e. all utilized entities are
   * expressed in the same FrameOfReference.
   */
  using CoordinateTupel = Eigen::Vector3d;
  
  /** \brief A tupel of 4 values, representing a set of quaternions.
   * 
   * Quaternions are normally defined as either <b>q = w + ix + jy + kz</b> or 
   * <b>q = q0 + iq1 + jq2 + kq3</b>. 
   * 
   * Eigen, the underlying matrix library has an own type for quaternions:
   * Eigen::Quaternion::Coefficients. However, those coefficients are
   * internally stored in an unconeventional manner, which could lead to 
   * confusion. Eigen stores the quaternion tupel as [x,y,z,w], which could lead 
   * to an ordering permutation when used in the "flight mechanics way" as
   * 
   * \code
   * double q0=Eigen::Quaternion::Coefficient[3]; // w, the real part
   * double q1=Eigen::Quaternion::Coefficient[0]; // x, the 1st imaginary part
   * double q2=Eigen::Quaternion::Coefficient[1]; // y, the 2nd imaginary part
   * double q3=Eigen::Quaternion::Coefficient[2]; // z, the 3rd imaginary part
   * \endcode
   * 
   * In order to avoid this, an  AMG::QuaternionTupel is defined as the 
   * [q0,q1,q2,q3] tupel of quaternions: 
   * 
   * \code
   * double q0=AMG::QuaternionTupel[0]; // w, the real part
   * double q1=AMG::QuaternionTupel[1]; // x, the 1st imaginary part
   * double q2=AMG::QuaternionTupel[2]; // y, the 2nd imaginary part
   * double q3=AMG::QuaternionTupel[3]; // z, the 3rd imaginary part
   * \endcode
   * 
   * \warning As a result of this 
   *  AMG::QuaternionTupel != Eigen::Quaternion::Coefficients ! (It is a 
   *  permutation of it.)
   */
  using QuaternionTupel = Eigen::Vector4d ;
  
  /** \brief A 3-tupel holding 3-2-1 Euler angles.
   * 
   * This is a 3-tupel holding the roll/phi, pitch/theta, and yaw/psi 
   * Euler-angles.
   * 
   * \note The tupel does not keep track of the used units, i.e. degrees or 
   *  radians! (That's a todo...)
   * 
   * \code
   * AMG::EulerAngleTupel euler;
   * 
   * double phi   = euler[0]; // roll angle
   * double theta = euler[1]; // pitch angle
   * double psi   = euler[2]; // yaw angle
   * \endcode
   * 
   * \todo Make this a template taking Units::degree or Units::radian
   */
  using EulerAngleTupel = Eigen::Vector3d;
  
  
  
}


#endif // TYPEDEFS_H