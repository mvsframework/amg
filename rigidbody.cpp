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


#include "rigidbody.h"

#include "units.hpp"
#include "commoncoordinatesystems.h"
#include "gis.h"
#include <Eigen/Geometry>

namespace AMG {


RigidBody::RigidBody(const AMG::FrameOfReference* parentFrame)
  : FrameOfReference(parentFrame)
  , mass(1.0) // initialize to 1kg
  , inertiaMatrix(Eigen::Matrix3d::Zero()) // initialize to a point mass
  , m_x(0.0)
  , m_y(0.0)
  , m_z(0.0)
  , m_u(0.0)
  , m_v(0.0)
  , m_w(0.0)
  , m_phi(0.0)
  , m_theta(0.0)
  , m_psi(0.0)
  , m_p(0.0)
  , m_q(0.0)
  , m_r(0.0)
  , x(m_x)
  , y(m_y)
  , z(m_z)
  , u(m_u)
  , v(m_v)
  , w(m_w)
  , phi(m_phi)
  , theta(m_theta)
  , psi(m_psi)
  , p(m_p)
  , q(m_q)
  , r(m_r)
{
  inertial=false; // this reference frame is not an inertial one as it could accelerate
}


RigidBody::RigidBody(const double& bodyMass                    ,
                     const Eigen::Matrix3d & bodyInertiaMatrix ,
                     const FrameOfReference* parentFrame       )
  : FrameOfReference(parentFrame)
  , mass(bodyMass)
  , inertiaMatrix(bodyInertiaMatrix)
  , m_x(0.0), x(m_x)
  , m_y(0.0), y(m_y)
  , m_z(0.0), z(m_z)
  , m_u(0.0), u(m_u)
  , m_v(0.0), v(m_v)
  , m_w(0.0), w(m_w)
  , m_phi(0.0), phi(m_phi)
  , m_theta(0.0), theta(m_theta)
  , m_psi(0.0), psi(m_psi)
  , m_p(0.0), p(m_p)
  , m_q(0.0), q(m_q)
  , m_r(0.0), r(m_r)
{
  inertial=false; // this reference frame is not an inertial one as it could accelerate
}


// NOTE: not implemented to dissallow copying.
// RigidBody::RigidBody(const RigidBody& other)
// {
// }


RigidBody::~RigidBody()
{

};




Vector RigidBody::position() const
{
  return FrameOfReference::positionOfOrigin();
}


Vector RigidBody::velocity() const
{
  return Vector(u,v,w,this); //FIXME: well, this is just wrong ;)
}

Vector RigidBody::angularRates() const
{
  return Vector(p,q,r,this); //FIXME: well, this is just wrong ;)
}


void RigidBody::movePositionBy(const Vector& relativeMotion)
{
  translate(relativeMotion.coordsIn(this));
}


void RigidBody::setPositionTo(const Vector& absoluteTargetLocation)
{
  translate(absoluteTargetLocation.absoluteCoordsIn(this));
}


void RigidBody::setVelocityTo(const Vector& velocity)
{
  CoordinateTupel bodySpeeds(velocity.coordsIn(this));
  
  m_u = bodySpeeds[0];
  m_v = bodySpeeds[1];
  m_w = bodySpeeds[2];
}

void RigidBody::setRatesTo(const Vector& rotationalRates)
{
  CoordinateTupel bodyRates(rotationalRates.coordsIn(this));
  
  m_phi = bodyRates[0];
  m_theta = bodyRates[1];
  m_psi = bodyRates[2];
  
}






//NOTE: order matters here! this specialization to <radian> has to appear ABOVE
// the specialization to <degree> as that calls this.
template<>
void RigidBody::setPositionTo<AMG::Units::radian>(const double& longitude, const double& latitude, const double& elevation)
{
  assert( NULL != CoSy::getECEF() );
  
  double x_ecef,y_ecef,z_ecef;
  GIS::geodetic2ecef(longitude,latitude,elevation,x_ecef,y_ecef,z_ecef);
  
  setPositionTo(Vector(x_ecef,y_ecef,z_ecef,CoSy::getECEF()));
}

template<>
void RigidBody::setPositionTo<AMG::Units::degree>(const double& longitude, const double& latitude, const double& elevation)
{
 using Units::degree2radian;
 setPositionTo<Units::radian>(degree2radian(longitude),degree2radian(latitude),elevation);
}




} // end namespace AMG
