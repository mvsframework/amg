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


#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>

#include <Eigen/Dense> // this is most likely total overkill for an inclusion...

#include "typedefs.hpp"
#include "frameofreference.h"



/** \brief The Analytical Mechanics/Geometry (AMG) namespace.
 */
namespace AMG {

  
/** \brief A 3D vector, i.e. a coordinate tuple in a cartesian coordinate system.
 * 
 * \todo Implement a way to permanently change the frame of reference used 
 *  internally.
 */
class Vector 
{

public:
  //NOTE: to see what this is and why it's necessary, read
  // http://eigen.tuxfamily.org/dox-devel/TopicStructHavingEigenMembers.html
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  //
  // Constructors
  //
  /** \brief Default (empty) constructor.
   * 
   * Initializes the coordiante tuple to (0,0,0). The related referece frame is
   * invalid (NULL).
   */
  Vector();
  
  /** \brief Convenience constructor taking individual arguments.
   * 
   * \todo Look into the const-correctness: Shouldn't the paramter 
   *  \c referenceFrame be a pointer to a const FrameOfReference?
   *
   * 
   * @param[in] c1 Value of the 1st coordinate of the vector.
   * @param[in] c2 Value of the 2nd coordinate of the vector.
   * @param[in] c3 Value of the 3rd coordinate of the vector.
   * @param[in] referenceFrame Reference frame the coordinates are relative to.
   */
  Vector(const double& c1,
         const double& c2,
         const double& c3,
         const FrameOfReference * referenceFrame);
  
  /** \overload
   * 
   * @param[in] coords The coordinate values for this vector
   * @param[in] referenceFrame The frame of reference \c coords are expressed in
   */
  Vector(const CoordinateTupel& coords, const FrameOfReference * referenceFrame);
  
  
  /** \brief Copy constructor. */
  Vector(const Vector& other);

  
  //
  // Operators
  //
  /** \brief Assingment operator.
   * 
   * Usage example: 
   * \code
   *  Vector a,b;
   *  a = b; // identical to a.operator=(b); 
   * \endcode
   * @param[in] other The Vector beeing assinged to \c *this.
   * @return \c *this, by definition
   */
  virtual Vector& operator=(const Vector& other);
  
  /** \brief Additive compound assingment operator.
   * 
   * \note This operator ignores non-cogruent origins, i.e. it internally uses
   *  Vector::corrdsIn() and not Vector::absoluteCoordsIn()!
   * 
   * Usage example:
   * \code
   * Vector a,b;
   * a += b; // identical to a = a+b; or a.operator+=(b);
   * \endcode
   */
  virtual Vector& operator+=(const Vector& rhs);
  
  /** \brief Subtractive compound assingment operator.
   * 
   *  \note This operator ignores non-cogruent origins, i.e. it internally uses
   *  Vector::corrdsIn() and not Vector::absoluteCoordsIn()!
   *  
   * Usage example:
   * \code
   * Vector a,b;
   * a -= b; // identical to a = a-b; or a.operator-=(b);
   * \endcode
   */
  virtual Vector& operator-=(const Vector& rhs);
  
  
  /** \brief Multiplicative compound assingment operator.
   * 
   * \note This operator ignores non-cogruent origins, i.e. it internally uses
   *  Vector::corrdsIn() and not Vector::absoluteCoordsIn()!
   *  
   * Usage example:
   * \code
   * Vector a;
   * double b;
   * a *= b; // identical to a = a*b; or a.operator*=(b);
   * \endcode
   */
  virtual Vector& operator*=(const double& rhs);
  
  
  /** \brief Dividing compound assingment operator.
   * 
   * \note This operator ignores non-cogruent origins, i.e. it internally uses
   *  Vector::corrdsIn() and not Vector::absoluteCoordsIn()!
   * 
   * Usage example:
   * \code
   * Vector a;
   * double b;
   * a /= b; // identical to a = a/b; or a.operator/=(b);
   * \endcode
   * 
   * \pre 0.0 != rhs
   */
  virtual Vector& operator/=(const double& rhs);
  
  /** \brief Vector addition. 
   * 
   * \note This operator ignores non-cogruent origins, i.e. it internally uses
   *  Vector::corrdsIn() and not Vector::absoluteCoordsIn()!
   * 
   */
  virtual Vector  operator+(const Vector& rhs) const;
  
  /** \brief Vector subtraction.
   * 
   * \note This operator ignores non-cogruent origins, i.e. it internally uses
   *  Vector::corrdsIn() and not Vector::absoluteCoordsIn()!
   * 
   */
  virtual Vector  operator-(const Vector& rhs) const;
  
  /** \brief Scalar vector multiplicaion. */
  virtual Vector operator*(const double& rhs) const;
  
  /** \brief Scalar vector division. */
  virtual Vector operator/(const double& rhs) const;
    
  /** \brief Comparison operator for equality. */
  virtual bool    operator==(const Vector& other) const;
  
  /** \brief Comparison operator for inequality. */
  virtual bool    operator!=(const Vector& other) const;
  
  /** \brief Outputs the vector to the given stream. */
  friend std::ostream& operator<< ( std::ostream&  stream, const Vector& v );   
  
  //
  // Accessors
  //
  /** \brief Get a reference to the coordinate tupel. */
  CoordinateTupel& coords(void);
  CoordinateTupel const& coords(void) const;
  
  /** \brief Get one of the three coordinate values 
   * \pre 0<= dimension && dimension <=3
   *
   * \param[in] dimension 0-indexed indicator for the requested dimension.
   * \return The corresponding coordinate value.
   */
  const double& coords(const int& dimension) const;
  
  
  
  /** \brief Get a coordinate tupel expressing this vector in a different frame.
   * 
   * This function returns a coordinate tupel (i.e. \b not a reference) 
   * expressing \c this Vector in a different frame of reference.
   * 
   * \note The function \b does \b not take translations between the two frames 
   *  into account, only rotations! (Use ::absoluteCoordsIn() to do that.)
   * 
   * Note also that the returned CoordinateTupel does not carry any reference to
   * the targetFrame. It's your job to keep track of which FrameOfReference a 
   * CoordinateTupel belongs to.
   * 
   * @param[in] targetFrame The frame this vector's coordinates should be expressed in
   * @return The corresponding set of coordinates.
   */
  CoordinateTupel coordsIn(const FrameOfReference* targetFrame) const;
  
  /** \brief Get the coordinates of the point this vector points to wrt. a different frame
   * 
   * Unlike coordsIn() this function takes translations between the 
   * FrameOfReference this Vector is defined in and the targetFrame into 
   * account.
   * 
   * Note that the returned CoordinateTupel does not carry any reference to the 
   * targetFrame. It's your job to keep track of which FrameOfReference a 
   * CoordinateTupel belongs to.
   * 
   * \todo Rename to "tipCoordsIn()" or "headCoordsIn()" ???
   * 
   * @param[in] targetFrame The frame the point this vector points to should be expressed in.
   * @return The corresponding set of coordinates.
   */
  CoordinateTupel absoluteCoordsIn(const FrameOfReference* targetFrame) const;
  
  /** \brief Get a pointer to the frame of reference this vector is defined in.
   */
  const FrameOfReference* frame(void) const;
  
  /** \brief Set the frame of reference this vector is defined in.
   * 
   * Changes the frame of refernce this vector's coordinates are valid in.
   * 
   * \note This function does \b not perform any alterations in the coordinate tupel
   * of this vector if transfromCoords==false, (the default!) 
   * If transfromCoords is set to TRUE, the transformation to frame is done 
   * using absoluteCoordsIn().
   *  
   * \param[in] frame The new frame of reference this vector is defined in.
   * \param[in] transfromCoords (optional) FALSE if the coordinates should not 
   *   be changed, TRUE otherwise.
   */
  void setFrameOfReference(const FrameOfReference* frame, bool transfromCoords = false);
  
  
  
  
  /** \brief Return the (Euclidian) length of the vector.
   * \return The Euclidian norm of this vector, ||Vector||
   */
  double norm() const;
  
   
private:
  /** \brief The actual container holding the coordinates. */
  CoordinateTupel coordinateTupel; 
  /** \brief A pointer to the frame of reference this vector is defined in. */
  const FrameOfReference * frameOfReference; 
  
};



} // end namespace AMG

#endif // VECTOR_H
