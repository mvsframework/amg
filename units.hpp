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

#ifndef UNITS_H
#define UNITS_H

namespace AMG {
  
/** \brief The namespace for all units used in \c AMG 
 * This (headers only) namespace contains empty structs used in templates to 
 * indicate the units of the in-parameters and some related conversion 
 * functions.
 */
namespace Units { 
  
  //
  // angular
  //
  
  /** \brief An empty struct used in templates to indicate the use of degrees*/
  struct degree {};
  
  /** \brief An empty struct used in templates to indicate the use of radians*/
  struct radian {};

  //
  // linear
  //
  /** \brief An empty struct used in templates to indicate the use of meters*/
  struct meter {};
  
  /** \brief An empty struct used in templates to indicate the use of feet*/
  struct foot {};
  
  //
  // Conversions
  //
  /** \brief Convert a degree value to a radian value
   * The internal precision is double.
   */
  template<typename T>
  inline T degree2radian(const T& degreeValue)
  {
    return T(degreeValue*M_PIl/180.0);
  };
  
  /** \brief Convert a radian value to a degree value
   * The internal precision is double.
   */
  template<typename T>
  inline T radian2degree(const T& radianValue)
  {
    return T(radianValue*180.0*M_1_PIl);
  }; 
  
  /** \brief Convert a foot value to a meter value
   * The internal precision is double.
   */
  template<typename T>
  inline T feet2meter(const T& footValue)
  {
    return T(0.3048 * footValue);
  };
  
  /** \brief Convert a meter value to a foot value
   * The internal precision is double.
   */
  template<typename T>
  inline T meter2feet(const T& meterValue)
  {
    return T( meterValue/0.3048);
  };
  
  
} // end namespace Units  
} // end namespace AMG


#endif // UNITS