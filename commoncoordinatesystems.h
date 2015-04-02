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

#ifndef COMMONCOORDINATESYSTEMS_H
#define COMMONCOORDINATESYSTEMS_H

#include "frameofreference.h"
#include "units.hpp"

namespace AMG  {
  
/** \brief The Coordinate Systems related namespace */  
namespace CoSy {
  
  /** \brief Set a pointer to the Earth-centered, Earth-fixed frame of reference
   * Some functions in AMG::Cosy need a referece to the frame that represents
   * the ECEF coordinate system. 
   * (This information is stored in a static variable.)
   *
   * @param[in] frame The frame of reference used as the ECEF frame.
   * 
   * \see ecef_frame
   */
  void setECEF(FrameOfReference * frame);
  
  
  /** \brief Get a pointer to the Earth-centered, Earth-fixed frame of reference
   * 
   * The counterpart to setECEF().
   * 
   * @return frame The frame of reference used as the ECEF frame.
   */
  FrameOfReference* getECEF();
  
  
  /** \brief Set a pointer to the frame used as the Datum
   * (This information is stored in a static variable.)
   * 
   * @param[in] frame The frame of reference used as the Datum frame.
   * 
   * \see datum_frame
   */
  void setDatumFrame(FrameOfReference* frame);
  
  
  /** \brief Get the currently set Datum frame.
   * The counterpart to setDatumFrame().
   *
   * \pre datum_frame != NULL
   * 
   * @return A pointer to the frame of reference used as the Datum frame.
   */
  FrameOfReference* const getDatumFrame();
   
  
  /** \brief Create a NED-frame at the given location.
   * 
   * Create a new North-East-Down (NED) frame at the given location (at an 
   * implied altitude of \c 0m, i.e. on the surfarce of the WGS84 referece
   * ellipsoid.
   * 
   * \pre The static \c ecef_frame is not \c NULL
   * \pre -M_PIl   <= longitude <= M_PIl
   * \pre -M_PI_2l <= latitude  <= M_PI_2l
   * 
   * @param[in] longitude The longitude of the datum of the NED frame in <units>.
   * @param[in] latitude The latitude of the datum of the NED frame in <units>.
   * @return A copy by value of the newly created NED frame.
   */
  template<typename units=AMG::Units::degree>
  FrameOfReference createNorthEastDown(const double& longitude  ,
                                         const double& latitude ); 
  
  /** \brief Create a NED-frame at the given location.
   * 
   * Create a new North-East-Down (NED) frame at the given location (at an 
   * implied altitude of \c 0m, i.e. on the surfarce of the WGS84 referece
   * ellipsoid.
   * 
   * \pre The static \c ecef_frame is not \c NULL
   * \pre -M_PIl   <= longitude <= M_PIl
   * \pre -M_PI_2l <= latitude  <= M_PI_2l
   * 
   * @param[in] longitude The longitude of the datum of the NED frame in <units>.
   * @param[in] latitude The latitude of the datum of the NED frame in <units>.
   * @return A pointer to the the newly created NED frame.
   * 
   * \note This function calls \b new FrameOfReference, so watch out for memory
   *  leaks. It's your resonsibiliyt to clean up after yourself!
   */
  template<typename units=AMG::Units::degree>
  FrameOfReference* newNorthEastDown(const double& longitude  ,
                                       const double& latitude );
  
  
} // end namespace CoSy
} // end namespace AMG




#endif // COMMONCOORDINATESYSTEMS_H
