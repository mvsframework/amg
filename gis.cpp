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


#include "gis.h"

#include <assert.h>
#include <cmath>
#include "units.hpp"
#include "commoncoordinatesystems.h"
#include <boost/concept_check.hpp>

// Anoymous namespace to hide the chosen implementation for ecef2geodetic()
namespace {

using namespace AMG::GIS;  
  

void ecef2geodetic_iterativ(const double& x, const double& y, const double& z, double& lambda, double& phi, double& h)
{
  // NOTE: this function is a replication of ecef2geodetic.m from Matlab
  //
  //   function [phi, lambda, h] = ecef2geodetic(x, y, z, ellipsoid)
  //   %ECEF2GEODETIC Convert geocentric (ECEF) to geodetic coordinates
  //   %
  //   %   [PHI, LAMBDA, H] = ECEF2GEODETIC(X, Y, Z, ELLIPSOID) converts point
  //   %   locations in geocentric Cartesian coordinates, stored in the
  //   %   coordinate arrays X, Y, Z, to geodetic coordinates PHI (geodetic
  //   %   latitude in radians), LAMBDA (longitude in radians), and H (height
  //   %   above the ellipsoid). The geodetic coordinates refer to the
  //   %   reference ellipsoid specified by ELLIPSOID (a row vector with the
  //   %   form [semimajor axis, eccentricity]). X, Y, and Z must use the same
  //   %   units as the semimajor axis;  H will also be expressed in these
  //   %   units.  X, Y, and Z must have the same shape; PHI, LAMBDA, and H
  //   %   will have this shape also.
  //   %
  //   %   For a definition of the geocentric system, also known as
  //   %   Earth-Centered, Earth-Fixed (ECEF), see the help for GEODETIC2ECEF.
  //   %
  //   %   See also ECEF2LV, GEODETIC2ECEF, GEOCENTRIC2GEODETICLAT, LV2ECEF.
  // 
  //   % Copyright 2005-2009 The MathWorks, Inc.
  //   % $Revision: 1.1.6.4 $  $Date: 2009/04/15 23:34:43 $
  // 
  //   % Reference
  //   % ---------
  //   % Paul R. Wolf and Bon A. Dewitt, "Elements of Photogrammetry with
  //   % Applications in GIS," 3rd Ed., McGraw-Hill, 2000 (Appendix F-3).
  // 
  //   % Implementation Notes from Rob Comer
  //   % -----------------------------------
  //   % The implementation below follows Wolf and DeWitt quite literally,
  //   % with a few important exceptions required to ensure good numerical
  //   % behavior:
  //   %
  //   % 1) I used ATAN2 rather than ATAN in the formulas for beta and phi.  This
  //   %    avoids division by zero (or a very small number) for points on (or
  //   %    near) the Z-axis.
  //   %
  //   % 2) Likewise, I used ATAN2 instead of ATAN when computing beta from phi
  //   %    (conversion from geodetic to parametric latitude), ensuring
  //   %    stability even for points at very high latitudes.
  //   %
  //   % 3) Finally, I avoided dividing by cos(phi) -- also problematic at high
  //   %    latitudes -- in the calculation of h, the height above the ellipsoid.
  //   %    Wold and Dewitt give
  //   %
  //   %                   h = sqrt(X^2 + Y^2)/cos(phi) - N.
  //   %
  //   %    The trick is to notice an alternative formula that involves division
  //   %    by sin(phi) instead of cos(phi), then take a linear combination of the
  //   %    two formulas weighted by cos(phi)^2 and sin(phi)^2, respectively. This
  //   %    eliminates all divisions and, because of the identity cos(phi)^2 +
  //   %    sin(phi)^2 = 1 and the fact that both formulas give the same h, the
  //   %    linear combination is also equal to h.
  //   %
  //   %    To obtain the alternative formula, we simply rearrange
  //   %
  //   %                   Z = [N(1 - e^2) + h]sin(phi)
  //   %    into
  //   %                   h = Z/sin(phi) - N(1 - e^2).
  //   %
  //   %    The linear combination is thus
  //   %
  //   %        h = (sqrt(X^2 + Y^2)/cos(phi) - N) cos^2(phi)
  //   %            + (Z/sin(phi) - N(1 - e^2))sin^2(phi)
  //   %
  //   %    which simplifies to
  //   %
  //   %      h = sqrt(X^2 + Y^2)cos(phi) + Zsin(phi) - N(1 - e^2sin^2(phi)).
  //   %
  //   %    From here it's not hard to verify that along the Z-axis we have
  //   %    h = Z - b and in the equatorial plane we have h = sqrt(X^2 + Y^2) - a.
  
  // % Ellipsoid constants
  // WGS84::a   % Semimajor axis
  // WGS84::e2  % Square of first eccentricity
  // WGS84::ep2 % Square of second eccentricity
  // WGS84::f   % Flattening
  // WGS84::b   % Semiminor axis



  // Longitude
  lambda = std::atan2(y,x);

  // Distance from the z-axis
  double rho = std::hypot(x,y);
    
  // Bowring's formula for initial parametric (beta) and geodetic (phi) latitudes
  double beta = std::atan2(z, (1 - WGS84::f) * rho);
  phi = std::atan2(z   + WGS84::b*WGS84::ep2 *std::pow(std::sin(beta),3) ,
                   rho - WGS84::a*WGS84::e2  *std::pow(std::cos(beta),3) );

  // Fixed-point iteration with Bowring's formula
  // (typically converges within two or three iterations)
  double betaNew = std::atan2((1 - WGS84::f)*std::sin(phi), std::cos(phi));
  int count = 0;
  while( beta != betaNew && count < 10 )
  {
    beta = betaNew;
    phi = std::atan2(z   + WGS84::b * WGS84::ep2 * std::pow(std::sin(beta),3) ,
                     rho - WGS84::a * WGS84::e2  * std::pow(std::cos(beta),3) );
    betaNew = std::atan2((1 - WGS84::f)*std::sin(phi), std::cos(phi));
    ++count;
  }

  // Calculate ellipsoidal height from the final value for latitude
  double sinphi = std::sin(phi);
  double N = WGS84::a / std::sqrt(1.0 - WGS84::e2 * sinphi*sinphi);
  h = rho * std::cos(phi) + (z + WGS84::e2 * N * sinphi) * sinphi - N;

}


void ecef2geodetic_closedForm(const double& x, const double& y, const double& z, double& lambda, double& phi, double& h)
{
  //NOTE: This method is taken from
  //  
  // Markku Heikkinen, "Geschlossene Formeln zur Berechnung räumlicher 
  // geodätischer Koordinaten aud rechtwinkeligen Koordinaten", Zeitschrift für 
  // Vermessungswesen, Vol. 107, 5/1982, p. 207--2011, in German
  
  //NOTE: This method is only valid if the radius r3 = sqrt(x*x+y*y+z*z) is 
  // larger than ~43 km. (See the above mentioned paper for the related details.)
  assert( std::sqrt(x*x+y*y+z*z) > 45e3 );
  
  //
  // Ellipsoidal constants
  //
  long double a = WGS84::a;
  long double b = WGS84::b;
  long double e2 = WGS84::e2; // pos. definite
  long double ep2 = WGS84::ep2;// pos. definite
  long double E2 = e2*a*a; // alternative E2 = a*a-b*b; // pos. definite

  //
  // intermediaries
  //
  
  long double r  = std::sqrt( (x*x)+(y*y) );

  if( 0.0L == r ) // on N-S axis.
  { 
    lambda = 0.0;
   
    phi =  static_cast<double>(M_PI_2l);
    if( z < 0 )
    { phi *= -1.0; }
   
    h = std::abs(z)-b;
    return;
  }
  
  long double F  = 54.0L*(b*b)*(z*z); // pos. definite
  long double G  = (r*r)+(1.0L-e2)*(z*z)-e2*E2; // pos. definite
  long double c  = (e2*e2)*F*(r*r)/(G*G*G)  ;
  long double s  = std::pow( 1.0L+c+std::sqrt((c*c)+2.0L*c) ,(1.0L/3.0L)); // 3rd root 
  long double P  = F/( 3.0L*(1.0L+s+1.0L/s)*(1.0L+s+1.0L/s)*(G*G) );
  long double Q  = std::sqrt( 1.0L+2.0L*(e2*e2)*P );
  long double r0 = ((P*e2*r)/(1.0L+Q))
    +std::sqrt( 0.5L*(a*a)*(1.0L+1.0L/Q)-(P*(1.0L-e2)*(z*z))/(Q*(1.0L+Q))-0.5L*P*(r*r) );
  long double U  = std::sqrt( (z*z)+(r-e2*r0)*(r-e2*r0) );
  long double V  = std::sqrt( (r-e2*r0)*(r-e2*r0)+(1.0L-e2)*(z*z) );
  long double z0 = ((b*b)*z)/(a*V);

  //NOTE: I introduced these 3 variables in order to save some typing
  long double W  = z+ep2*z0; 
  long double sinphi = W/std::sqrt( (r*r)+(W*W) );
  long double cosphi = r/std::sqrt( (r*r)+(W*W) );
  
  //
  // final computations
  //
  h = double( U*(1.0L-(b*b)/(a*V)) );
  
  //NOTE: as both, sinphi and cosphi are given in the algorithm, I decided to 
  // invert the lower value, hoping that the higher slope in that region gives
  // a better numerical result. Obviously, this is just a gutt feel, no real
  // scienece involved... 
  if( sinphi<=cosphi )
  {
    phi = double(std::asin(sinphi));
  }
  else
  {
    phi = double(std::acos(cosphi));
  }

  //NOTE: same argument as above... Claus
  if( x<=y )
  {
    lambda = double(std::acos(x/r));
  }
  else
  {
    lambda = double(std::asin(y/r));
  }
  
}

}

namespace AMG { namespace GIS {
  
std::string convertToString ( const CompassPoint& cp )
{
  using CP=AMG::GIS::CompassPoint;
  switch( cp )
  {
    case CP::N    : return "N"   ; //"North";
    case CP::NbE  : return "NbE" ; //"North by East";
    case CP::NNE  : return "NNE" ; //"North-northeast";
    case CP::NEbN : return "NEbN"; //"Northeast by North";
    case CP::NE   : return "NE"  ; //"Northeast";
    case CP::NEbE : return "NEbE"; //"Northeast by East";
    case CP::ENE  : return "ENE" ; //"East-northeast";
    case CP::EbN  : return "EbN" ; //"East by North";
    case CP::E    : return "E"   ; //"East";
    case CP::EbS  : return "EbS" ; //"South";
    case CP::ESE  : return "ESE" ; //"East-southeast";
    case CP::SEbE : return "SEbE"; //"Southeast by East";
    case CP::SE   : return "SE"  ; //"Southeast";
    case CP::SEbS : return "SEbS"; //"Southeast by South";
    case CP::SSE  : return "SSE" ; //"South-southeast";
    case CP::SbE  : return "SbE" ; //"South by East";
    case CP::S    : return "S"   ; //"South";
    case CP::SbW  : return "SbW" ; //"South by West";
    case CP::SSW  : return "SSW" ; //"South-southwest";
    case CP::SWbS : return "SWbS"; //"Southwest by South";
    case CP::SW   : return "SW"  ; //"Southwest";
    case CP::SWbW : return "SWbW"; //"Southwest by West";
    case CP::WSW  : return "WSW" ; //"West-southwest";
    case CP::WbS  : return "WbS" ; //"West by South";
    case CP::W    : return "W"   ; //"West";
    case CP::WbN  : return "WbN" ; //"West by North";
    case CP::WNW  : return "WNW" ; //"West-northwest";
    case CP::NWbW : return "NWbW"; //"Northwest by West";
    case CP::NW   : return "NW"  ; //"Northwest";
    case CP::NWbN : return "NWbN"; //"Northwest by North";
    case CP::NNW  : return "NNW" ; //"North-norhtwest";
    case CP::NbW  : return "NbW" ; //"North by West";
    
    default: return std::string();
  }
  return std::string();
}

double ensureHeading ( const double& heading_deg )
{
 // check for circular multiplicity, i.e. reduce heading to within one circle
 double hdg_deg = std::fmod(heading_deg,360.0);
 
 if( hdg_deg > 180.0 )
 {hdg_deg -=360.0;}
 else if( hdg_deg <= -180.0 )
 {hdg_deg += 360; }
 
 return hdg_deg;
}

double opposingHeading ( const double& heading_deg )
{
  return ensureHeading(heading_deg+180.0);
}
  
double compassPointToHeading ( CompassPoint const & cp )
{
  double center_deg = static_cast<double>(cp)/100;
  return ensureHeading(center_deg);
}

bool isInLeftHalfCircle ( const double& reference_deg, const double& heading_deg )
{
  double ref_deg  = ensureHeading(reference_deg );
  double test_deg = ensureHeading(heading_deg);
  
  //NOTE: remember: left  is [ref-180,ref), i.e. the opposite heading is in 
  // both, the left and the right half circle!
  
  if( ref_deg == test_deg )
  { return false; }
  
  if( ref_deg == 0.0 )
  { return test_deg < 0  or test_deg == 180.0; } 
  
  if( ref_deg > 0.0 )
  { return  (ref_deg-180.0)<=test_deg && test_deg < ref_deg; } 
  
  if( ref_deg < 0.0 )
  { return !( ref_deg<test_deg && test_deg <(ref_deg+180.0) ) ; }
  //NOTE: ---------------------------------^
  // This is not a "<=" as we defined left as [ref-180,ref) !
  
  if( test_deg == opposingHeading (ref_deg) )
  { return true; }
  
    
  // catch all that shouldn't be reached.
  return false;
}

bool isInRightHalfCircle ( const double& reference_deg, const double& heading_deg )
{
  double ref_deg  = ensureHeading(reference_deg );
  double test_deg = ensureHeading(heading_deg);
  
  //NOTE: remember: right  is (ref,ref+180], i.e. the opposite heading is in 
  // both, the left and the right half circle!
  
  if( ref_deg == test_deg )
  { return false; }
  
  if( ref_deg == 0.0 )
  { return test_deg > 0  or test_deg == 180.0; } 
  
  if( ref_deg > 0.0 )
  { return  !( (ref_deg-180.0)<test_deg && test_deg < ref_deg ); } 
  //NOTE: --------------------^
  // This is not a "<=" as we defined right as (ref,ref+180] !
  
  if( ref_deg < 0.0 )
  { return ref_deg<test_deg && test_deg <=(ref_deg+180.0) ; }
  
  if( test_deg == opposingHeading (ref_deg) )
  { return true; }
  
  // catch all that shouldn't be reached.
  return false;
}
  
CompassPoint headingToCompassPoint ( const double& heading_deg, int resolution )
{
  if( resolution != 4 or resolution != 8 /*or resolution != 16*/ )
  { resolution = 8; }
  
  
  double hdg_deg = ensureHeading(heading_deg);
  
  using CP=AMG::GIS::CompassPoint;  
  CP cp = CP::N;
  
  auto inbetween = [hdg_deg](CP left,CP right)
  { return compassPointToHeading(left)<hdg_deg && hdg_deg <= compassPointToHeading(right); };

  if( inbetween(CP::SW, CP ::NW ) )
  { cp = CP::W; 
     
    if( resolution == 4 ){ /*break;*/ }
    else if( inbetween(CP::SSW, CP::WSW) )
    { cp = CP::SW; }
    else if( inbetween(CP::WNW,CP::NNW) )
    { cp = CP::NW; }
      
  }
  else if( inbetween(CP::NW, CP::NE) )
  { cp = CP::N;
    
    if( resolution == 4 ){ /*break;*/ }
    else if( inbetween(CP::WNW,CP::NNW) )
    { cp = CP::NW; }
    else if( inbetween(CP::NNE,CP::ENE) )
    { cp = CP::NE; }
  
  }
  else if( inbetween(CP::NE,CP::SE) )
  { cp = CP::E;
    
    if( resolution == 4 ){ /*break;*/ }
    else if( inbetween(CP::NNE,CP::ENE) )
    { cp = CP::NE; }
    else if( inbetween(CP::ESE,CP::SSE) )
    { cp = CP::SE; }
  
  }
  else
  { cp = CP::S;
    
    if( resolution == 4 ){ /*break;*/ }
    else if( inbetween(CP::ESE,CP::SSE) )
    { cp = CP::SE; }
    else if( inbetween(CP::SSW,CP::WSW) )
    { cp = CP::SW; }
  
  }  
  
    
  return cp; 
}

std::string geodeticPositionString ( const AMG::Vector& position )
{
  CoordinateTupel ecefCoords = position.absoluteCoordsIn( AMG::CoSy::getECEF() );
                  
  double lat_rad,lon_rad,alt_m;
  GIS::ecef2geodetic(ecefCoords[0],ecefCoords[1],ecefCoords[2],
                          lon_rad,lat_rad,alt_m );

  CompassPoint northSouth = CompassPoint::N;
  if ( lat_rad < 0)
  { // Southern hemisphere
    northSouth = CompassPoint::S;
          lat_rad *= -1.0;
  }

  CompassPoint eastWest = CompassPoint::E;
  if ( lon_rad < 0)
  { // Western hemisphere
    eastWest = CompassPoint::W;
          lon_rad *= -1.0;
  }
   
  std::ostringstream stream;
  stream << AMG::Units::radian2degree( lat_rad ) << "°" << convertToString(northSouth) << ", "
         << AMG::Units::radian2degree( lon_rad ) << "°" << convertToString(eastWest)   << ", "
         << alt_m << " m WGS84";
  return stream.str();
}



void ecef2geodetic(const double& x, const double& y, const double& z, double& lambda, double& phi, double& h)
{
  ecef2geodetic_iterativ(x,y,z,lambda,phi,h);
//   ecef2geodetic_closedForm(x,y,z,lambda,phi,h); // forward call to closed form algorithm.
}
  


void geodetic2ecef(const double& lambda, const double& phi, const double& h, double& x, double& y, double& z)
{
  // NOTE: this function is a replication of ecef2geodetic.m from Matlab
  //
  //   function [x, y, z] = geodetic2ecef(phi, lambda, h, ellipsoid)
  //   %GEODETIC2ECEF Convert geodetic to geocentric (ECEF) coordinates
  //   %
  //   %   [X, Y, Z] = GEODETIC2ECEF(PHI, LAMBDA, H, ELLIPSOID) converts geodetic
  //   %   point locations specified by the coordinate arrays PHI (geodetic
  //   %   latitude in radians), LAMBDA (longitude in radians), and H (ellipsoidal
  //   %   height) to geocentric Cartesian coordinates X, Y, and Z.  The geodetic
  //   %   coordinates refer to the reference ellipsoid specified by ELLIPSOID (a
  //   %   row vector with the form [semimajor axis, eccentricity]).  H must use
  //   %   the same units as the semimajor axis;  X, Y, and Z will be expressed in
  //   %   these units also.
  //   %
  //   %   The geocentric Cartesian coordinate system is fixed with respect to the
  //   %   Earth, with its origin at the center of the ellipsoid and its X-, Y-,
  //   %   and Z-axes intersecting the surface at the following points:
  //   %
  //   %                PHI  LAMBDA
  //   %      X-axis:    0     0      (Equator at the Prime Meridian)
  //   %      Y-axis:    0   pi/2     (Equator at 90-degrees East)
  //   %      Z-axis:  pi/2    0      (North Pole)
  //   %
  //   %   A common synonym is Earth-Centered, Earth-Fixed coordinates, or ECEF.
  //   %
  //   %   See also ECEF2GEODETIC, ECEF2LV, GEODETIC2GEOCENTRICLAT, LV2ECEF.
  //   
  //   % Copyright 2004-2009 The MathWorks, Inc.
  //   % $Revision: 1.1.6.4 $  $Date: 2009/04/15 23:34:46 $
  //   
  //   % Reference
  //   % ---------
  //   % Paul R. Wolf and Bon A. Dewitt, "Elements of Photogrammetry with
  //   % Applications in GIS," 3rd Ed., McGraw-Hill, 2000 (Appendix F-3).
  
  // % Ellipsoid constants
  // WGS84::a   % Semimajor axis
  // WGS84::e2  % Square of first eccentricity
  // WGS84::ep2 % Square of second eccentricity
  // WGS84::f   % Flattening
  // WGS84::b   % Semiminor axis
  
  double sinphi = std::sin(phi);
  double cosphi = std::cos(phi);
  
  double N  = WGS84::a / std::sqrt(1.0 - WGS84::e2 * sinphi*sinphi);
//   double N  = WGS84::e2 ;
//   N *= sinphi;
//   N *= sinphi;
//   N = 1.0 - N;
//   N = std::sqrt(N);
//   N = WGS84::a / N;
  
  x = (N + h) * cosphi * std::cos(lambda);
  y = (N + h) * cosphi * std::sin(lambda);
  z = (N*(1 - WGS84::e2) + h) * sinphi;
}


void GudermannFunction(const double& yGudermann, double& phiGeocentric)
{
  phiGeocentric = std::asin( std::tanh( yGudermann ) );
}



void InverseGudermannFunction(const double& phiGeocentric, double& yGudermann)
{
  assert( (-M_PI_2l < phiGeocentric) && (phiGeocentric < M_PI_2l) );
  yGudermann = std::asinh( std::tan( phiGeocentric ) );
}



void geocentric2mercator(const double& lambdaGeocentric,
                                    const double& phiGeocentric,
                                    const double& hGeocentric, 
                                          double& xGudermann, 
                                          double& yGudermann, 
                                          double& zGudermann)
{
  // longitude
  xGudermann = lambdaGeocentric;
  // latitude
  InverseGudermannFunction(phiGeocentric,yGudermann);
  // elevation
  zGudermann = hGeocentric;
}

void mercator2geocentric(const double& xGudermann,
                                    const double& yGudermann,
                                    const double& zGudermann, 
                                          double& lambdaGeocentric,
                                          double& phiGeocentric,
                                          double& hGeocentric)
{
  // longitude
  lambdaGeocentric = xGudermann;
  // latitude
  GudermannFunction(yGudermann,phiGeocentric);
  // elevation
  hGeocentric=zGudermann;  
}


void cartesian2spherical(const double& x,
                                   const double& y,
                                   const double& z,
                                         double& azimuth,
                                         double& elevation,
                                         double& distance)
{
  double xy = std::hypot(x,y); // length of the projection in the xy-plane
  distance  = std::hypot(z,xy);
  azimuth   = std::atan2(y,x);
  elevation = std::atan2(z,xy);
}


void spherical2cartesian(const double& azimuth,
                                   const double& elevation,
                                   const double& distance,
                                         double& x,
                                         double& y,
                                         double& z)
{
  assert(-M_PIl <= azimuth && azimuth <= M_PIl);
  assert(-M_PI_2l <= elevation && elevation <= M_PI_2l);
  assert( 0.0 <= distance );
  
  
  double xy = std::cos(elevation) * distance; // length of the projection in the xy-plane
  z = std::sin(elevation) * distance;
  y = std::sin(azimuth) * xy;
  x = std::cos(azimuth) * xy;
}

void ecef2geocentric(const double& x, const double& y, const double& z, double& lambdaGeocentric, double& phiGeocentric, double& hGeocentric)
{
  // treat it as a regular cartesian2spherical conversion
  double radialDistance(0.0);
  cartesian2spherical(x,y,z,lambdaGeocentric,phiGeocentric,radialDistance);
  
  // subtract the raddius of the referecen sphere from the distance
  hGeocentric = radialDistance-WGS84::a;
}


void geocentric2ecef(const double& lambdaGeocentric, const double& phiGeocentric, const double& hGeocentric, double& x, double& y, double& z)
{
  double radialDistance = hGeocentric+WGS84::a;
  spherical2cartesian(lambdaGeocentric,phiGeocentric,radialDistance,x,y,z);
}


void geocentric2geodetic(const double& lambdaGeocentric, const double& phiGeocentric, const double& hGeocentric, double& lambda, double& phi, double& h)
{
 double xTemp,yTemp,zTemp;
 
 geocentric2ecef(lambdaGeocentric,phiGeocentric,hGeocentric,xTemp,yTemp,zTemp);
 ecef2geodetic(xTemp,yTemp,zTemp,lambda,phi,h);
}


void geodetic2geocentric(const double& lambda, const double& phi, const double& h, double& lambdaGeocentric, double& phiGeocentric, double& hGeocentric)
{
  double xTemp,yTemp,zTemp;
 
  geodetic2ecef(lambda, phi, h, xTemp, yTemp, zTemp);
  ecef2geocentric(xTemp,yTemp,zTemp,lambdaGeocentric,phiGeocentric,hGeocentric);
}


void mercator2ecef(const double& xGudermann, const double& yGudermann, const double& zGudermann, double& x, double& y, double& z)
{
  double lambdaGeocentric, phiGeocentric, hGeocentric;
  
  mercator2geocentric(xGudermann,yGudermann,zGudermann,lambdaGeocentric,phiGeocentric,hGeocentric);
  geocentric2ecef(lambdaGeocentric,phiGeocentric,hGeocentric,x,y,z);
  
}

void ecef2mercator(const double& x, const double& y, const double& z, double& xGudermann, double& yGudermann, double& zGudermann)
{
  double lambdaGeocentric, phiGeocentric, hGeocentric;
  
  ecef2geocentric(x,y,z,lambdaGeocentric,phiGeocentric,hGeocentric);
  geocentric2mercator(lambdaGeocentric,phiGeocentric,hGeocentric,xGudermann,yGudermann,zGudermann);
}


void mercator2geodetic(const double& xGudermann, const double& yGudermann, const double& zGudermann, double& lambda, double& phi, double& h)
{
  double lambdaGeocentric, phiGeocentric, hGeocentric;
  
  mercator2geocentric(xGudermann,yGudermann,zGudermann,lambdaGeocentric,phiGeocentric,hGeocentric);
  geocentric2geodetic(lambdaGeocentric,phiGeocentric,hGeocentric,lambda,phi,h);

}

void geodetic2mercator(const double& lambda, const double& phi, const double& h, double& xGudermann, double& yGudermann, double& zGudermann)
{
  double lambdaGeocentric, phiGeocentric, hGeocentric;
  
  geodetic2geocentric(lambda,phi,h,lambdaGeocentric,phiGeocentric,hGeocentric);
  geocentric2mercator(lambdaGeocentric,phiGeocentric,hGeocentric,xGudermann,yGudermann,zGudermann);
}




} // end namespace GIS
} // end namespace AMG





std::ostream& operator<< ( std::ostream& out, const CompassPoint& cp )
{
  return out << AMG::GIS::convertToString(cp);
}

