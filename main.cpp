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


#include <Eigen/Geometry>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>


#include "amg.hpp"

using namespace Eigen;


int main(int argc, char **argv) {
  
  using namespace AMG;
    
  FrameOfReference* ecef = new FrameOfReference;
  CoSy::setECEF(ecef);
  
  { std::cout << "\nWGS 84:\n" 
              <<   "=======" << std::endl;
    
    using namespace AMG::GIS;
              
    std::cout
    << "WGS84\n"
    << "|- a     : " << WGS84::a << "\n"
    << "|- finv  : " << WGS84::finv << "\n"
    << "|- omega : " << WGS84::omega << "\n"
    << "|- GM    : " << WGS84::GM << "\n"
    << "|\n"
    << "|- f     : " << WGS84::f << "\n"
    << "|- b     : " << WGS84::b << "\n"
    << "|- e2    : " << WGS84::e2 << "\n"
    << "|- ep2   : " << WGS84::ep2 << "\n"
    << "|- g0    : " << WGS84::g0 << "\n"
    << std::endl;
    
  }

  { std::cout << "\nGeodetic to ECEF conversions:\n" 
              <<   "=============================" << std::endl;
              
    auto g2e = [] (double lat_deg, double lon_deg, double alt_m)
    {
      double x_ecef, y_ecef, z_ecef;
      
      GIS::geodetic2ecef(Units::degree2radian(lon_deg),
                        Units::degree2radian( lat_deg),
                        alt_m,x_ecef,y_ecef,z_ecef);
      
      std::cout << "Lat,Lon,Alt :" << lat_deg << ", "
                                   << lon_deg << ", "
                                   << alt_m   << "\n"
                << "X,Y,Z       :" << x_ecef << ", " << y_ecef << ", " << z_ecef 
      << std::endl;
    };
    
    g2e(0,0,0);
    g2e(0,0,1000);
    g2e(45,45,0);
    g2e(90,90,0);
    g2e(33.772022,-84.396188,0);
  }
  
  { std::cout << "\nECEF to Geodetic conversions:\n" 
              <<   "=============================" << std::endl;
              
    auto e2g = [] (double x_ecef, double y_ecef, double z_ecef)
    {
      double lat_rad, lon_rad, alt_m;
      
      GIS::ecef2geodetic( x_ecef,y_ecef,z_ecef,lon_rad,lat_rad,alt_m);
      
      std::cout << "X,Y,Z       :" << x_ecef << ", " << y_ecef << ", " << z_ecef << "\n"
                << "Lat,Lon,Alt :" << Units::radian2degree(lat_rad) << ", "
                                   << Units::radian2degree(lon_rad) << ", "
                                   << alt_m 
      << std::endl;
    };
    


    e2g(518259,-5.28199e+06, 3.52545e+06);


  }
  
  
  FrameOfReference* datum = CoSy::newNorthEastDown<Units::degree>( -84.396188, 33.772022 );
  CoSy::setDatumFrame(datum);
  datum->outputRelativePositionToParent();
  std::cout
    << "Datum attitude (roll, pitch, yaw; in degrees) : "
    << datum->attitude<Units::degree>().transpose() << "\n"
    << std::endl;
  
 
    
    
    
    
    
  double x_ecef,y_ecef,z_ecef;

  GIS::geodetic2ecef(Units::degree2radian(-84.396188),
                      Units::degree2radian( 33.772022),
                      0,x_ecef,y_ecef,z_ecef);
  Vector southCorner(x_ecef,y_ecef,z_ecef,CoSy::getECEF());
  
  double lat,lon,alt;
  GIS::ecef2geodetic(southCorner.coords(0),southCorner.coords(1),southCorner.coords(2),
                     lon,lat,alt);
  
  std::cout
  << "South Corner in Lat, Lon, Alt  : " << Units::radian2degree(lat) << ", "
                                         << Units::radian2degree(lon) << ", "
                                         << alt << "\n"
  << "South Corner in ECEF           : "
  << southCorner.coords().transpose() << "\n"
  << "South Corner in absolute Datum : "
  << southCorner.absoluteCoordsIn(datum).transpose() << "\n"
  << std::endl;
  
  std::cout << AMG::GIS::geodeticPositionString(southCorner) << std::endl;

  GIS::geodetic2ecef(Units::degree2radian(-84.396188),
                      Units::degree2radian( 33.772111),
                      0,x_ecef,y_ecef,z_ecef);
  Vector northCorner(x_ecef,y_ecef,z_ecef,CoSy::getECEF());
  std::cout
  << "North Corner in ECEF           : "
  << northCorner.coords().transpose() << "\n"
  << "North Corner in absolute Datum : "
  << northCorner.absoluteCoordsIn(datum).transpose() << "\n"
  << std::endl;
  
  std::cout << AMG::GIS::geodeticPositionString(northCorner) << std::endl;

  Vector delta = northCorner - southCorner;

  std::cout
  << "Delta Corner in ECEF           : "
  << delta.coords().transpose() << "\n"
  << "Delta Corner in relative Datum : "   
  << delta.coordsIn(datum).transpose() << "\n"
  << std::endl;  
    
    
    
    
    
      
  AMG::RigidBody body;
  body.outputRelativePositionToParent();

  body.setPositionTo<Units::degree>(-84.396188 , 33.772111 , 0.0 );
  body.outputRelativePositionToParent();
  
  body.rotate<Units::degree>(90.0,0,0);
  body.outputRelativePositionToParent();
  
  std::cout
    << "Body Euler attitude (phi, theat, psi):\n"
    << "(degree) : " << body.attitude<Units::degree>().transpose() << "\n"
    << "(radian) : " << body.attitude<Units::radian>().transpose() << "\n"
    << "(quaternions) : " << body.attitude().transpose() << "\n"
    << "yaw (deg) : " << body.yaw<Units::degree>() << "\n"
    << "pitch (deg) : " << body.pitch<Units::degree>() << "\n"
    << "roll (deg) : " << body.roll<Units::degree>() << "\n"
    << std:: endl;

  {  
    using namespace AMG::GIS;
    
    std:: cout << "isInLeftHalfCircle( 0,   1) = " << isInLeftHalfCircle(0,    1) << std::endl;  
    std:: cout << "isInLeftHalfCircle( 0,   0) = " << isInLeftHalfCircle(0,    0) << std::endl;  
    std:: cout << "isInLeftHalfCircle( 0,  -1) = " << isInLeftHalfCircle(0,   -1) << std::endl;  
    std:: cout << "isInLeftHalfCircle(10, 190) = " << isInLeftHalfCircle(10, 190) << std::endl;  
    std:: cout << "isInLeftHalfCircle(10,-190) = " << isInLeftHalfCircle(10,-190) << std::endl;  
    std:: cout << "isInLeftHalfCircle(10,-170) = " << isInLeftHalfCircle(10,-170) << std::endl;  
    std:: cout << "isInLeftHalfCircle(10,-175) = " << isInLeftHalfCircle(10,-175) << std::endl;                
    std:: cout << "isInLeftHalfCircle(10,-165) = " << isInLeftHalfCircle(10,-165) << std::endl;                
    std:: cout << "isInLeftHalfCircle(10, 370) = " << isInLeftHalfCircle(10, 370) << std::endl;                
    std:: cout << "isInLeftHalfCircle(10, 365) = " << isInLeftHalfCircle(10, 365) << std::endl;                
      
  }
  
  { std::cout << "Testing \"headingToCompassPoint()\""<< std::endl;
    for( double hdg = 0; hdg <=360; hdg += 5 )
    { std::cout << hdg << "° => "<< AMG::GIS::headingToCompassPoint(hdg) << std::endl; }
  }
  
  return 0;
}


