#include "sun_position.hpp"
#include "../environment/unit_test.hpp"

namespace flick {
  begin_test_case(sun_position_test_A) {
    double to_radians = std::numbers::pi/180;
    double to_degrees = 180/std::numbers::pi;

    time_point tp_utc = {2024,1,1,22,0,0};
    check_close(tp_utc.day_of_year(),22./24,1e-7_pct);
    check_close(tp_utc.hour_of_day(),22,1e-7_pct);
    
    earth_orbit eo = {2024,tp_utc.day_of_year()};
    check_close(eo.distance(),1.47102016e11,0.01_pct);
    check_close(eo.declination()*to_degrees,-23.03,0.1_pct);
 
    // Bergen
    double latitude = 60.3925*to_radians;
    double longitude = 5.323333*to_radians; 
    sun_position sp(tp_utc,latitude ,longitude);
    check_close(sp.zenith_angle()*to_degrees,(90+48.45),1_pct);
    check_close(sp.azimuth_angle()*to_degrees,322.5-180,1_pct);

    // Negative azimuth before noon
    sun_position sp2(time_point{2026,6,1,10,0,0},latitude ,longitude);
    double aa = sp2.azimuth_angle()*to_degrees;
    check(aa < 0 && aa > -90);

    // Positive azimuth after noon
    sun_position sp3(time_point{2026,6,1,13,0,0},latitude ,longitude);
    aa = sp3.azimuth_angle()*to_degrees;
    check(aa > 0 && aa < 90);
  } end_test_case()

  begin_test_case(sun_position_test_B) {
    double to_radians = std::numbers::pi/180;
    double to_degrees = 180/std::numbers::pi;
    double earth_tilt = 23.4;
 
    // North pole at northern summer solstice
    time_point tp_utc = {2023,6,21,3,57,0};
    double latitude = 90*to_radians;
    double longitude = 0*to_radians; 
    sun_position sp(tp_utc,latitude ,longitude);
    check_close(sp.zenith_angle()*to_degrees,90-earth_tilt,0.1_pct);

    // South pole at northern winter solstice
    tp_utc = {2023,12,22,3,27,0};
    latitude = -90*to_radians;
    longitude = 0*to_radians; 
    sp = sun_position{tp_utc,latitude ,longitude};
    check_close(sp.zenith_angle()*to_degrees,90-earth_tilt,0.1_pct);
  } end_test_case()

  begin_test_case(sun_position_test_C) {
    // Compare with satellite prodcut values
    double to_radians = std::numbers::pi/180;
    double to_degrees = 180/std::numbers::pi;
    time_point tp_utc = {2022,5,20,10,52,0};
 
    // Hardangerfjord
    double latitude = 59.87*to_radians;
    double longitude = 5.68*to_radians;

    // Sentinel-3 OLCI geometry
    double SZA = 40.55;  // zero for overhead sun
    double SAA = 164.6;  // zero for sun in north, positive clockwise
    double OZA = 14.2;   // zero when looking vertically down
    double OAA = -67.85; // zero when looking south, positive clockwise
      
    // Convert to Flick geometry
    double SZA_F = SZA;     // zero for overhead sun
    double SAA_F = SAA-180; // zero for sun in south, positive clockwise
    double OZA_F = 180-OZA; // zero when looking vertically up
    double OAA_F = OAA-SAA_F; // zero when looking towards the sun, positive clockwise   
      
    sun_position sp(tp_utc,latitude ,longitude);
    check_close(sp.zenith_angle()*to_degrees, SZA_F, 1_pct);
    check_close(sp.azimuth_angle()*to_degrees, SAA_F, 1.5_pct);
    check(OZA_F > 90 && OZA_F <= 180);
    check(OAA_F < 0 && OAA_F > -90);
   
  } end_test_case()
}
