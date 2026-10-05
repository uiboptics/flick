#ifndef flick_toa_sun
#define flick_toa_sun

#include <algorithm>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include "../astronomy/sun_position.hpp"
#include "../radiator/toa_solar.hpp"
#include "../radiator/filter.hpp"

namespace flick {
  class toa_sun {
    time_point tp_;
    sun_position sp_;
    double latitude_;
    double longitude_;
    double spectral_width_;
    pp_function spectrum_;
  public:
    toa_sun(const std::string& s) {
      std::istringstream is(s);
      is >> tp_;
      is >> latitude_;
      is >> longitude_;
      is >> spectral_width_;
      if (!is)
	throw std::invalid_argument("toa_sun: expected date, time, latitude, longitude, and spectral width");
      constexpr double to_radians = std::numbers::pi/180;
      sp_ = sun_position(tp_, latitude_*to_radians, longitude_*to_radians);
      spectrum_ = radiator::toa_solar().spectrum();
      spectrum_.scale_y(distance_factor()*angle_factor());
    }
    pp_function spectrum() const {
      return smooth(spectrum_);
    }
    pp_function multiply_with(const pp_function& f) const {
      return smooth(multiply(spectrum_,f,f.x()));
    }
    double zenith_angle() const {
      return sp_.zenith_angle();
    }
  private:
    double distance_factor() const {
      earth_orbit eo(tp_.year(), tp_.day_of_year());
      return pow(constants::au/eo.distance(),2);
    }
    double angle_factor() const {
      return std::max(0.0, cos(sp_.zenith_angle()));
    }
    pp_function smooth(const pp_function& f) const {
      return smooth<filter::gaussian,pp_function>(f,spectral_width_);
    } 
  };
}

#endif
