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
    filter::sentinel3 sentinel3_srf_;
    bool use_sentinel3_srf_ = false;
  public:
    toa_sun() = default;
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
    void use_sentinel3_srf(bool use) {
      use_sentinel3_srf_ = use;
    }
    pp_function spectrum() const {
      return band_average(spectrum_);
    }
    pp_function multiply_with(const pp_function& f) const {
      const auto& wl = spectrum_.x();
      auto first = std::lower_bound(wl.begin(), wl.end(), f.x().front());
      auto last  = std::upper_bound(wl.begin(), wl.end(), f.x().back());
      return band_average(multiply(spectrum_,f,std::vector(first,last)));
    }
    double zenith_angle() const {
      return sp_.zenith_angle();
    }
    double azimuth_angle() const {
      return sp_.azimuth_angle();
    }
  private:
    double distance_factor() const {
      earth_orbit eo(tp_.year(), tp_.day_of_year());
      return pow(constants::au/eo.distance(),2);
    }
    double angle_factor() const {
      return std::max(0.0, cos(sp_.zenith_angle()));
    }
    pp_function band_average(const pp_function& f) const {
      if (use_sentinel3_srf_) {
	const auto& g = moving_average(f,sentinel3_srf_);
	auto& cwl = sentinel3_srf_.center_wavelengths();
	double to_mW_per_nm = 1e-6;
	return pp_function(cwl, g.value(cwl)).scale_y(to_mW_per_nm);
      }
      return moving_average(f,filter::gaussian(spectral_width_));
    } 
  };
}

#endif
