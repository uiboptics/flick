#ifndef flick_filter
#define flick_filter

#include "../environment/input_output.hpp"
#include "../numeric/linalg/matrix.hpp"
#include "../numeric/function.hpp"
#include "../numeric/flist.hpp"
#include "../numeric/range.hpp"
#include "../numeric/distribution.hpp"
#include "../environment/exception.hpp"
#include "radiator.hpp"
#include <algorithm>

namespace flick {
  namespace filter {     
    class filter {      
    public:
      virtual ~filter() = default;
      virtual double transmittance(double wavelength) const = 0;
      virtual std::vector<double> extract_range(const std::vector<double>& wl) const {
	return wl;
      }
      template<typename Function>
      Function transmission(const Function& radiation_spectrum) const {
	Function ts = transmittance_spectrum<Function>(radiation_spectrum.x());
	return multiply(ts, radiation_spectrum, ts.x());
      } 
      template<typename Function>
      double weighted_average(const Function& radiation_spectrum) const {
	const Function &f = radiation_spectrum;
	return transmission(f).integral() /
	  transmittance_spectrum<Function>(f.x()).integral();
      }
    private:
      template<typename Function>
      Function transmittance_spectrum(const std::vector<double>& wl) const {
	auto wle = extract_range(wl);
	Function t;
	for (size_t i=0; i < wle.size(); ++i) {
	  t.append({wle[i], transmittance(wle[i])});
	}
	return t;	
      }
    };

    class band_filter : public filter {
    protected:
      double wl0_{500e-9};
      double fwhm_{10e-9};
    public:
      band_filter(double center_wavelength, double fwhm)
	: wl0_{center_wavelength}, fwhm_{fwhm} {
	if (fwhm <= 0)
	  throw std::invalid_argument("band_filter: fwhm must be positive");
      }
    private:
      std::vector<double> extract_range(const std::vector<double>& wl) const override {
	double dwl = 3*fwhm_/2;
	auto first = std::lower_bound(wl.begin(), wl.end(), wl0_-dwl);
	auto last  = std::upper_bound(wl.begin(), wl.end(), wl0_+dwl);
	if (first != wl.begin())
	  --first;
	if (last != wl.end())
	  ++last;
	return std::vector<double>(first, last);
      }
    };
    
    class cut_ends : public filter {
      double l_;
      double u_;
    public:
      cut_ends(double lower_edge, double upper_edge)
	: l_{lower_edge},u_{upper_edge} {}
      double transmittance(double wavelength) const override {
	if (wavelength < l_ or wavelength > u_)
	  return 0;	
	return 1;
      }
    };

    class erythema : public filter
    //  CIE/ISO-standard erythemal action spectrum. Schmalwieser,
    //  A.W., Wallisch, S. and Diffey, B., 2002. A library of action
    //  spectra for erythema and pigmentation. Photochemical &
    //  Photobiological Sciences, 1, pp.251-268.
    {
    public:
      double transmittance(double wavelength) const override {
	if (wavelength < 298e-9)
	  return 1;
	else if (wavelength < 328e-9)
	  return pow(10, 0.094*(298-wavelength*1e9));
	else
	  return pow(10, 0.015*(140-wavelength*1e9));
      }
    };

    class photons : public filter
    // Converts to number of photons per wavelength
    {
    public:
      double transmittance(double wavelength) const override {
	using namespace constants;
	return wavelength / (h * c);
      }
    };
    
    class tabulated : public filter {
      pl_function t_;
      double wavelength_shift_ = 0;
    public:
      tabulated(const pl_function& filter_transmittance)
	: t_{filter_transmittance} {
      }
      double transmittance(double wavelength) const override {
	return t_.value(wavelength-wavelength_shift_);
      }
      void shift(double wavelength) {
	wavelength_shift_ = wavelength;
      }
    };

    template<int Lms_no>
    class cone_lms : public filter
    // Normalized human cone cell spectral sensitivity
    {
      std::string p = path()+"/radiator/filter_data/cie_photometry";
      pl_flist flist = read<pl_flist>(p+"/cie_lms_cf_2deg_1nm.txt");
      pl_function f = flist(Lms_no);
    public:
      cone_lms() {
	f.scale_x(1e-9);
      }
      const std::vector<double>& wavelength_grid() const {
	return f.x();
      }
      double transmittance(double wavelength) const override {
	return f.value(wavelength);
      }
      friend std::ostream& operator<<(std::ostream &os, const cone_lms<Lms_no>& c) {
	os << c.f;
      return os;
      }
    };
    
    template<int Xyz_no>
    class xyz_bar : public filter
    // Color matching functions
    {
      std::string p = path()+"/radiator/filter_data/cie_photometry";
      pl_flist flist = read<pl_flist>(p+"/cie_xyz_1931_2deg.txt");
      pl_function f = flist(Xyz_no);
    public:
      xyz_bar() {
	f.scale_x(1e-9);
      }
      void use_cone_fundamentals() {
	cone_lms<0> L;
	cone_lms<1> M;
	cone_lms<2> S;
	std::vector<std::vector<double>> m =
	  {{1.94735469, -1.41445123, 0.36476327},
	   {0.68990272, 0.34832189, 0},
	   {0, 0, 1.93485343}};
	const std::vector<double>& x = L.wavelength_grid();
	f.clear();
	for (size_t i = 0; i < x.size(); i++) {
	  double y = L.transmittance(x[i]) * m[Xyz_no][0] +
	    M.transmittance(x[i]) * m[Xyz_no][1] +
	    S.transmittance(x[i]) * m[Xyz_no][2];
	  f.append(point(x[i]*1e-9,y));
	}
      }
      double transmittance(double wavelength) const override {
	return f.value(wavelength);
      }
      friend std::ostream& operator<<(std::ostream &os, const xyz_bar<Xyz_no>& c) {
	os << c.f;
	return os;
      }
    };

    class sentinel3 : public filter
    // Use closest normalized sentinel3 OLIC spectral response
    // function around a given user wavelength
    {
      pl_function centers_;
      std::shared_ptr<tabulated> srf_;
      double user_center_wavelength_;
    public:
      sentinel3(double user_center_wavelength)
	: user_center_wavelength_{user_center_wavelength} {
	std::string p = path()+"/radiator/filter_data/sentinel3/srf";
	centers_ = read<pl_function>(p+"/center_wavelength.txt");
	centers_.scale_x(1e-9);
	size_t n = closest_srf();
	pl_function f = read<pl_function>(p+"/band_"+std::to_string(n)+".txt");
	f.scale_x(1e-9);
	f.scale_y(1/f.integral());
	f.add_zero_extrapolation();
	srf_ = std::make_shared<tabulated>(tabulated(f));
	srf_->shift(user_center_wavelength_-centers_.x()[n]);
      }
      double transmittance(double wavelength) const override {
	return srf_->transmittance(wavelength);
      }
      size_t closest_srf() const {
	double n = centers_.value(user_center_wavelength_);
	if (n > centers_.size()-1)
	  return centers_.size()-1;
	if (n < 0)
	  return 0;
	return std::round(n);
      }
    };

      class gaussian : public band_filter {
    public:
      using band_filter::band_filter;
      double transmittance(double wavelength) const override {
	using namespace constants;
	double sigma = fwhm_/(2*sqrt(2*log(2)));
	return 1/(sigma*sqrt(2*pi))*exp(-0.5*pow((wavelength-wl0_)/sigma,2));
      }
    };

    class triangular : public band_filter {
    public:
      using band_filter::band_filter;
      double transmittance(double wavelength) const override {
	return distribution::triangular(wl0_-fwhm_,wl0_+fwhm_,wl0_).pdf(wavelength);
      }
    };

    class square : public band_filter {
    public:
      using band_filter::band_filter;
      double transmittance(double wavelength) const override {
	return cut_ends(wl0_-fwhm_/2, wl0_+fwhm_/2).transmittance(wavelength)/fwhm_;
      }
    };
  }
  
  template<typename Function>
  inline double n_photons(const Function& radiation_spectrum, double wl_low,
			  double wl_high) {
    return filter::photons().transmission(radiation_spectrum).integral(wl_low, wl_high);
  }
  inline double uv_index(const pl_function& radiation_spectrum) {
    return 40*filter::erythema().transmission(radiation_spectrum).integral();
  }
  inline double uva_index(const pl_function& radiation_spectrum) {
    return 40*filter::erythema().transmission(radiation_spectrum).integral(315e-9,400e-9);
  }
  inline double uvb_index(const pl_function& radiation_spectrum) {
    return 40*filter::erythema().transmission(radiation_spectrum).integral(280e-9,315e-9);
  }
  template<typename Function>
  inline std::vector<double> chromaticity(const Function& radiation_spectrum) {
    double wl1 = 380e-9;
    double wl2 = 780e-9;
    const Function& s = radiation_spectrum;
    std::vector<double> xyz(3);
    xyz[0] = filter::xyz_bar<0>().transmission(s).integral(wl1,wl2);
    xyz[1] = filter::xyz_bar<1>().transmission(s).integral(wl1,wl2);
    xyz[2] = filter::xyz_bar<2>().transmission(s).integral(wl1,wl2);
    double sum = 0;
    for (size_t i = 0; i < xyz.size(); i++)
      sum += xyz[i];
    for (size_t i = 0; i < xyz.size(); i++) {
      xyz[i] = xyz[i]/sum;
    }
    return xyz;
  }
  template<typename Function>
  inline std::vector<double> rgb(const Function& radiation_spectrum) {
    using namespace linalg;
    linalg::matrix sRGB_D65 =  // White for CIE D65 spectrum 
      {{3.2404542,-1.5371385,-0.4985314},
       {-0.9692669,1.8760108,0.0415560},
       {0.0556434,-0.2040259,1.0572252}};
    linalg::matrix xyz = {chromaticity(radiation_spectrum)};
    std::vector<double> rgb = linalg::t(sRGB_D65*linalg::t(xyz))[0];
    double max = *std::max_element(rgb.begin(), rgb.end());
    double gamma = 1/2.2;
    for (int i = 0; i < 3; i++) {
      rgb[i] /= max;
      rgb[i] = pow(std::clamp<double>(rgb[i], 0, 1),gamma);
    }
    return rgb;
  }

  template<typename Band_filter, typename Function>
  inline Function smooth(const Function& radiation_spectrum, double fwhm) {
    if (fwhm > 0) {
      auto& f = radiation_spectrum;
      auto& wl = f.x();
      std::vector<double> y(f.size());
      for (size_t i=0; i<wl.size(); ++i) {
	y[i] = Band_filter(wl[i],fwhm).weighted_average(f);
      }
      return Function(wl,y);
    }
    return radiation_spectrum;
  }
}

#endif
