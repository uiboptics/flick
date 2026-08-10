#ifndef flick_material_ocean
#define flick_material_ocean

#include "water/pure_water.hpp"
#include "water/cdom.hpp"
#include "water/phytoplankton.hpp"
#include "water/nap.hpp"
#include "marine_particles/marine_particles.hpp"
#include "marine_cdom/marine_cdom.hpp"
#include "mixture.hpp"
#include "spheres.hpp"
#include "../environment/configuration.hpp"
#include <algorithm>
#include <cctype>

namespace flick {
namespace material {
  struct ocean : public mixture<pl_function> {
    struct configuration : public mixture::configuration {
      configuration() {
	add<double>("bottom_depth", 200, R"(
Total depth of the water column [m].
)");
	
	add<double>("concentration_relative_depths", {0,1}, R"(
Space-separated list of depth fractions (df), ranging from the surface
(df = 0) to the bottom (df = 1), defining depths of the material
scaling factor profile. See the `concentration_scaling_factors`
variable.
)");
	
	add<double>("concentration_scaling_factors", {1,1}, R"(
Space-separated list of scaling factors that scale the concentration
of all ocean materials except pure water, sea ice with brines and
bubbles, and those listed by the concentration_exception_names
variable.
)");      
	
	add<double>("concentration_exception_scaling_factors", {1,1}, R"(
Space-separated list of scaling factors that scale the concentration
of all ocean materials listed by the `concentration_exception_materials`
variable.
)");

	add<std::string>("concentration_exception_names", " ", R"(
Space-separated list of material names that is scaled by the factors
listed by the concentration_exception_scaling_factors. Valid material
names are: cdom, phytoplankton, nap, water_bubbles, and any names
listed by the variables `mp_names` or `mcdom_names`.
)");
	
	add<double>("cdom_440", 0.0, R"(
CDOM absorption coefficient at 440 nm [1/m].
)");

	add<double>("cdom_slope", 0.017, R"(
Slope of the CDOM absorption spectrum 1/nm. Note that
this unit is an exception to the SI mks unit convention.
)");
		    	
	add<double>("chl_concentration", 0, R"(
Chlorophyll concentration in the water column [kg/m^3]. A
concentration of e.g., 10.0 mg/m^3 may be written as 10.0e-6 kg/m^3
for clarity.
)");

	add<double>("nap_concentration", 0, R"(
Dry mass concentration of nonalgal particles in the water column
[kg/m^3]. A concentration of e.g., 10.0 g/m^3 may be written as
10.0e-3 kg/m^3 for clarity.
)");
	
	add<double>("bubble_volume_fraction", 0, R"(
Bubble volume fraction in water [unitless].
)");
	
	add<std::string>("bubble_calcualator", {"full_mie"}, R"(
Select `full_mie` or `parameterized_mie` for calculation of water
bubble IOPs with a full Mie code or a fast parameterized Mie code
optimized for large spheres.
)");

	add<double>("bubble_radius", 1e-7, R"(
Median radius of the log-normal size distribution [m]. This median
radius is exp(mu), where mu is the mean of the natural logarithm of the radius
distribution.
)");
	
	add<double>("bubble_sigma", 0.1, R"(
Size-distribution sigma. Sigma is the standard deviation of the
natural logarithm of the radius distribution.
)");

	add<double>("water_temperature", 290, R"(
Temperature in the water column [K].
)");

	add<double>("water_salinity", 30, R"(
Salinity of the water column [PSU].
)");

	add<std::string>("mp_names", "SD16_VF17", R"(
Space-separated list of names of measured marine particles with
inherent optical properties tabulated in separate ASCII files stored in
the Flick directory material/marine_particles/iop_tables.
)");
	
	add<double>("mp_concentrations", 0, R"(
Space-separated list of dry mass concentrations [kg/m^3]
of measured marine particles with inherent optical properties
tabulated in separate ASCII files stored in the Flick directory
material/marine_particles/iop_tables. One concentration value must be
provided for each material specified in `mp_names`.

Concentrations may be written in scientific notation for clarity, e.g.,
10.0 g/m^3 as 10.0e-3 kg/m^3.
)");
	
	add<double>("mp_scattering_scaling_factors", 1, R"(
Space-separated list of scaling factors [unitless] for
manual scaling of the scattering coefficient of marine particles. One
scaling factor must be provided for each marine particle specified in
`mp_names`.
)");
	
	add<double>("mp_bleaching_factors", 0, R"(
Space-separated list of factors [unitless] representing
the degree of particle bleaching. One factor must be provided for each
marine particle specified in mp_names.

A value of '0' corresponds to full absorption, while a value of '1'
corresponds to absorption after bleaching. Values greater than '1'
reduce the absorption beyond the bleached state.
)");
	
	add<std::string>("mcdom_names", "ECOSENS_HF22_D1", R"(
Space-separated list of names of measured marine CDOM with absorption
coefficients tabulated in separate ASCII files stored in the Flick
directory material/marine_cdom/iop_tables.
)");
	
	add<double>("mcdom_scaling_factors", 0, R"(
Space-separated list of scaling factor for measured marine CDOM
absorption coefficients listed in separated ASCII files in the Flick
directory material/marine_cdom/iop_table, one concentration value for
each CDOM spectra given in mcdom_names.
)");
	
	add<int>("ice_depths", 0, R"(
Number of depths, defined by the concentration_relative_depths
variable, that include sea ice.
)");
	
	add<double>("ice_bubble_fraction", {0.01,0.01}, R"(
Space-separated list of the volume fraction of bubble inclusions in
sea ice at each depth defined by the concentration_relative_depths
variable.
)");
	
	add<double>("ice_bubble_radius", 200e-6, R"(
Radius of sea ice bubble inclusions [m].
)");
	
	add<double>("ice_brine_fraction", {0.02,0.02}, R"(
Space-separated list of the volume fraction of saline brine pocket
inclusions in the sea ice at each of the depths defined by the
concentration_relative_depths variable.
)");
	
	add<double>("ice_brine_radius", 500e-6, R"(
Radius of sea ice brine pocket inclusions [m].
)");
      }
    };
  private:
    basic_configuration c_;
  public:
    ocean(const basic_configuration& c=ocean::configuration())
      : mixture(angle_range(c.get<size_t>("n_angles")),height_grid(c)) {
      c_ = c;    
      auto_update_iops(false);
      add_sea_ice();
      add_pure_water();
      add_cdom();
      add_phytoplankton();
      add_nap();
      add_water_bubbles();
      add_marine_particles();
      add_marine_cdom(); 
      auto_update_iops(true);
    }
    static stdvector height_grid(const basic_configuration& c) {
      double epsilon = 1e-6;      
      double depth = c.get<double>("bottom_depth");
      stdvector z = absolute_depth(c);
      stdvector h = {-depth};
      for (size_t i = 0; i < z.size(); ++i) {
	if (z[i] < -2*epsilon and z[i] > -depth+epsilon)
	  h.push_back(z[i]);
      }
      h.push_back(-epsilon);
      return h;
    }
  private:
    void add_sea_ice() {
      size_t n_ice_depths = c_.get<int>("ice_depths");
      stdvector bubble_fraction = c_.get_vector<double>("ice_bubble_fraction");
      stdvector brine_fraction = c_.get_vector<double>("ice_brine_fraction");
      if (n_ice_depths > 0) {
	size_t n_total = c_.get_vector<double>("concentration_relative_depths").size();
	ensure(n_ice_depths <= n_total,"ice_depths");
	ensure(bubble_fraction.size()==n_ice_depths,"bubble_fraction size");
	ensure(brine_fraction.size()==n_ice_depths,"brine_fraction size");

	auto m1 = std::make_shared<pure_ice>();
	add_profile(m1, stdvector(n_ice_depths,1.0), "pure ice"); 

	double width = 0;
	double r_bu = c_.get<double>("ice_bubble_radius");
	using bubbles = bubbles_in_ice<parameterized_monodispersed_mie>;
	auto m2 = std::make_shared<bubbles>(1,log(r_bu),width);
	add_profile(m2, bubble_fraction, "ice bubbles");
	
	double r_br = c_.get<double>("ice_brine_radius");
	using brines = brines_in_ice<parameterized_monodispersed_mie>;
	double salinity = 100;
	auto m3 = std::make_shared<brines>(1,log(r_br),width,salinity);
	add_profile(m3, brine_fraction, "ice brines");
      }
    }
    void add_pure_water() {
      double S = c_.get<double>("water_salinity");
      double T = c_.get<double>("water_temperature");
      size_t n_ice_depths = c_.get<int>("ice_depths");
      size_t n_total = c_.get_vector<double>("concentration_relative_depths").size();
      stdvector vol_frac(n_ice_depths, 0.0);
      vol_frac.resize(n_total, 1.0);
      auto m = std::make_shared<pure_water>(S,T);
      add_profile(m,vol_frac,"pure water");
    }
    void add_cdom() {
      double a440 = c_.get<double>("cdom_440");
      if (a440 > 0) {
	auto m = std::make_shared<cdom>(a440, c_.get<double>("cdom_slope"));
	add_concentration_profile(m, "cdom");
      }
    }
    void add_profile(const std::shared_ptr<base>& m, stdvector factor,
		     const std::string& name) {
      stdvector z = absolute_depth(c_);
      factor.resize(z.size(), 0.0);
      std::reverse(factor.begin(),factor.end());
      add_material(make_scaled_z_profile<pl_function>(m,z,factor),name);
    }
    void add_concentration_profile(const std::shared_ptr<base>& m,
				   const std::string& name) {
      auto f = c_.get_vector<double>("concentration_scaling_factors");
      auto f_e = c_.get_vector<double>("concentration_exception_scaling_factors");
      auto n_e = c_.get_vector<std::string>("concentration_exception_names");
      if (std::ranges::find(n_e, name) != n_e.end())
	add_profile(m,f_e,name);
      else
	add_profile(m,f,name);
    }  
    static stdvector absolute_depth(const basic_configuration& c) {
      double bottom_depth = c.get<double>("bottom_depth");
      stdvector relative_depth = c.get_vector<double>("concentration_relative_depths");
      std::reverse(relative_depth.begin(),relative_depth.end());
      return (-1)*relative_depth*bottom_depth;
    }
    void add_phytoplankton() {
      double chl = c_.get<double>("chl_concentration");
      if (chl > 0) {
	add_concentration_profile(std::make_shared<phytoplankton>(chl),"phytoplankton");
      }
    }
    void add_nap() {
      double con = c_.get<double>("nap_concentration");
      if (con > 0) {
	add_concentration_profile(std::make_shared<nap>(con),"nap");
      }
    }
    void add_water_bubbles() {
      double volume_fraction = c_.get<double>("bubble_volume_fraction");
      if (volume_fraction > 0) {
	std::string calculator = c_.get<std::string>("bubble_calculator");
	double mu = log(c_.get<double>("mie_bubble_radius"));
	double sigma = c_.get<double>("mie_bubble_sigma");
	double S = c_.get<double>("water_salinity");
	double T = c_.get<double>("water_temperature");
	using full = bubbles_in_water<monodispersed_mie>;
	using param = bubbles_in_water<parameterized_monodispersed_mie>;
	if (calculator == "full_mie")	  
	  add_concentration_profile(std::make_shared<full>(volume_fraction,mu,sigma,S,T),"water_bubbles");
	else if (calculator == "parameterized_mie")
	  add_concentration_profile(std::make_shared<param>(volume_fraction,mu,sigma,S,T),"water_bubbles");
	else
	  ensure(false,"bubble_calculator");
      }
    }
    void add_marine_particles() {
      std::vector<std::string> names = c_.get_vector<std::string>("mp_names");
      std::vector<double> concentrations = c_.get_vector<double>("mp_concentrations");
      std::vector<double> scattering_scaling_factors = c_.get_vector<double>("mp_scattering_scaling_factors");
      std::vector<double> bleaching_factors = c_.get_vector<double>("mp_bleaching_factors");
      for (size_t i = 0; i < names.size(); i++) {
	if (concentrations.at(i) > 0) {
	  auto m = std::make_shared<marine_particles>(names.at(i),
						      at_or_last(concentrations,i),
						      at_or_last(scattering_scaling_factors,i),
						      at_or_last(bleaching_factors,i));
	  auto name = "marine_particles_"+names[i];
	  add_concentration_profile(m,name);
	}
      }
    }
    void add_marine_cdom() {
      std::vector<std::string> names = c_.get_vector<std::string>("mcdom_names");
      stdvector scaling_factors = c_.get_vector<double>("mcdom_scaling_factors");
      ensure(names.size()==scaling_factors.size(), "marine cdom");
      for (size_t i = 0; i<names.size(); i++) {
	if (scaling_factors.at(i) > 0) {
	  auto m = std::make_shared<marine_cdom>(names[i], at_or_last(scaling_factors,i));
	  auto name = "marine_cdom_"+names[i];
	  add_concentration_profile(m,name);
	}
      }
    }
  private:
    template<class T>
    T at_or_last(const std::vector<T>& v, size_t i) {
      if (i > v.size()-1)
	return v.back();
      return v[i];
    }
    void ensure(bool b, const std::string& s) {
      if (not b)
	throw std::runtime_error("ocean "+s+" error");
    }
  };
}
}

#endif
