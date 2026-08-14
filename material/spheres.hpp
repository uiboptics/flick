#ifndef flick_material_spheres
#define flick_material_spheres

#include "material.hpp"
#include "water/pure_water.hpp"
#include "ice/pure_ice.hpp"
#include "../mie/polydispersed_mie.hpp"
#include "../mie/monodispersed_mie.hpp"
#include "../mie/parameterized_monodispersed_mie.hpp"
#include "../polarization/mueller.hpp"

namespace flick {
namespace material {
  template<class Size_distribution, class Host_material,
	   class Sphere_material, class Monodispersed_mie>
  class spheres : public base {
    double volume_fraction_;
    Size_distribution size_distribution_;
    Host_material host_material_;
    Sphere_material sphere_material_;
    const std::vector<int> row_{0,0,1,1,2,2,3,3};
    const std::vector<int> col_{0,1,0,1,2,3,2,3};
    mutable std::vector<pl_function> scattering_matrix_elements_;
    mutable Monodispersed_mie mono_mie_;
    mutable polydispersed_mie<Monodispersed_mie,Size_distribution> poly_mie_;
    mutable bool has_changed_{true};
  public:
    spheres(const spheres&) = delete;
    spheres& operator=(const spheres&) = delete;
    spheres(double volume_fraction,
	    const Size_distribution& sd,
	    const Host_material& hm,
	    const Sphere_material& sm)
      : volume_fraction_{volume_fraction}, size_distribution_{sd},
	host_material_{hm}, sphere_material_{sm},
	scattering_matrix_elements_(row_.size()),
	mono_mie_{host_material_.refractive_index(),
      sphere_material_.refractive_index(), wavelength()},
	poly_mie_{mono_mie_, size_distribution_}
    {
    }
    void set_wavelength(double wl) override {
      base::set_wavelength(wl);
      host_material_.set_wavelength(wl);
      sphere_material_.set_wavelength(wl);
      mono_mie_.set_wavelength(wl);
      mono_mie_.set_refractive_indices(host_material_.refractive_index(),
				       sphere_material_.refractive_index());
      poly_mie_.set_wavelength(wl);
      poly_mie_.set_refractive_indices(host_material_.refractive_index(),
				       sphere_material_.refractive_index());
      has_changed_ = true;
    }
    void set_angles(const stdvector& angles) override {
      base::set_angles(angles);
      mono_mie_.angles(angles);
      poly_mie_.angles(angles);
      has_changed_ = true;
    }
    void percentage_accuracy(double p) {
      poly_mie_.percentage_accuracy(p);
    }
    double absorption_coefficient() const override {
      return poly_mie_.absorption_cross_section()
	* size_distribution_.particles_per_volume(volume_fraction_);
    }
    double scattering_coefficient() const override {
      return poly_mie_.scattering_cross_section()
	      * size_distribution_.particles_per_volume(volume_fraction_);
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      if (has_changed_) {
	for (size_t i=0; i < row_.size(); ++i) {
	  scattering_matrix_elements_[i] =
	    pl_function(mono_mie_.angles(),
			poly_mie_.scattering_matrix_element(row_[i],col_[i]));
	}
	has_changed_ = false;
      }
      double theta = angle(scattering_direction);
      mueller m;
      double c = poly_mie_.scattering_cross_section();
      for (size_t i=0; i<scattering_matrix_elements_.size(); ++i) {
	double s = scattering_matrix_elements_[i].value(theta);
	m.add(row_[i],col_[i],s/c);
      }
      return m;  
    }
    double real_refractive_index() const override {
      return 1;
    }   
  };

  template<class Monodispersed_mie>
  struct bubbles_in_ice : public spheres<log_normal_distribution,
					 material::pure_ice,
					 material::vacuum,
					 Monodispersed_mie> {
    bubbles_in_ice(double volume_fraction,double mu, double sigma) :
      spheres<log_normal_distribution,
	      material::pure_ice,
	      material::vacuum,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::pure_ice(),
				 material::vacuum()) {}
  };

  template<class Monodispersed_mie>
  struct bubbles_in_water : public spheres<log_normal_distribution,
					 material::pure_water,
					 material::vacuum,
					 Monodispersed_mie> {
    bubbles_in_water(double volume_fraction,double mu, double sigma,
		     double salinity, double temperature) :
      spheres<log_normal_distribution,
	      material::pure_water,
	      material::vacuum,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::pure_water(salinity,temperature),
				 material::vacuum()) {}
  };

  template<class Monodispersed_mie>
  struct brines_in_ice : public spheres<log_normal_distribution,
					 material::pure_ice,
					 material::pure_water,
					 Monodispersed_mie> {
    brines_in_ice(double volume_fraction, double mu, double sigma,
		   double salinity) :
      spheres<log_normal_distribution,
	      material::pure_ice,
	      material::pure_water,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::pure_ice(),
				 material::pure_water(salinity,273.15)) {}
  };
  
  template<class Monodispersed_mie>
  struct water_cloud : public spheres<log_normal_distribution,
				      material::vacuum,
				      material::pure_water,
				      Monodispersed_mie> {
    water_cloud(double volume_fraction, double mu, double sigma) :
      spheres<log_normal_distribution,
	      material::vacuum,
	      material::pure_water,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::vacuum(),
				 material::pure_water()) {}
  };
  
  template<class Monodispersed_mie>
  struct ice_cloud : public spheres<log_normal_distribution,
				      material::vacuum,
				      material::pure_ice,
				      Monodispersed_mie> {
    ice_cloud(double volume_fraction, double mu, double sigma) :
      spheres<log_normal_distribution,
	      material::vacuum,
	      material::pure_ice,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::vacuum(),
				 material::pure_ice()) {}
  }; 
}
}

#endif
