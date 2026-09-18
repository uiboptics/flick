#ifndef flick_material_pure_water_refractive_index
#define flick_material_pure_water_refractive_index

#include "../../../environment/input_output.hpp"
#include "../../../numeric/table.hpp"
#include "../../../numeric/constants.hpp"

namespace flick {
namespace material {
namespace water { 
  class refractive_index
  // Real refractive index of pure sea water
  {
    mutable double S_; // [PSU]
    mutable double T_; // [Kelvin]
    mutable double Tc_; // [Celsius]
    pp_table wang_etal_;
    std::vector<double> a = {1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6,
      15.868, 0.01155, -0.00423, -4382, 1.1455e6};
    const std::string path_{"/material/water/refractive_index"};
  public:
    refractive_index(double salinity, double temperature) {
      set_S_T(salinity, temperature);
      wang_etal_ = read<pp_table>(path_+"/wang_etal_real.txt");    
    }
    void set_S_T(double S, double T) const {
      S_ = S;
      T_ = T;
      Tc_ = constants::to_celsius(T_);
    }
    double value(double wavelength) const {
      double wl_min = 280e-9;
      double wl_max = 700e-9;
      double wl = wavelength;
      if (wl < wl_min) {
	return wang_etal(wavelength) * salinity_factor(wl_min);
      } else if (wl > wl_max) {
	return wang_etal(wavelength) * salinity_factor(wl_max);
      } else {
	return quan_and_fry(wavelength);
      }
    } 
    double dn_dS(double wavelength) const {
      double wl = wavelength*1e9; // [nm]
      return (a[1]+a[2]*Tc_+a[3]*pow(Tc_,2)+a[6]/wl)*
	air_refractive_index(wavelength);
    }  
    double density_variation(double wavelength) const {
      // Proutiere, A., Megnassan, E. and Hucteau, H.,
      // 1992. Refractive index and density variations in pure
      // liquids: A new theoretical relation. The Journal of Physical
      // Chemistry, 96(8), pp.3485-3489.
      double n = value(wavelength);
      double n2 = pow(n,2);
      return (n2-1)*(1+2./3*(n2+2)*pow(n/3-1/(3*n),2));
    }   
  private:
    double salinity_factor(double wavelength) const {
      return quan_and_fry(wavelength) / wang_etal(wavelength);
    }
    double wang_etal(double wavelength) const {
      double wl_um = wavelength*1e6;
      return wang_etal_.value(wl_um, T_);
    }
    double quan_and_fry(double wavelength) const {
      // Quan, X. and Fry, E.S., 1995. Empirical equation for the
      // index of refraction of seawater. Applied optics, 34(18),
      // pp.3477-3480.
      double lambda = wavelength*1e9; // [nm]
      return (a[0]+(a[1]+a[2]*Tc_+a[3]*pow(Tc_,2))*S_+a[4]*pow(Tc_,2)+
	      (a[5]+a[6]*S_+a[7]*Tc_)/lambda+a[8]/pow(lambda,2)+a[9]/pow(lambda,3))*
	air_refractive_index(wavelength);
    }
    double air_refractive_index(double wavelength) const
    // Ciddor, P.E., 1996. Refractive index of air: new equations for
    // the visible and near infrared. Applied optics, 35(9),
    // pp.1566-1573.
    {
      std::vector<double> k = {238.0185, 5792105, 57.362, 167917}; // [microns^-2]
      double nu = 1/(wavelength*1e6); // [microns^-1]
      return 1+(k[1]/(k[0]-pow(nu,2))+k[3]/(k[2]-pow(nu,2)))/1e8;
    }      
  };
}
}
}  

#endif
