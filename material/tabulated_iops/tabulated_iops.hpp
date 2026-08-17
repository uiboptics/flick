#ifndef flick_material_tabulated_iops
#define flick_material_tabulated_iops

#include "../material.hpp"
#include "../../numeric/table.hpp"
#include "../../numeric/flist.hpp"
#include "../../numeric/function.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class tabulated_iops : public base {
    double volume_fraction_;
    pl_flist ab_;
    pe_table s_11_;
    const double to_nm_{1e9};
  public:
    tabulated_iops(const std::string& name, double volume_fraction)
      : volume_fraction_{volume_fraction} {
      std::string p = "/material/tabulated_iops/iop_tables/"+name+'/';
      ab_ = read<pl_flist>("ab.txt", p);
      matrix<size_t> m = read<matrix<size_t>>("include_elements.txt", p);
      for (size_t i=0; i < m.n_rows(); i++) {
	std::string fname = "s_"+std::to_string(m.element(i,0))+
	  std::to_string(m.element(i,1))+".txt";
	s_11_ = read<pe_table>(fname, p);
      //a_.add_constant_extrapolation();
      //b_.add_constant_extrapolation();
      //m11_.add_constant_extrapolation();
      }
    }
    double absorption_coefficient() const override {
      return ab_(0).value(wavelength()*to_nm_)*volume_fraction_;
    }
    double scattering_coefficient() const override {
      return ab_(1).value(wavelength()*to_nm_)*volume_fraction_;
    }
  private:   
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      mueller m;
      double cos_ang = std::clamp<double>(scattering_direction.mu(),-1,1);
      double ang_deg = 180/constants::pi * acos(cos_ang);
      return m.add(0,0,s_11_.value(ang_deg, wavelength()));
    }
  };
}
}

#endif
