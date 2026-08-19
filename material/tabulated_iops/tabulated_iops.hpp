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
    std::vector<pe_table> s_;
    matrix<size_t> e_;
    const double to_nm_{1e9};
  public:
    tabulated_iops(const std::string& name, double volume_fraction)
      : volume_fraction_{volume_fraction} {
      std::string p = "/material/tabulated_iops/iop_tables/"+name+'/';
      ab_ = read<pl_flist>("ab.txt", p);
      e_ = read<matrix<size_t>>("include_elements.txt", p);
      s_.resize(e_.n_rows());
      for (size_t i=0; i < s_.size(); i++) {
	std::string fname = "s_"+std::to_string(e_.element(i,0))+
	  std::to_string(e_.element(i,1))+".txt";
	s_[i] = read<pe_table>(fname, p);
      }
    }
    double absorption_coefficient() const override {
      return ab_(0).value(wavelength()*to_nm_)*volume_fraction_;
    }
    double scattering_coefficient() const override {
      return ab_(1).value(wavelength()*to_nm_)*volume_fraction_;
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      mueller m;
      double cos_ang = std::clamp<double>(scattering_direction.mu(),-1,1);
      double ang_deg = 180/constants::pi * acos(cos_ang);
      for (size_t i=0; i < s_.size(); ++i)
	m.add(e_.element(i,0)-1,e_.element(i,1)-1,
	      s_[i].value(wavelength()*to_nm_, ang_deg));
      return m;
    }
  };
}
}


#endif
