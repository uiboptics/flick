#ifndef flick_layered_iops
#define flick_layered_iops

#include "material.hpp"
#include "ab_functions.hpp"

namespace flick {
  class layered_iops {
    std::shared_ptr<material::base> m_;
    stdvector boundaries_;
    size_t n_terms_;
    stdvector oda_;
    stdvector ods_;
    std::vector<std::vector<stdvector>> alpha_;
    std::vector<std::vector<stdvector>> beta_;
    stdvector refidx_;
    stdvector s_height_;
  public:
    layered_iops(std::shared_ptr<material::base> m,
		 const stdvector& boundaries, size_t n_terms)
      : m_{m}, boundaries_{boundaries},
	n_terms_{n_terms},
	oda_(n_layers()),
	ods_(n_layers()),
	alpha_(4, std::vector<stdvector>(n_layers(), stdvector(n_terms))),
	beta_(2, std::vector<stdvector>(n_layers(), stdvector(n_terms))),
	refidx_(n_layers()), s_height_(n_layers())
    {
      if (boundaries.size() < 2 or boundaries[1] < boundaries[0] or
	  (not std::is_sorted(boundaries.begin(), boundaries.end())))
	throw std::runtime_error("boundary error in layered_iops");
      if (n_terms < 3)
	throw std::runtime_error("number of terms error in layered_iops");
    }
    void set_wavelength(double wl) {
      m_->set_wavelength(wl);
      update();
    }
    size_t n_layers() const {
      return boundaries_.size() - 1;
    }
    stdvector scattering_optical_depth() const {
      return ods_;
    }
    stdvector absorption_optical_depth() const {
      return oda_;
    }
    stdvector single_scattering_albedo() const {
      return (ods_ + oda_) / ods_;
    }
    std::vector<stdvector> alpha_terms(size_t n) const {
      return alpha_.at(n);
    }
    std::vector<stdvector> beta_terms(size_t n) const {
      return beta_.at(n);
    }
    stdvector refractive_index() const {
      return refidx_;
    }
    stdvector absorption_coefficient() const {
      return oda_ / layer_thicknesses();
    }
    stdvector scattering_coefficient() const {
      return ods_ / layer_thicknesses();
    }
    double delta_fit_scaling_factor(size_t layer_no) const {
      return 4*std::numbers::pi*alpha_terms(0)[layer_no][0];
    }
    double asymmetry_factor(size_t layer_no) const {
      return 4*std::numbers::pi/3*alpha_terms(0)[layer_no][1];
    }
    double phase_function(size_t layer_no, double angle) const {
      double h = m_->pose().position().z();
      m_->set_position({0,0,s_height_.at(layer_no)});
      double p = material::phase_function(*m_).value(cos(angle));
      m_->set_position({0,0,h});
      return p;
    }

    friend std::ostream& operator<<(std::ostream &os,
				    const layered_iops& iops) {
      auto h = iops.boundaries_;
      auto oda = iops.absorption_optical_depth();
      auto ods = iops.scattering_optical_depth();
      auto real_n = iops.refractive_index();
      os << std::defaultfloat << std::setprecision(4) << "\n";
      os << "Wavelength [m]: " << iops.m_->wavelength() << "\n";
      os << "column 1: Layer bottom boundary height [m]" << "\n";
      os << "column 2: Layer geometrical thickness [m]" << "\n";
      os << "column 3: Absorption optical thickness" << "\n";
      os << "column 4: Scattering optical thickness without delta-fit" << "\n";
      os << "column 5: Real refractive index" << "\n";
      os << "column 6: Delta-fit scattering scaling factor" << "\n";
      os << "column 7: Asymmetry factor with delta-fit scaling" << "\n";
      os << "column 8: Average layer scattering height [m]" << "\n";
      os << "column 9: Phase function forward at scattering height, without delta-fit [1/sr]" << "\n";
      os << "column 10: Phase function at pi/2 at scattering height without delta-fit [1/sr]" << "\n";
      os << "column 11: Phase function backward at scattering height without delta-fit [1/sr]" << "\n";
      size_t n_last = oda.size()-1;
      for (size_t i = 0; i<oda.size(); i++) {
	size_t n = n_last - i;
	double dh = h[n_last];
	if (i > 0)
	  dh = h[n+1]-h[n];
	os << std::scientific << std::setprecision(3) 
	   << "(1)"<< h[n] << " (2)" << dh << " (3)"<< oda[n] << " (4)"
	   << ods[n] << " (5)"
	   << real_n[n] << " (6)"
	   << iops.delta_fit_scaling_factor(n)<< " (7)"
	   << iops.asymmetry_factor(n) << " (8)"
	   << iops.s_height_[n] << " (9)"
	   << iops.phase_function(n, 0) << " (10)"
	   << iops.phase_function(n, std::numbers::pi/2) << " (11)"
	   << iops.phase_function(n, std::numbers::pi)
	   << std::endl;
      }
      return os;
    }  
  private:
    void update() {
      m_->set_direction({0,0});
      m_->set_position({0,0,boundaries_[0]});
      for (size_t i=0; i < n_layers(); i++) {	
	double h = layer_thickness(i);
	s_height_[i] = average_scattering_height(boundaries_[i],boundaries_[i+1]);
	double dh = s_height_[i] - boundaries_[i];
	move(dh);
	set_alpha_beta(i);
	refidx_[i] = m_->real_refractive_index();
	move(-dh);
	set_optical_depth(i);
	move(h);
      }
    }
    double average_scattering_height(double h_low, double h_high) {
      stdvector h = range(h_low,h_high,100).linspace();
      stdvector s(h.size());
      stdvector sh(h.size());
      vector p0 = m_->pose().position();
      for (size_t i=0; i<sh.size(); i++) {
	m_->set_position({0,0,h[i]});
	s[i] = m_->scattering_coefficient();
	sh[i] = s[i]*h[i];
      }
      m_->set_position(p0);
      pl_function shf(h,sh);
      pl_function sf(h,s);
      double h_avg = shf.integral()/sf.integral();
      if (std::isfinite(h_avg))
      	return h_avg;
      return h_low + (h_high-h_low)/2;
    }
    void set_alpha_beta(size_t i) {
      auto [alpha, beta] = material::fitted_mueller_alpha_beta(*m_,n_terms_);
      for (size_t n=0; n < alpha_.size(); n++) {
	alpha_[n][i] = alpha[n];
      }
      for (size_t n=0; n < beta_.size(); n++)
	beta_[n][i] = beta[n];
    }
    void set_optical_depth(size_t i) {
      oda_[i] = m_->absorption_optical_depth(layer_thickness(i));
      ods_[i] = m_->scattering_optical_depth(layer_thickness(i));      
    }
    double layer_thickness(size_t i) const {
      return boundaries_.at(i+1)-boundaries_.at(i);
    }
    void move(double distance) const {
      vector new_position = m_->pose().position() +
	m_->pose().direction()*distance;
      m_->set_position(new_position);
    }
    stdvector layer_thicknesses() const {
      stdvector t(n_layers());
      for (size_t i=0; i < t.size(); i++)
	t[i] = layer_thickness(i);
      return t;
    }
  };
}

#endif
