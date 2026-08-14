#ifndef flick_polydispersed_mie
#define flick_polydispersed_mie

#include "basic_monodispersed_mie.hpp"
#include "../numeric/legendre/legendre.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/physics_function.hpp"
#include "../environment/input_output.hpp"
#include <tuple>

namespace flick {
  class basic_quantity {
  protected:
    basic_monodispersed_mie& bm_;
    const size_distribution& sd_;
    stdvector center_quantity_;
    double alpha_{0};
    size_t size_{1};
    bool do_center_subtraction_ = false;
  public:
    basic_quantity(basic_monodispersed_mie& bm,
		   const size_distribution& sd)
      : bm_{bm}, sd_{sd} {
    }
    virtual ~basic_quantity() = default;
    virtual stdvector value(double x) const = 0;

    stdvector center_quantity() const {
      return center_quantity_;
    }
    stdvector center_approximation_total() const {
      return center_quantity_ * sd_.weighted_integral(alpha_);
    }
    stdvector center_subtraction_total() const {
      if (do_center_subtraction_)
	return center_approximation_total();
      return stdvector(size(),0);
    }
    const size_distribution& sd() const {
      return sd_;
    }
    double alpha() const {
      return alpha_;
    }
    size_t size() const {
      return center_quantity_.size();
    }
    stdvector transformed_value(const stdvector& quantity) const {
      double r = bm_.radius();
      if (do_center_subtraction_)
	return (quantity - center_quantity_ * pow(r,alpha_)) * r * sd_.value(r);
      return quantity * r * sd_.value(r);
    }
    stdvector transformed_center(const stdvector& quantity) {
      return quantity/pow(bm_.radius(),alpha_);
    }
  };
 
  struct absorption_quantity : public basic_quantity {
    absorption_quantity(basic_monodispersed_mie& bm,
			     const size_distribution& sd)
      : basic_quantity(bm,sd) {
      alpha_ = 3;
      bm_.radius(sd_.center());
      center_quantity_ = transformed_center({bm_.absorption_cross_section()});
    }  
    stdvector value(double x) const {
      bm_.radius(exp(x));
      return transformed_value({bm_.absorption_cross_section()});
    }
  };

  struct scattering_quantity : public basic_quantity {
    scattering_quantity(basic_monodispersed_mie& bm,
			const size_distribution& sd)
      : basic_quantity(bm,sd) {
      alpha_ = 2;
      bm_.radius(sd_.center());
      center_quantity_ = transformed_center({bm_.scattering_cross_section()});
    }  
    stdvector value(double x) const {
      bm_.radius(exp(x));
      return transformed_value({bm_.scattering_cross_section()});
    }
  };

  struct smatrix_quantity : public basic_quantity {
    size_t row_;
    size_t col_;
    smatrix_quantity(basic_monodispersed_mie& bm,
		     const size_distribution& sd,
		     size_t row, size_t col)
      : basic_quantity(bm,sd), row_{row}, col_{col} {
      alpha_ = 2;
      bm_.radius(sd_.center());
      center_quantity_ = transformed_center({bm_.scattering_matrix_element(row_,col_)});
    }  
    stdvector value(double x) const {
      bm_.radius(exp(x));
      return transformed_value(bm_.scattering_matrix_element(row_,col_));
    }
  };
 
  template<class Monodispersed_mie, class Size_distribution>
  class polydispersed_mie {
    const double epsilon_ = std::numeric_limits<double>::epsilon();
    Monodispersed_mie mm_;
    const Size_distribution sd_;
    double accuracy_{0.05};
    bool keep_integration_points_ = true;
    pl_function xy_points_;

    stdvector integral(const basic_quantity& bq) {
      double max_step_factor = 0.25;
      double step_factor = max_step_factor;
      double x0 = log(bq.sd().center());
      accumulated_integral_vector ai(bq, 100*accuracy_);
      if (bq.sd().width() > 0.5 && accuracy_ >= 0.05)
	return bq.center_approximation_total();
      stdvector a0 = bq.center_subtraction_total();
      ai.set_total(a0);
      ai.keep_integration_points(keep_integration_points_);
      double width = 1.2*bq.sd().width();
      double x1, x2;
      double dx = step_factor * width;
      for (size_t i=0; i<2; i++) {
	x1 = x0;
	ai.reset_convergence();
	while(ai.significant_added_value()) {
	  if (i==0) {
	    x2 = x1 - dx;
	    ai.add_value(x2,x1,width);
	  } else {
	    x2 = x1 + dx;
	    ai.add_value(x1,x2,width);
	  }
	  if (not ai.has_converged()) {
	    step_factor /= 2;
	    ai.set_total(ai.previous_total());
	  } else {
	    x1 = x2;
	  }
	}	
      }
      if (keep_integration_points_)
	xy_points_ = ai.integration_points();
      return ai.total();
    }
  public:
    polydispersed_mie(const Monodispersed_mie& mm, const Size_distribution& sd)
      : mm_{mm}, sd_{sd} {
    }
    void percentage_accuracy(double p) {
      accuracy_ = p/100;
    }
    void angles(const stdvector& angles) {
      mm_.angles(angles);
    }
    void set_wavelength(double wl) {
      mm_.set_wavelength(wl);
    }
    void set_refractive_indices(stdcomplex m_host, stdcomplex m_sphere) {
      mm_.set_refractive_indices(m_host, m_sphere);
    }
    pl_function xy_points() {
      return xy_points_;
    }
    double absorption_cross_section() {
      if (sd_.width() < epsilon_) {
	mm_.radius(sd_.center());
	return mm_.absorption_cross_section();
      } else {
	return integral(absorption_quantity(mm_,sd_))[0];
      }
    }
    double scattering_cross_section() {
      if (sd_.width() < epsilon_) {
	mm_.radius(sd_.center());
	return mm_.scattering_cross_section();
      } else {
	return integral(scattering_quantity(mm_,sd_))[0];
      }
    }
    stdvector scattering_matrix_element(size_t row, size_t col) {
      if (sd_.width() < epsilon_) {
	mm_.radius(sd_.center());
	return mm_.scattering_matrix_element(row,col);
      } else {
	return integral(smatrix_quantity(mm_,sd_,row,col)); 
      }
    }
    stdvector normalized_scattering_matrix_element(size_t row, size_t col)
    /* Note that the normalization is such that the upper left corner
       element is normalized to 4*pi, not one, which is normally the
       flick convention for the phase function */
    {
      double k = 4*constants::pi/scattering_cross_section();
      return k*scattering_matrix_element(row,col);
    }
    std::tuple<std::vector<stdvector>, std::vector<stdvector>, stdvector> ab_functions() {
      std::vector<stdvector> a(4);
      a[0] = normalized_scattering_matrix_element(0,0);
      a[1] = normalized_scattering_matrix_element(1,1);
      a[2] = normalized_scattering_matrix_element(2,2);
      a[3] = normalized_scattering_matrix_element(3,3);
      std::vector<stdvector> b(2);
      b[0] = normalized_scattering_matrix_element(0,1);
      b[1] = normalized_scattering_matrix_element(2,3);
      stdvector x = vec::cos(mm_.angles());
      std::reverse(x.begin(),x.end());
      for (size_t i=0; i<a.size(); ++i)
	std::reverse(a[i].begin(),a[i].end());
      for (size_t i=0; i<b.size(); ++i)
	std::reverse(b[i].begin(),b[i].end());
      return {a,b,x};
    }    
    double scattering_efficiency() {
      return scattering_cross_section()/sd_.average_area();
    }
    double absorption_efficiency() {
      return absorption_cross_section()/sd_.average_area();
    }
  };
}

#endif
