#ifndef flick_mueller
#define flick_mueller
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/std_operators.hpp"
#include <numbers>
#include <algorithm>
#include <vector>

namespace flick {
  class mueller
  // Holds non-zero Mueller matrix element values for one angle
  {
    struct element {
      size_t row;
      size_t col;
      double value;
    };
    std::vector<element> elements_;
  public:
     mueller& add(size_t row, size_t col, double value) {
      elements_.push_back(element{row,col,value});
      return *this;
    }
    size_t size() const {
      return elements_.size();
    }
    const element& operator()(size_t n) const {
      return elements_.at(n);
    }
    double value(size_t row, size_t col) const {
      auto& e = elements_;
      for (size_t n = 0; n < e.size(); ++n) {
	if (e[n].row == row && e[n].col == col) {
	  return e[n].value;
	}
      }
      return 0;
    }
  };

  std::ostream& operator<<(std::ostream &os, const mueller& m) {
    size_t n = 0;
    os << '\n';
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
	if (n < m.size() && m(n).row == i && m(n).col == j) {
	  os << m(n).value << " ";
	  n++;
	} else {
	  os << 0 << " ";
	}
      }
      os << '\n';
    }
    return os;
  }

  std::vector<double> hg_importance_sampling(double g, size_t n_points) {
    henyey_greenstein hg(g);
    std::vector<double> fraction = range(0,1,n_points).linspace();
    std::vector<double> x(n_points); 
    for (size_t i=0; i<n_points; ++i) {
      x[i] = hg.inverted_accumulated_angle(fraction[i]); 
    }
    return x;
  }

  template<class Function>
  Function resample(const Function& f, const std::vector<double>& new_x) {
    std::vector<double> new_y(new_x.size());
    for (size_t i=0; i<new_x.size(); ++i) {
      new_y[i] = f.value(new_x[i]);
    }
    return Function{new_x, new_y};
  }

  template<class Function>
  Function to_cos_x(const Function& f) {
    Function fc;
    size_t n = f.size();
    for (size_t i=0; i<n; ++i) {
      fc.append({cos(f.x()[n-i-1]), f.y()[n-i-1]});
    }
    return fc;
  }
  
  template<class Function>
  Function to_acos_x(const Function& f) {
    Function fc;
    size_t n = f.size();
    for (size_t i=0; i<n; ++i) {
      fc.append({acos(std::clamp<double>(f.x()[n-i-1],-1,1)), f.y()[n-i-1]});
    }
    return fc;
  }
  
  template<class Function>
  Function hg_phase_function(double g, size_t n_points,
			     double sampling_skew_factor = 0.8) {
    henyey_greenstein hg(g);
    std::vector<double> x = hg_importance_sampling(sampling_skew_factor, n_points);
    x = vec::cos(x);
    std::vector<double> y(x.size());
    for (size_t i=0; i<x.size(); ++i) {
      y[i] = hg.value(x[i]);
    }
    return Function{x,y};    
  }

  class angle_converter
  // Converts back and forth between mu = cos(angle) and forward
  // streached half-angles (fsha)
  {
    double q_;
  public:
    explicit angle_converter(double q) : q_{q} {
    }
    double to_fsha(double mu) const {
      return std::pow((1 - mu) / 2, 1 / q_);
    }
    double to_mu(double fsha) const {
      return 1 - 2 * std::pow(fsha, q_);
    }
    std::vector<double> to_fsha(const std::vector<double>& mu_values) const {
      std::vector<double> fsha_values;
      fsha_values.reserve(mu_values.size());
      for (double mu : mu_values) {
	fsha_values.push_back(to_fsha(mu));
      }
      std::reverse(fsha_values.begin(), fsha_values.end());
      return fsha_values;
    }
    std::vector<double> to_mu(const std::vector<double>& fsha_values) const {
      std::vector<double> mu_values;
      mu_values.reserve(fsha_values.size());
      for (double fsha : fsha_values) {
	mu_values.push_back(to_mu(fsha));
      }
      std::reverse(mu_values.begin(), mu_values.end());
      return mu_values;
    }
    static std::vector<double> reversed(const std::vector<double>& values) {
      return {values.rbegin(), values.rend()};
    }
  };

  class tabulated_phase_function {
    pl_function p_;
    double q_;
    pe_function p_mu_;
    static constexpr double pi_ = std::numbers::pi;
  public:
    tabulated_phase_function() = default;
    tabulated_phase_function(const pe_function& p_mu, double q=4)
      : q_{q}, p_mu_{p_mu} {
      std::vector<double> x = angle_converter(q).to_fsha(p_mu.x());
      std::vector<double> y = angle_converter::reversed(p_mu.y());
      if (x.front() > 0) {
	x.insert(x.begin(), 0);
	y.insert(y.begin(), y.front());
      }
      if (x.back() < 1) {
	x.push_back(1);
	y.push_back(y.back());
      }
      p_ = pl_function(x,y);
    }
    void normalize() {
      double scale = 1/integral_4pi();
      p_.scale_y(scale);
      p_mu_.scale_y(scale);
    }
    const std::vector<double>& x_mu() const {
      return p_mu_.x();
    }    
    const std::vector<double>& y_mu() const {
      return p_mu_.y();
    }
    bool empty() const {
      return p_.size() == 0;
    }
    double value(double mu) const {
      double fsha = angle_converter{q_}.to_fsha(mu);
      return p_.value(fsha);
    }
    double asymmetry_factor() const {
      double deficit = 8*pi_*q_*weighted_integral(2*q_-1);
      return 1-deficit/integral_4pi();
    }
    double integral_4pi() const {
      return 4*pi_*q_*weighted_integral(q_-1);
    }  
  private:
    double weighted_integral(double power) const {
      double integral = 0;
      const std::vector<double>& x = p_.x();
      const std::vector<double>& y = p_.y();
      for (size_t i=0; i+1<x.size(); ++i) {
	double x1 = x[i];
	double x2 = x[i+1];
	double slope = (y[i+1]-y[i])/(x2-x1);
	double intercept = y[i]-slope*x1;
	integral += slope/(power+2)*
	  (pow(x2,power+2)-pow(x1,power+2));
	integral += intercept/(power+1)*
	  (pow(x2,power+1)-pow(x1,power+1));
      }
      return integral;
    }
  };
   
  class angular_mueller
  // Keeps Mueller matrix values at discrete angles.
  // Elements are kept as values divided by phase function
  {
    tabulated_phase_function p_;
    struct element {
      size_t row;
      size_t col;
      pl_function f;
    };
    std::vector<element> elements_;
    double epsilon = 1e-6;
  public:
    angular_mueller() = default;
    angular_mueller(const tabulated_phase_function& p)
      : p_{p} {
      add(0,0,pl_function{p.x_mu(),p.y_mu()});
    }
    angular_mueller& append(size_t row, size_t col,
			    double mu, double value) {
      ensure(mu > -1-epsilon && mu < 1+epsilon);
      auto& e = elements_;
      for (size_t i=0; i<e.size(); ++i) {
	if (e[i].row == row && e[i].col == col) {
	  e[i].f.append({mu,value});
	  return *this;
	}
      }
      element new_e;
      new_e.row = row;
      new_e.col = col;
      new_e.f = pl_function{};
      e.push_back(new_e);
      e.back().f.append({mu,value});
      return *this;
    }
    angular_mueller& add(size_t row, size_t col,
			 const pl_function& f) {//should be tabulated phf
      elements_.push_back(element{row,col, normalize(f)});
      return *this;
    }
    size_t size() const {
      return elements_.size();
    }
    const element& operator()(size_t n) const {
      return elements_.at(n);
    }
    double value(size_t row, size_t col, double mu) const {
      ensure(mu > -1-epsilon && mu < 1+epsilon);
      auto& e = elements_;
      for (size_t n = 0; n < e.size(); ++n) {
	if (e[n].row == row && e[n].col == col) {
	  return e[n].f.value(mu)*p_.value(mu);
	}
      }
      return 0;
    }
    angular_mueller& add(const angular_mueller& am, double fractional_weight) {
      if (p_.empty()) {
	*this = am; 
	return *this;
      }
      double w = fractional_weight;
      auto& e1 = elements_;
      auto& e2 = am.elements_;
      std::vector<element> new_e;
      std::vector<bool> added1(e1.size(),false);
      std::vector<bool> added2(e2.size(),false);
      for (size_t i=0; i<e1.size(); ++i) {
	for (size_t j=0; j<e2.size(); ++j) {
	  if (e1[i].row == e2[j].row && e1[i].col == e2[j].col) {
	    added1[i]=true;
	    added2[j]=true;
	    element e;
	    e.row = e1[i].row;
	    e.col = e1[i].col;
	    e.f = add_elements(e1[i].f, p_, e2[j].f, am.p_, w);
	    new_e.push_back(e);
	  }
	}
      }
      for (size_t i=0; i<e1.size(); ++i) {
	if (!added1[i]) {
	  element e;
	  e.row = e1[i].row;
	  e.col = e1[i].col;
	  e.f = add_elements(e1[i].f, p_, {0}, am.p_, w);
	  new_e.push_back(e);
	}
      }
      for (size_t i=0; i<e2.size(); ++i) {
	if (!added2[i]) {
	  element e;
	  e.row = e2[i].row;
	  e.col = e2[i].col;
	  e.f = add_elements({0}, p_, e2[i].f, am.p_, w);
	  new_e.push_back(e);
	}	  
      }
      elements_ = new_e;
      update_phase_function(am, w);
      return *this;
    }
    void print(double mu) {
      std::cout << std::endl;
      for (size_t i = 0; i < 4; ++i) {
	for (size_t j = 0; j < 4; ++j) {
	  std::cout << value(i,j,mu) << " ";
	}
	std::cout << std::endl;
      }
    }
  private:
    pl_function add_elements(const pl_function& f1,
			     const tabulated_phase_function& p1,
			     const pl_function& f2,
			     const tabulated_phase_function& p2,
			     double w) {
      std::vector<double> x = p_.x_mu();
      std::vector<double> y(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	double mixed_p = (1-w)*p1.value(x[i]) + w*p2.value(x[i]);
	y[i] = ((1-w)*f1.value(x[i])*p1.value(x[i]) +
		w*f2.value(x[i])*p2.value(x[i])) / mixed_p;
      }
      return pl_function{x,y};
    }
    void update_phase_function(const angular_mueller& am,
				double fractional_weight) {
      double w = fractional_weight;
      pe_function p = pe_function{p_.x_mu(),p_.y_mu()};
      pe_function amp = pe_function{am.p_.x_mu(),am.p_.y_mu()};
      pl_function p1{p.x(),p.y()};
      pl_function p2{amp.x(),amp.y()};
      std::vector<double> x = p1.x();
      if (p2.size() > p1.size())
	x = p2.x();
      std::vector<double> y(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	y[i] = (1-w)*p1.value(x[i]) + w*p2.value(x[i]);
      }      
      p_ = tabulated_phase_function(pe_function{x,y});
    }
    pl_function normalize(const pl_function& f) {
      pl_function f_new;
      for (size_t i=0; i < f.size(); ++i) {
	double x = f.x()[i];
	ensure(x >= -1-epsilon && x <= 1+epsilon);
	double y = f.y()[i] / p_.value(x);
	f_new.append({x,y});
      }
      return f_new;
    }
    void ensure(bool b) const {
      if (!b)
	throw std::runtime_error("angular_mueller");
    }
  };
}

#endif
