#ifndef flick_wigner_fit
#define flick_wigner_fit

#include "../std_operators.hpp"
#include "../function.hpp"
#include "../constants.hpp"
#include "../linalg/solve_with_eigen.hpp"
#include "wigner_d.hpp"
#include <numbers>
 
namespace flick {
  enum class fit {
    absolute,
    relative,
    scaling,
    relative_and_cutoff
  };
  struct default_scaling_function {
    double value() {
      return 1;
    }
  };
  size_t wigner_n_sampling_points(size_t n_terms) {
    return pow(n_terms, 1.6);
  }
  double wigner_theta_min(size_t n_angles) {
    return 1./n_angles;
  }
  stdvector wigner_x_uniform_angle(size_t n_angles) {
    double theta_min = wigner_theta_min(n_angles);
    stdvector theta = range(theta_min, std::numbers::pi, n_angles).linspace();
    return -1*vec::cos(theta);
  }
  stdvector wigner_uniform_x(size_t n_angles) {
    double x_max = cos(wigner_theta_min(n_angles));
    return range(-1, x_max, n_angles).linspace();
  }
  stdvector wigner_x_values(size_t n_angles) {
    return wigner_uniform_x(n_angles);
  }  
  template<class Function>
  class wigner_fit {
    stdvector coefficients_;
    int m_, n_;
    double x_cutoff_;
  public:
    wigner_fit(const Function& f, int m, int n, int n_terms, fit fit,
	       const pe_function& scaling_function = pe_function{{-1,1},{1,1}})
      : m_{m}, n_{n}, coefficients_(n_terms) {
      x_cutoff_ = cos(3./180*constants::pi);
      size_t n_angles = wigner_n_sampling_points(n_terms);
      stdvector x = wigner_x_values(n_angles);
      linalg::matrix mat(x.size(), std::vector<double>(n_terms));
      linalg::vector v(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	stdvector d = wigner_d(x[i], m_, n_, n_terms).terms();
	if (fit == fit::absolute)
	  v[i] = f.value(x[i]);
	else if (fit == fit::relative)
	  v[i] = 1;
	else if (fit == fit::scaling)
	  v[i] = f.value(x[i])*scaling_function.value(x[i]);
	else if (fit == fit::relative_and_cutoff) {
	  v[i] = 1;
	  if (x[i]>x_cutoff_) {
	    v[i] = 0;
	  }
	}  
	for (size_t j=0; j<n_terms; ++j) {
	  if (fit == fit::absolute)
	    mat[i][j] = d[j];
	  else if (fit == fit::relative)
	    mat[i][j] = d[j]/f.value(x[i]);
	  else if (fit == fit::scaling)
	    mat[i][j] = d[j]*scaling_function.value(x[i]);
	  else if (fit == fit::relative_and_cutoff) {
	    mat[i][j] = d[j]/f.value(x[i]);
	    if (x[i]>x_cutoff_) {
	      mat[i][j] = 0;
	    }
	  }
	}
      }
      size_t n_zeros = wigner_d::leading_zeros(m_,n_);
      if (n_zeros > 0) {
	for (size_t i=0; i<n_zeros; i++) {
	  mat = linalg::remove_column(0, mat);
	}
      }
      coefficients_ = linalg::solve(mat,v);
      if (n_zeros > 0) {
	std::vector<double> zeros(n_zeros,0.0);
	coefficients_.insert(coefficients_.begin(),zeros.begin(),zeros.end());
      }
    }
    stdvector coefficients() {
      return coefficients_;
    }
    double value(double x) {
      stdvector d = wigner_d(x, m_, n_, coefficients_.size()).terms();
      return vec::sum(d*coefficients_);
    }
  };

  stdvector wigner_evaluate(const stdvector& coefficients, const stdvector& x, int m, int n) {
    stdvector v(x.size());
    for (size_t i=0; i<x.size(); ++i)
      v[i] = vec::sum(wigner_d(x[i],m,n,coefficients.size()).terms() * coefficients);
    return v;
  }
}

#endif
