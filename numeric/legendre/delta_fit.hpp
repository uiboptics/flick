#ifndef flick_delta_fit
#define flick_delta_fit

#include "../constants.hpp"
#include "legendre.hpp"
#include <armadillo>
 
namespace flick {
  template<class Function>
  class delta_fit {
    std::vector<double> coefficients_;
  public:
    delta_fit(const Function& f, int n_terms)
      : coefficients_(n_terms) {
      int n_points = pow(n_terms,1.6);
      std::vector<double> x = range(-1,1,n_points).linspace();      
      arma::mat m(x.size(), n_terms);
      legendre p(n_terms, x);
      for (size_t i=0; i<x.size(); ++i) {
	for (size_t j=0; j<n_terms; ++j) {
	  m(i,j) = p.value(j,i)/f.value(x[i]); 
	}
      }
      arma::vec v = arma::ones(x.size());      
      arma::vec c = arma::solve(m,v);
      coefficients_ = arma::conv_to<std::vector<double>>::from(c);
    }
    
    std::vector<double> coefficients() {
      return coefficients_;
    }
    std::vector<double> function_values(const std::vector<double>& x) {
      //return flick::abs(legendre_evaluation(coefficients_).values(x));
      return legendre_evaluation(coefficients_).values(x);
    }
    double scaling_factor() const {
      return coefficients_.at(0)*4*constants::pi;
    }
  };  
}

#endif