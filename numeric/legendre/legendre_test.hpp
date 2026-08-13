#include "legendre.hpp"
#include "../physics_function.hpp"
#include "../constants.hpp"

namespace flick {
  struct test_f {
    double g = 0.3;
    double value(double mu) const {
      return henyey_greenstein(g).value(mu);
    }
  };
  begin_test_case(legendre_test_A) {
    auto p = legendre(129,{0.5});
    check_close(p.value(64,0), -0.0755712307, 1e-8_pct);
    check_close(p.value(128,0), -0.01953466424, 1e-8_pct);
  } end_test_case()
  
  begin_test_case(legendre_test_B) {
    size_t n_terms = 12;
    size_t log2_n_points = 7;
    test_f f;
    check_close(gl_integral(f,log2_n_points).value(-1,1),
    		1/(2*constants::pi),0.2_pct);
    std::vector<double> terms = legendre_expansion(std::make_shared<test_f>(),n_terms,log2_n_points); 
    for (size_t i=0; i < terms.size(); ++i) {
      check_close(terms[i], pow(test_f().g,i)*(2.*i+1)/2/(2*constants::pi), 0.01_pct);
    }
    std::vector<double> x = read_quadrature(log2_n_points).column(0);
    std::vector<double> v = legendre_evaluation(terms).values(x);
    for (size_t i=0; i < v.size(); ++i) {
      check_close(v[i], test_f().value(x[i]), 0.02_pct);
    }
  } end_test_case()
}
