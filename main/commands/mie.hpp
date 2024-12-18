#ifndef flick_command_mie
#define flick_command_mie

#include "basic_command.hpp"
#include "../../mie/monodispersed_mie.hpp"
#include "../../mie/polydispersed_mie.hpp"
#include "../../environment/input_output.hpp"
#include <sstream>

namespace flick {
  namespace command {
    stdcomplex stoc(std::string s) {
      if(s.back()=='i')
	s.pop_back();
      double re = 0;
      double im = 0;
      stdcomplex m;
      std::stringstream ss(s);
      ss >> re  >> im;
      m = stdcomplex{re,im};
      return m;
    }
    class mie : public basic_command {
      const double pi = constants::pi;
      std::vector<double> wls_;
    public:
      mie():basic_command("mie"){};
      void run() {
	double wl = std::stod(a(1));
	stdcomplex m_host = stoc(a(2));
	stdcomplex m_sphere = stoc(a(3));	
	double median_r = std::stod(a(4));	
	double sigma = std::stod(a(5));	
	double paccuracy = std::stod(a(6));	
	std::string output_kind = a(7);
	log_normal_distribution sd(log(median_r),sigma);
	monodispersed_mie mono_mie(m_host,m_sphere,wl);
	int n_out = 6;
	if (-log10(paccuracy) > n_out)
	  n_out = -log10(paccuracy);
	if (output_kind=="absorption_cross_section") {
	  polydispersed_mie pm(mono_mie,sd);
	  pm.percentage_accuracy(paccuracy);
	  std::cout << std::setprecision(n_out)
		    << pm.absorption_cross_section() << std::endl;
	}
	else if (output_kind=="scattering_cross_section") {
	  polydispersed_mie pm(mono_mie,sd);
	  pm.percentage_accuracy(paccuracy);
	  std::cout << std::setprecision(n_out)
		    << pm.scattering_cross_section() << std::endl;
	}
	else if (output_kind=="scattering_matrix_element") {
	  int row = std::stoi(a(8));
	  int col = std::stoi(a(9));
	  int n_angs = std::stoi(a(10));
	  auto angs = range(0,constants::pi,n_angs).linspace();
	  mono_mie.angles(angs);
	  polydispersed_mie pm(mono_mie,sd);
	  pm.percentage_accuracy(paccuracy);
	  stdvector F = pm.scattering_matrix_element(row,col); 
	  std::cout << std::setprecision(n_out)
		    << pl_function{mono_mie.angles(),F};
	} else
	  error();
      }
    };
  }
}

#endif
