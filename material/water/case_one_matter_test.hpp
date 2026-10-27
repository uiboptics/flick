#include "case_one_matter.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(case_one_matter_test) {
    using namespace units;
    material::case_one_matter m;
    m.set_wavelength(440e-9);

    material::phytoplankton p;
    p.set_wavelength(440e-9);
    check_close(m.absorption_coefficient(),
                1.2*p.absorption_coefficient(), 1e-10);
    check_close(m.scattering_coefficient(),
                p.scattering_coefficient(), 1e-10);

  } end_test_case()
}
