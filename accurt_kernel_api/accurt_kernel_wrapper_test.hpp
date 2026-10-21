#include "../environment/unit_test.hpp"
#include "accurt_kernel_wrapper.hpp"
#include <numbers>

namespace flick {
  begin_test_case(accurt_kernel_wrapper_test_A) {
    // Compare with van de Hulst 1980, vol 1, chapter 9, table 12,
    // p258, FLUX U and AVERAGE N for mu_0 = 1.0
    accurt_kernel_wrapper a;
    a.streams_upper_slab(6);
    a.append_output_optical_depth(0);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    check(a.number_of_layers() == 0);
    a.insert_layer_at_bottom(1.0, 1.0, {1.0, 0.0, 0.0, 0.0, 0.0});
    check(a.number_of_layers() == 1);
    
    a.solve_rte();
    check_close(a.plane_irradiance_up()[0], 0.34133, 0.01_pct);
    check_close(a.scalar_irradiance_up()[0]/2, 0.37869, 0.01_pct);
  } end_test_case()

  begin_test_case(accurt_test_B) {
    // Compare with van de Hulst 1980, vol 2, chapter 13, table 35,
    // p431, FLUX
    double ssalb = 0.9;
    double od = 1;
    std::vector<double> p0_75 =
      {
	1.0,
	0.75,
	0.5625,
	0.421875,
	0.31640625,
	0.2373046875,
	0.177978515625
      };
    accurt_kernel_wrapper a;
    a.streams_upper_slab(6);
    a.append_output_optical_depth(od);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    a.insert_layer_at_bottom(od, ssalb, p0_75);
    a.solve_rte();
    check_close(a.plane_irradiance_down()[0], 0.82792, 0.02_pct);
  } end_test_case()
    
  begin_test_case(accurt_test_C) {
    // Van de Hulst, chapter 13, page 418 and 419, mu_0=0.9
    double ssalb = 1;
    double od = 1;
    std::vector<double> p0_5 =
      {
	1.0,
	0.5,
	0.25,
	0.125,
	0.0625,
	0.03125,
	0.015625
      };
    accurt_kernel_wrapper a;
    double ang = std::numbers::pi-acos(0.9);
    a.source_zenith_angle(ang);
    a.streams_upper_slab(6);
    a.append_output_optical_depth(0);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    a.insert_layer_at_bottom(od, ssalb, p0_5);
    a.solve_rte();
    double albedo = a.plane_irradiance_up()[0]/a.plane_irradiance_down()[0];
    check_close(albedo, 0.20197, 0.04_pct);
  } end_test_case()
  
  begin_test_case(accurt_test_D) {
    // Van de Hulst, chapter 13, page 418 and 419, mu_0=0.9
    accurt_kernel_wrapper a;
    a.append_output_optical_depth(0);
    a.bottom_albedo(1);
    a.refractive_index(1.0);
    a.solve_rte();
    double albedo = a.plane_irradiance_up()[0]/a.plane_irradiance_down()[0];
    check_close(albedo, 1, 0.001_pct);
    check_close(a.scalar_irradiance_up()[0], 2.0, 1e-5_pct);
    check_close(a.scalar_irradiance_down()[0], 1.0, 1e-5_pct);
    check_close(a.beam_plane_irradiance_down()[0], 1, 1e-5_pct);
    check_small(a.beam_plane_irradiance_up()[0], 1e-9);
  } end_test_case()
    
  begin_test_case(accurt_test_E) {
    // Transparent refractive slab with Lambert white bottom
    // Test against AI-generated transmittance to bottom.
    accurt_kernel_wrapper a;
    double od1 = 1e-12;
    double od2 = od1;
    double ssa = 1;
    a.insert_layer_at_bottom(od1,ssa,{1,0});
    a.insert_layer_at_bottom(od2,ssa,{1,0});
    a.layers_upper_slab(1);
    a.source_zenith_angle(std::numbers::pi);
    a.append_output_optical_depth(0);
    a.append_output_optical_depth(od1+od2);
    a.streams_upper_slab(12);
    a.bottom_albedo(1);
    a.refractive_index(1.3);
    a.solve_rte();
    double T_bottom = 1.7694;
    check_close(a.plane_irradiance_up()[0], 1, 0.001_pct);
    check_close(a.plane_irradiance_down()[1], T_bottom, 0.001_pct);
  } end_test_case()

  begin_test_case(accurt_test_F) {
    // Radiance
    accurt_kernel_wrapper a;
    double od1 = 1e-12;
    double od2 = od1;
    double ssa = 1;
    a.insert_layer_at_bottom(od1,ssa,{1,0});
    a.insert_layer_at_bottom(od2,ssa,{1,0});
    a.layers_upper_slab(1);
    a.source_zenith_angle(std::numbers::pi);
    a.append_output_optical_depth(0);
    a.streams_upper_slab(16);
    a.bottom_albedo(1);
    a.refractive_index(1);
    a.solve_rte();
    fourier_field f(a.radiance_fourier_components()[0],
                    a.quadrature_cosines_upper_slab());
    check_close(f.radiance(0,0), 1/std::numbers::pi, 1e-8_pct);
    check_close(f.radiance(1,1), 1/std::numbers::pi, 1e-8_pct);
    check_small(f.radiance(3.14,0), 1e-9);
  } end_test_case()
  
  begin_test_case(accurt_test_G) {
    // Conical scalar and plane irradiance with rotations
    accurt_kernel_wrapper a;
    double od1 = 1e-12;
    double od2 = od1;
    double ssa = 1;
    const double pi = std::numbers::pi;
    a.insert_layer_at_bottom(od1,ssa,{1,0});
    a.insert_layer_at_bottom(od2,ssa,{1,0});
    a.layers_upper_slab(1);
    a.source_zenith_angle(3*pi/4);
    a.append_output_optical_depth(0);
    a.streams_upper_slab(16);
    a.bottom_albedo(1);
    a.refractive_index(1);
    a.solve_rte();
    fourier_field f(a.radiance_fourier_components()[0],
                    a.quadrature_cosines()[0],
		    a.dirac_delta_beams()[0]);

    check_close(f.conical_scalar_irradiance(pi/2,pi,0,2*pi),
		a.scalar_irradiance_down()[0], 1e-9_pct);
    check_small(f.conical_plane_irradiance(pi/2,pi,0,2*pi), 1e-12);

    check_close(f.conical_scalar_irradiance(0,pi/2,0,2*pi),
		a.scalar_irradiance_up()[0], 1e-9_pct);
    check_close(f.conical_plane_irradiance(0,pi/2,0,2*pi),
		a.plane_irradiance_up()[0], 1e-9_pct);

    check_close(f.beam_scalar_irradiance_lower_heimisphere(),
		a.beam_irradiance_down()[0], 1e-9_pct);
    check_close(f.beam_plane_irradiance_lower_heimisphere(),
		a.beam_plane_irradiance_down()[0], 1e-9_pct);
    check_small(f.beam_scalar_irradiance_upper_heimisphere(), 1e-12);
    check_small(f.beam_plane_irradiance_upper_heimisphere(), 1e-12);

    f.rotate_frame_about_x_axis(pi);
    check_close(f.conical_plane_irradiance(0,pi/2,0,2*pi),
		a.plane_irradiance_down()[0], 1e-9_pct);
    check_close(f.beam_plane_irradiance_upper_heimisphere(),
		a.beam_plane_irradiance_down()[0], 1e-9_pct);
    f.rotate_frame_about_y_axis(pi);
    check_close(f.conical_plane_irradiance(0,pi/2,0,2*pi),
		a.plane_irradiance_up()[0], 1e-9_pct);
  } end_test_case()

  begin_test_case(accurt_test_H) {
    // Beam angles
    accurt_kernel_wrapper a;
    double od1 = 1e-12;
    double od2 = od1;
    double ssa = 1;
    const double pi = std::numbers::pi;
    a.insert_layer_at_bottom(od1,ssa,{1,0});
    a.insert_layer_at_bottom(od2,ssa,{1,0});
    a.layers_upper_slab(1);
    double source_angle = 3*pi/4;
    a.source_zenith_angle(source_angle);
    a.append_output_optical_depth(0);
    a.append_output_optical_depth(od1+od2);
    a.streams_upper_slab(4);
    a.refractive_index(1.3);
    check(a.streams_lower_slab() == 8);
    a.solve_rte();
    check_close(a.dirac_delta_beam_down()[0].irradiance,
		a.beam_irradiance_down()[0],1e-6_pct);
    check_close(a.dirac_delta_beam_down()[0].polar_angle,source_angle,1e-6_pct);
    check_small(a.dirac_delta_beam_down()[0].azimuth_angle,1e-9);
    check(a.dirac_delta_beam_down()[1].polar_angle > source_angle);
    check_close(a.dirac_delta_beam_up()[0].polar_angle,pi-source_angle,1e-5_pct);
    check_small(a.dirac_delta_beam_up()[0].azimuth_angle,1e-9);
  } end_test_case()

  begin_test_case(accurt_test_I) {
    // Single scattering by a weakly scattering, isotropic layer. For normal
    // illumination and black boundaries, the upward radiance at the top is
    //   omega/(4*pi) * (1 - exp(-tau*(1 + 1/mu))) / (1 + mu).
    // Keeping omega small makes higher-order scattering negligible.
    const double pi = std::numbers::pi;
    const double optical_depth = 1.0;
    const double single_scattering_albedo = 1e-5;

    accurt_kernel_wrapper a;
    a.streams_upper_slab(22);
    a.source_zenith_angle(pi);
    a.append_output_optical_depth(0.0);
    a.bottom_albedo(0.0);
    a.refractive_index(1.0);
    a.insert_layer_at_bottom(
        optical_depth, single_scattering_albedo, {1.0});
    a.solve_rte();

    const auto quadrature = a.quadrature_cosines()[0];
    const std::size_t upward_index = 3 * quadrature.size() / 4;
    const double mu = quadrature[upward_index];
    const double expected = single_scattering_albedo / (4.0 * pi) *
        (1.0 - std::exp(-optical_depth * (1.0 + 1.0 / mu))) /
        (1.0 + mu);
    fourier_field field(a.radiance_fourier_components()[0], quadrature,
                        a.dirac_delta_beams()[0]);

    check_close(field.radiance(std::acos(mu), 0.0), expected, 0.001_pct);
  } end_test_case()

  begin_test_case(accurt_test_J) {
    // Reusing a wrapper in flux-only mode must not retain Fourier orders
    // calculated by an earlier azimuth-dependent solve.
    accurt_kernel_wrapper a;
    a.streams_upper_slab(6);
    a.source_zenith_angle(3.0 * std::numbers::pi / 4.0);
    a.append_output_optical_depth(0.0);
    a.bottom_albedo(0.0);
    a.insert_layer_at_bottom(
        1.0, 0.9, {1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125});

    a.invoke_azimuthal_radiance_dependence(true);
    a.solve_rte();
    double largest_higher_order = 0.0;
    auto components = a.radiance_fourier_components();
    for (const auto& direction : components[0]) {
      for (std::size_t order = 1; order < direction.size(); ++order)
        largest_higher_order =
            std::max(largest_higher_order, std::abs(direction[order]));
    }
    check(largest_higher_order > 1e-8);

    a.invoke_azimuthal_radiance_dependence(false);
    a.solve_rte();
    largest_higher_order = 0.0;
    components = a.radiance_fourier_components();
    for (const auto& direction : components[0]) {
      for (std::size_t order = 1; order < direction.size(); ++order)
        largest_higher_order =
            std::max(largest_higher_order, std::abs(direction[order]));
    }
    check_small(largest_higher_order, 1e-15);
  } end_test_case()

}
