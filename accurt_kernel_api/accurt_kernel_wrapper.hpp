#ifndef ACCURT_KERNEL
#define ACCURT_KERNEL

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace flick {
  extern "C" {
    void kernel_(int&, int&, double&, double*, double*, double*, 
		    double&, double&, int&, int&, double*, int&, 
		    int&, int&, int&, double*, int&, double*, 
		    int&, double*, double&, double&, double&, double&, 
		    int&, double&, double&, double&, double&, 
		    int&, int&, int&, int*, 
		    int&, int&, int&, int&, int&, 
		    double*, double*, double*, double*, double*, double*,  
		    double*, double*, double*, double*, int&,
		    double*, double*, double*, double*, 
		    double*, int&,
		    int&, int&, double&, double&, int&, double*);
  }
  
  class accurt_kernel_wrapper
  /* Wrapper for the AccuRT FORTRAN kernel that solves the scalar
     radiative transfer equation, using the discrete ordinate
     method. See Thomas and Stamnes, Radiative transfer in the
     atmosphere and ocean, Cambridge, 2002.
  */
  {
  private:
    int n_upper_layers = 0;
    int n_total_layers = 0;
    double refractive_index_ = 1.0;
    std::vector<double> optical_depths;
    std::vector<double> single_scattering_albedos;
    std::vector<std::vector<double>> phase_function_moments;

    double wvmlo = 0;
    double wvmhi = 0;
    int do_user_optical_depths = 1;
    int n_user_optical_depths = 0;
    std::vector<double> user_optical_depths;
    int n_streams_upper = 4;
    
    int n_streams_lower = 16;
    int do_user_angles = 0;
    int n_cos_polar_angles_upper = 0;
    std::vector<double> cos_polar_angles_upper;
    int n_cos_polar_angles_lower = 0;
    std::vector<double> cos_polar_angles_lower;

    int n_azimuth_angles = 0;
    std::vector<double> azimuth_angles;
    double incident_beam_irradiance = 1;
    double cos_source_zenith_angle = 1;
    double solar_azimuth_angle = 0;
    double fisot = 0;

    int lamber = 1;
    double bottom_albedo_ = 0;
    double btemp = 0;
    double ttemp = 0;
    double temis = 0;

    int do_delta_m = 0;
    int no_planck = 1;
    int do_flux_only = 1;
    std::vector<int> do_printing = std::vector<int>(7, 0);
    
    int max_n_total_layers = n_total_layers;
    int max_n_user_optical_depths = n_user_optical_depths;
    int max_n_user_cos_polar_angles = n_cos_polar_angles_lower;
    int max_n_cos_polar_angles = n_cos_polar_angles_lower;
    int max_n_azimuth_angles = n_azimuth_angles;

    std::vector<double> beam_plane_irradiance_down_output;
    std::vector<double> diffuse_plane_irradiance_down;
    std::vector<double> beam_plane_irradiance_up_output;
    std::vector<double> diffuse_plane_irradiance_up;
    std::vector<double> dfdt;
    std::vector<double> scalar_irradiance;
    
    std::vector<double> user_radiance;
    std::vector<double> user_azimuthally_averaged_radiance;
    std::vector<double> scalar_irradiance_up_output;
    std::vector<double> scalar_irradiance_down_output;
    int onlyim = 0;

    std::vector<double> beam_sine_irradiance_down;
    std::vector<double> diffuse_sine_irradiance_down;
    std::vector<double> beam_sine_irradiance_up;
    std::vector<double> diffuse_sine_irradiance_up;

    std::vector<double> geometrical_depth; // Input
    int do_pseudo_spherical_correction = 0;
    int do_2d_rough_surface = 0;
    int do_1d_rough_surface = 0;
    double wind_speed = 0;
    double wind_direction = 0;

    int do_lpick = 0;
    std::vector<double> radiance_fourier_components;

    static int checked_size(const std::vector<double>& values) {
      if (values.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("Array is too large for the Fortran API");
      }
      return static_cast<int>(values.size());
    }

    static std::vector<double> flatten(
        const std::vector<std::vector<double>>& matrix) {
      std::vector<double> result;
      for (const auto& row : matrix)
        result.insert(result.end(), row.begin(), row.end());
      return result;
    }

    void update_layer_counts() {
      n_total_layers = checked_size(optical_depths);
      // New layers initially belong to the upper slab. Call
      // layers_upper_slab() after insertion to place trailing layers in the
      // lower slab.
      n_upper_layers = n_total_layers;
    }

    std::size_t required_phase_moment_count() const noexcept {
      return static_cast<std::size_t>(
          std::max(n_streams_upper, n_streams_lower) + 1);
    }

    void resize_phase_function_moments() {
      const auto required_moments = required_phase_moment_count();
      for (auto& moments : phase_function_moments)
        moments.resize(required_moments, 0.0);
    }

    void validate_inputs() const {
      if (single_scattering_albedos.size() != optical_depths.size() ||
          phase_function_moments.size() != optical_depths.size())
        throw std::logic_error("Layer input arrays have inconsistent sizes");

      const auto required_moments = required_phase_moment_count();
      for (const auto& moments : phase_function_moments) {
        if (moments.size() != required_moments)
          throw std::logic_error(
              "Each layer must have n_streams + 1 phase moments");
      }
    }

    void prepare_buffers(int kernel_layer_count) {
      n_user_optical_depths = checked_size(user_optical_depths);
      if (do_user_angles) {
        n_cos_polar_angles_upper = checked_size(cos_polar_angles_upper);
        n_cos_polar_angles_lower = checked_size(cos_polar_angles_lower);
      } else {
        // With user angles disabled, the kernel returns its quadrature
        // angles through these arrays and still expects nonzero dimensions.
        n_cos_polar_angles_upper = n_streams_upper;
        n_cos_polar_angles_lower = n_streams_lower;
        cos_polar_angles_upper.resize(
            static_cast<std::size_t>(n_cos_polar_angles_upper));
        cos_polar_angles_lower.resize(
            static_cast<std::size_t>(n_cos_polar_angles_lower));
      }
      n_azimuth_angles = checked_size(azimuth_angles);

      max_n_total_layers = kernel_layer_count;
      max_n_user_optical_depths = n_user_optical_depths;
      max_n_user_cos_polar_angles = std::max(
          {n_cos_polar_angles_upper, n_cos_polar_angles_lower,
           n_streams_upper, n_streams_lower});
      max_n_cos_polar_angles = std::max(n_streams_upper, n_streams_lower);
      max_n_azimuth_angles = n_azimuth_angles;

      const auto ntau = static_cast<std::size_t>(n_user_optical_depths);
      const auto n_angles = static_cast<std::size_t>(max_n_cos_polar_angles);
      const auto n_phi = static_cast<std::size_t>(n_azimuth_angles);

      beam_plane_irradiance_down_output.resize(ntau);
      diffuse_plane_irradiance_down.resize(ntau);
      beam_plane_irradiance_up_output.resize(ntau);
      diffuse_plane_irradiance_up.resize(ntau);
      dfdt.resize(ntau);
      scalar_irradiance.resize(ntau);
      scalar_irradiance_up_output.resize(ntau);
      scalar_irradiance_down_output.resize(ntau);

      beam_sine_irradiance_down.resize(ntau);
      diffuse_sine_irradiance_down.resize(ntau);
      beam_sine_irradiance_up.resize(ntau);
      diffuse_sine_irradiance_up.resize(ntau);

      geometrical_depth.resize(static_cast<std::size_t>(kernel_layer_count));
      user_radiance.resize(ntau * n_angles * n_phi);
      user_azimuthally_averaged_radiance.resize(ntau * n_angles);
      radiance_fourier_components.resize(
          static_cast<std::size_t>(max_n_user_cos_polar_angles) * ntau *
          (static_cast<std::size_t>(max_n_cos_polar_angles) + 1));
    }

  public:
    void source_zenith_angle(double angle) {
      if (!std::isfinite(angle) || angle < 0.0 ||
          angle > std::acos(-1.0)) {
        throw std::invalid_argument(
            "Source zenith angle must be finite and in [0, pi]");
      }
      // Accept either a conventional source zenith angle (zero is vertically
      // downward) or its propagation-direction supplement. The kernel uses a
      // positive downward direction cosine.
      cos_source_zenith_angle = std::abs(std::cos(angle));
      if (cos_source_zenith_angle <= 0.0)
        throw std::invalid_argument(
            "Source direction must point into the upper boundary");
    }

    void insert_user_optical_depth(double optical_depth) {
      if (optical_depth < 0.0)
        throw std::invalid_argument("User optical depth cannot be negative");
      user_optical_depths.push_back(optical_depth);
    }

    void refractive_index(double value) {
      if (value < 1.0)
        throw std::invalid_argument("Refractive index must be at least 1");
      refractive_index_ = value;
    }

    void bottom_albedo(double value) {
      if (value < 0.0 || value > 1.0)
        throw std::invalid_argument("Bottom albedo must be in [0, 1]");
      bottom_albedo_ = value;
    }

    void streams_upper_slab(int number_of_streams) {
      if (number_of_streams < 2 || number_of_streams % 2 != 0)
        throw std::invalid_argument(
            "The number of streams must be an even integer of at least 2");
      if (number_of_streams > 20)
        throw std::invalid_argument(
            "Upper-slab streams cannot exceed 20 because the lower-slab "
            "stream count is its square and the kernel limit is 400");

      n_streams_upper = number_of_streams;
      n_streams_lower = number_of_streams * number_of_streams;
      resize_phase_function_moments();
    }

    int streams_upper_slab() const noexcept {
      return n_streams_upper;
    }

    void layers_upper_slab(int number_of_layers) {
      if (number_of_layers < 0 || number_of_layers > n_total_layers)
        throw std::invalid_argument(
            "The number of upper-slab layers must be between 0 and the "
            "total number of layers");

      n_upper_layers = number_of_layers;
    }

    std::size_t phase_moment_count(std::size_t layer) const {
      return phase_function_moments.at(layer).size();
    }

    void insert_layer_at_bottom(
        double optical_thickness,
        double single_scattering_albedo,
        std::vector<double> phase_moments) {
      if (!std::isfinite(optical_thickness) || optical_thickness < 0.0)
        throw std::invalid_argument(
            "Optical thickness must be finite and non-negative");
      if (!std::isfinite(single_scattering_albedo) ||
          single_scattering_albedo < 0.0 || single_scattering_albedo > 1.0)
        throw std::invalid_argument(
            "Single-scattering albedo must be finite and in [0, 1]");

      const auto required_moments = required_phase_moment_count();
      if (phase_moments.empty())
        throw std::invalid_argument(
            "At least one phase-function moment is required");
      phase_moments.resize(required_moments, 0.0);

      optical_depths.push_back(optical_thickness);
      single_scattering_albedos.push_back(single_scattering_albedo);
      phase_function_moments.push_back(std::move(phase_moments));
      update_layer_counts();
    }

    int number_of_layers() const noexcept {
      return n_total_layers;
    }

    void solve_rte() {
      validate_inputs();
      // The AccuRT kernel is formulated for coupled upper and lower media.
      // Its boundary equations require at least one lower-medium layer even
      // when the public model contains only upper-medium layers. Add a
      // transparent sentinel below the real layers; it does not change their
      // optical thickness and keeps the implementation detail out of the API.
      auto kernel_optical_depths = optical_depths;
      auto kernel_single_scattering_albedos = single_scattering_albedos;
      auto kernel_phase_function_moments = phase_function_moments;
      std::vector<double> transparent_phase_moments(
          required_phase_moment_count(), 0.0);
      transparent_phase_moments.front() = 1.0;

      int kernel_upper_layers = n_upper_layers;
      if (kernel_upper_layers == 0) {
        // The kernel's L1 == 0 path is singular and does not inject the
        // direct beam correctly. Put a transparent upper sentinel before
        // any real lower-slab layers. Its zero optical thickness preserves
        // all user optical-depth coordinates.
        kernel_optical_depths.insert(kernel_optical_depths.begin(), 0.0);
        kernel_single_scattering_albedos.insert(
            kernel_single_scattering_albedos.begin(), 0.0);
        kernel_phase_function_moments.insert(
            kernel_phase_function_moments.begin(),
            transparent_phase_moments);
        kernel_upper_layers = 1;
      }

      kernel_optical_depths.push_back(0.0);
      kernel_single_scattering_albedos.push_back(0.0);
      kernel_phase_function_moments.push_back(
          std::move(transparent_phase_moments));

      int kernel_total_layers = checked_size(kernel_optical_depths);
      prepare_buffers(kernel_total_layers);
      auto flattened_phase_function_moments =
          flatten(kernel_phase_function_moments);

      kernel_(kernel_upper_layers,
			    kernel_total_layers,
				    refractive_index_,
			    kernel_optical_depths.data(),
			    kernel_single_scattering_albedos.data(),
			    flattened_phase_function_moments.data(),
		    
			    wvmlo,
			    wvmhi,
			    do_user_optical_depths,
			    n_user_optical_depths,
			    user_optical_depths.data(),
			    n_streams_upper,
		    
			    n_streams_lower,
			    do_user_angles,
			    n_cos_polar_angles_upper,
			    cos_polar_angles_upper.data(),
			    n_cos_polar_angles_lower,
			    cos_polar_angles_lower.data(),
		    
			    n_azimuth_angles,
			    azimuth_angles.data(),
			    incident_beam_irradiance,
			    cos_source_zenith_angle,
			    solar_azimuth_angle,
			    fisot,                 
			    
			    lamber,
				    bottom_albedo_,
			    btemp,
			    ttemp,
			    temis,                     
			    do_delta_m,
			    
			    no_planck,
			    do_flux_only,
			    do_printing.data(),                         
			    max_n_total_layers,
			    max_n_user_optical_depths,
			    max_n_user_cos_polar_angles,
		    
			    max_n_cos_polar_angles,
			    max_n_azimuth_angles,
			    beam_plane_irradiance_down_output.data(),
			    diffuse_plane_irradiance_down.data(),
			    beam_plane_irradiance_up_output.data(),
			    diffuse_plane_irradiance_up.data(),
			    dfdt.data(),
			    scalar_irradiance.data(),
		    
			    user_radiance.data(),
			    user_azimuthally_averaged_radiance.data(),
			    scalar_irradiance_up_output.data(),
			    scalar_irradiance_down_output.data(),
			    onlyim,
		    
			    beam_sine_irradiance_down.data(),
			    diffuse_sine_irradiance_down.data(),
			    beam_sine_irradiance_up.data(),
			    diffuse_sine_irradiance_up.data(),
		    
			    geometrical_depth.data(),
			    do_pseudo_spherical_correction,                      
		    
			    do_2d_rough_surface,
			    do_1d_rough_surface,
			    wind_speed,
			    wind_direction,
		    
			    do_lpick, radiance_fourier_components.data()); 
    }
    std::vector<double> plane_irradiance_up() const {
      std::vector<double> result(beam_plane_irradiance_up_output.size());
      std::transform(beam_plane_irradiance_up_output.begin(),
                     beam_plane_irradiance_up_output.end(),
                     diffuse_plane_irradiance_up.begin(), result.begin(),
                     std::plus<double>());
      return result;
    }

    std::vector<double> plane_irradiance_down() const {
      std::vector<double> result(beam_plane_irradiance_down_output.size());
      std::transform(beam_plane_irradiance_down_output.begin(),
                     beam_plane_irradiance_down_output.end(),
                     diffuse_plane_irradiance_down.begin(), result.begin(),
                     std::plus<double>());
      return result;
    }

    const std::vector<double>& scalar_irradiance_up() const noexcept {
      return scalar_irradiance_up_output;
    }

    const std::vector<double>& beam_plane_irradiance_up() const noexcept {
      return beam_plane_irradiance_up_output;
    }

    const std::vector<double>& beam_plane_irradiance_down() const noexcept {
      return beam_plane_irradiance_down_output;
    }

    const std::vector<double>& scalar_irradiance_down() const noexcept {
      return scalar_irradiance_down_output;
    }
  };
}
#endif
