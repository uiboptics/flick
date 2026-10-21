#ifndef ACCURT_KERNEL
#define ACCURT_KERNEL

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace flick {
  struct dirac_delta_beam {
    double polar_angle;
    double azimuth_angle;
    double irradiance;
  };

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
    
    int n_streams_lower = 4;
    int do_user_angles = 0;
    int n_cos_polar_angles_upper = 0;
    std::vector<double> cos_polar_angles_upper_output;
    int n_cos_polar_angles_lower = 0;
    std::vector<double> cos_polar_angles_lower_output;

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
    std::vector<double> radiance_fourier_components_output;

    static int checked_size(const std::vector<double>& values) {
      if (values.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("Array is too large for the Fortran API");
      }
      return static_cast<int>(values.size());
    }

    static int lower_stream_count(int upper_stream_count,
                                  double refractive_index) {
      constexpr int kernel_stream_limit = 400;
      const double required = static_cast<double>(upper_stream_count) *
                              refractive_index * refractive_index;
      if (!std::isfinite(required) || required > kernel_stream_limit)
        throw std::invalid_argument(
            "The refractive-index-scaled lower-slab stream count cannot "
            "exceed the kernel limit of 400");

      // Round upward to avoid undersampling and then to an even count, as
      // required by the double-Gauss quadrature used by the kernel.
      return 2 * static_cast<int>(std::ceil(required / 2.0));
    }

    static std::vector<double> flatten(
        const std::vector<std::vector<double>>& matrix) {
      std::vector<double> result;
      for (const auto& row : matrix)
        result.insert(result.end(), row.begin(), row.end());
      return result;
    }

    double upper_slab_optical_depth() const noexcept {
      double result = 0.0;
      for (int layer = 0; layer < n_upper_layers; ++layer)
        result += optical_depths[static_cast<std::size_t>(layer)];
      return result;
    }

    bool output_is_in_upper_slab(double optical_depth) const noexcept {
      const double boundary = upper_slab_optical_depth();
      const double scale = std::max({1.0, std::abs(boundary),
                                     std::abs(optical_depth)});
      const double tolerance =
          16.0 * std::numeric_limits<double>::epsilon() * scale;
      return optical_depth <= boundary + tolerance;
    }

    double downward_polar_angle(bool upper_slab) const noexcept {
      if (upper_slab)
        return std::acos(-cos_source_zenith_angle);

      const double sin_upper = std::sqrt(std::max(
          0.0, 1.0 - cos_source_zenith_angle * cos_source_zenith_angle));
      const double sin_lower = sin_upper / refractive_index_;
      const double cos_lower = std::sqrt(std::max(
          0.0, 1.0 - sin_lower * sin_lower));
      return std::acos(-cos_lower);
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
        n_cos_polar_angles_upper =
            checked_size(cos_polar_angles_upper_output);
        n_cos_polar_angles_lower =
            checked_size(cos_polar_angles_lower_output);
      } else {
        // With user angles disabled, the kernel returns its quadrature
        // angles through these arrays and still expects nonzero dimensions.
        n_cos_polar_angles_upper = n_streams_upper;
        n_cos_polar_angles_lower = n_streams_lower;
        cos_polar_angles_upper_output.resize(
            static_cast<std::size_t>(n_cos_polar_angles_upper));
        cos_polar_angles_lower_output.resize(
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
      radiance_fourier_components_output.assign(
          static_cast<std::size_t>(max_n_user_cos_polar_angles) * ntau *
              (static_cast<std::size_t>(max_n_cos_polar_angles) + 1),
          0.0);
    }

  public:
    void source_zenith_angle(double angle) {
      const double pi = std::acos(-1.0);
      if (!std::isfinite(angle) || angle <= pi / 2.0 || angle > pi) {
        throw std::invalid_argument(
            "Source zenith angle must be finite and in (pi/2, pi]");
      }

      // The public angle describes the light's propagation direction: pi is
      // vertically downward. The kernel expects its positive downward cosine.
      cos_source_zenith_angle = -std::cos(angle);
    }

    void append_output_optical_depth(double optical_depth) {
      if (optical_depth < 0.0)
        throw std::invalid_argument("User optical depth cannot be negative");
      user_optical_depths.push_back(optical_depth);
    }

    void refractive_index(double value) {
      if (!std::isfinite(value) || value < 1.0)
        throw std::invalid_argument(
            "Refractive index must be finite and at least 1");
      const int lower_streams = lower_stream_count(n_streams_upper, value);
      refractive_index_ = value;
      n_streams_lower = lower_streams;
      resize_phase_function_moments();
    }

    void bottom_albedo(double value) {
      if (value < 0.0 || value > 1.0)
        throw std::invalid_argument("Bottom albedo must be in [0, 1]");
      bottom_albedo_ = value;
    }
    
    void enable_azimuthal_radiance_dependence(bool enabeled) {
      do_flux_only = enabeled ? 0 : 1;
    }

    void streams_upper_slab(int number_of_streams) {
      if (number_of_streams < 2 || number_of_streams % 2 != 0)
        throw std::invalid_argument(
            "The number of streams must be an even integer of at least 2");
      if (number_of_streams > 400)
        throw std::invalid_argument(
            "Upper-slab streams cannot exceed the kernel limit of 400");

      const int lower_streams =
          lower_stream_count(number_of_streams, refractive_index_);
      n_streams_upper = number_of_streams;
      n_streams_lower = lower_streams;
      resize_phase_function_moments();
    }

    int streams_upper_slab() const noexcept {
      return n_streams_upper;
    }

    int streams_lower_slab() const noexcept {
      return n_streams_lower;
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
			    cos_polar_angles_upper_output.data(),
			    n_cos_polar_angles_lower,
			    cos_polar_angles_lower_output.data(),
		    
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
		    
			    do_lpick, radiance_fourier_components_output.data()); 
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

    std::vector<double> beam_irradiance_down() const {
      std::vector<double> result(beam_plane_irradiance_down_output.size());
      std::transform(beam_plane_irradiance_down_output.begin(),
                     beam_plane_irradiance_down_output.end(),
                     beam_sine_irradiance_down.begin(), result.begin(),
                     [](double normal, double tangential) {
                       return std::hypot(normal, tangential);
                     });
      return result;
    }

    std::vector<double> beam_irradiance_up() const {
      std::vector<double> result(beam_plane_irradiance_up_output.size());
      std::transform(beam_plane_irradiance_up_output.begin(),
                     beam_plane_irradiance_up_output.end(),
                     beam_sine_irradiance_up.begin(), result.begin(),
                     [](double normal, double tangential) {
                       return std::hypot(normal, tangential);
                     });
      return result;
    }

    std::vector<dirac_delta_beam> dirac_delta_beam_down() const {
      const auto irradiances = beam_irradiance_down();
      std::vector<dirac_delta_beam> result;
      result.reserve(irradiances.size());
      for (std::size_t i = 0; i < irradiances.size(); ++i) {
        result.push_back({
            downward_polar_angle(output_is_in_upper_slab(
                user_optical_depths.at(i))),
            solar_azimuth_angle,
            irradiances[i]});
      }
      return result;
    }

    std::vector<dirac_delta_beam> dirac_delta_beam_up() const {
      const double pi = std::acos(-1.0);
      const auto irradiances = beam_irradiance_up();
      std::vector<dirac_delta_beam> result;
      result.reserve(irradiances.size());
      for (std::size_t i = 0; i < irradiances.size(); ++i) {
        const double downward =
            downward_polar_angle(output_is_in_upper_slab(
                user_optical_depths.at(i)));
        result.push_back(
            {pi - downward, solar_azimuth_angle, irradiances[i]});
      }
      return result;
    }

    std::vector<std::vector<dirac_delta_beam>> dirac_delta_beams() const {
      const auto downward = dirac_delta_beam_down();
      const auto upward = dirac_delta_beam_up();
      std::vector<std::vector<dirac_delta_beam>> result;
      result.reserve(user_optical_depths.size());
      for (std::size_t i = 0; i < user_optical_depths.size(); ++i) {
        result.push_back({downward.at(i)});
        if (output_is_in_upper_slab(user_optical_depths[i]))
          result.back().push_back(upward.at(i));
      }
      return result;
    }

    const std::vector<double>& scalar_irradiance_down() const noexcept {
      return scalar_irradiance_down_output;
    }

    const std::vector<double>&
    quadrature_cosines_upper_slab() const noexcept {
      return cos_polar_angles_upper_output;
    }

    const std::vector<double>&
    quadrature_cosines_lower_slab() const noexcept {
      return cos_polar_angles_lower_output;
    }

    std::vector<std::vector<double>> quadrature_cosines() const {
      std::vector<std::vector<double>> result;
      result.reserve(user_optical_depths.size());
      for (double optical_depth : user_optical_depths) {
        result.push_back(output_is_in_upper_slab(optical_depth)
                             ? cos_polar_angles_upper_output
                             : cos_polar_angles_lower_output);
      }
      return result;
    }

    std::vector<std::vector<std::vector<double>>>
    radiance_fourier_components() const {
      const auto n_depths =
          static_cast<std::size_t>(max_n_user_optical_depths);
      const auto n_angles =
          static_cast<std::size_t>(max_n_user_cos_polar_angles);
      const auto n_orders =
          static_cast<std::size_t>(max_n_cos_polar_angles) + 1;

      std::vector<std::vector<std::vector<double>>> result(
          n_depths,
          std::vector<std::vector<double>>(
              n_angles, std::vector<double>(n_orders)));

      // Fortran stores UUM(angle, depth, order), with the first index
      // contiguous. Present it in the more convenient C++ order
      // [output depth][polar angle][Fourier order].
      for (std::size_t order = 0; order < n_orders; ++order) {
        for (std::size_t depth = 0; depth < n_depths; ++depth) {
          for (std::size_t angle = 0; angle < n_angles; ++angle) {
            const auto flat_index =
                angle + n_angles * (depth + n_depths * order);
            result[depth][angle][order] =
                radiance_fourier_components_output[flat_index];
          }
        }
      }
      return result;
    }
  };

  class fourier_field {
    static constexpr double pi = std::numbers::pi;
    std::vector<dirac_delta_beam> dirac_delta_beams_;
    std::vector<std::vector<double>> components_;
    std::vector<double> quadrature_cosines_;
    std::vector<double> quadrature_weights_;
    std::array<double, 9> frame_rotation_ =
        {1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0};

    bool frame_is_unrotated() const {
      constexpr std::array<double, 9> identity =
          {1.0, 0.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 0.0, 1.0};
      for (std::size_t i = 0; i < identity.size(); ++i) {
        if (std::abs(frame_rotation_[i] - identity[i]) > 1e-14)
          return false;
      }
      return true;
    }

    void rotate_frame(const std::array<double, 9>& rotation) {
      std::array<double, 9> product{};
      for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
          for (std::size_t k = 0; k < 3; ++k)
            product[3 * row + column] +=
                frame_rotation_[3 * row + k] * rotation[3 * k + column];
        }
      }
      const double two_pi = 2.0 * std::acos(-1.0);
      for (auto& beam : dirac_delta_beams_) {
        const double sine = std::sin(beam.polar_angle);
        const std::array<double, 3> old_direction = {
            sine * std::cos(beam.azimuth_angle),
            sine * std::sin(beam.azimuth_angle),
            std::cos(beam.polar_angle)};
        std::array<double, 3> new_direction{};
        // rotation maps new-frame coordinates to old-frame coordinates, so
        // its transpose maps a stored beam into the new frame.
        for (std::size_t column = 0; column < 3; ++column) {
          for (std::size_t row = 0; row < 3; ++row)
            new_direction[column] +=
                rotation[3 * row + column] * old_direction[row];
        }
        beam.polar_angle =
            std::acos(std::clamp(new_direction[2], -1.0, 1.0));
        beam.azimuth_angle =
            std::atan2(new_direction[1], new_direction[0]);
        if (beam.azimuth_angle < 0.0)
          beam.azimuth_angle += two_pi;
      }
      frame_rotation_ = product;
    }

    static bool azimuth_is_in_interval(double azimuth,
                                       double phi1, double phi2) {
      const double two_pi = 2.0 * std::acos(-1.0);
      const double tolerance =
          100.0 * std::numeric_limits<double>::epsilon();
      if (phi2 - phi1 >= two_pi - tolerance)
        return true;
      const double turns = std::ceil((phi1 - azimuth) / two_pi);
      const double equivalent = azimuth + turns * two_pi;
      return equivalent <= phi2 + tolerance;
    }

    double dirac_delta_integral(double theta1, double theta2,
                                double phi1, double phi2,
                                bool cosine_weighted) const {
      const double tolerance =
          100.0 * std::numeric_limits<double>::epsilon();
      double result = 0.0;
      for (const auto& beam : dirac_delta_beams_) {
        if (beam.polar_angle < theta1 - tolerance ||
            beam.polar_angle > theta2 + tolerance ||
            !azimuth_is_in_interval(beam.azimuth_angle, phi1, phi2))
          continue;
        result += beam.irradiance *
            (cosine_weighted
                 ? std::max(0.0, std::cos(beam.polar_angle))
                 : 1.0);
      }
      return result;
    }

    double beam_hemisphere_integral(bool upper_hemisphere,
                                    bool cosine_weighted) const noexcept {
      double result = 0.0;
      for (const auto& beam : dirac_delta_beams_) {
        const double mu = std::cos(beam.polar_angle);
        // A horizontal beam belongs to the upper hemisphere, matching the
        // public methods' documented <= pi/2 and > pi/2 split.
        if ((upper_hemisphere && mu < 0.0) ||
            (!upper_hemisphere && mu >= 0.0))
          continue;
        result += beam.irradiance *
                  (cosine_weighted ? std::abs(mu) : 1.0);
      }
      return result;
    }

    static double gauss_legendre_weight(double mu,
                                        std::size_t direction_count) {
      const std::size_t order = direction_count / 2;
      const double x = 2.0 * std::abs(mu) - 1.0;
      if (order == 1)
        return 1.0;
      double p_previous = 1.0;
      double p_current = x;
      for (std::size_t degree = 2; degree <= order; ++degree) {
        const double p_next =
            ((2.0 * degree - 1.0) * x * p_current -
             (degree - 1.0) * p_previous) /
            static_cast<double>(degree);
        p_previous = p_current;
        p_current = p_next;
      }
      const double derivative = static_cast<double>(order) *
          (x * p_current - p_previous) / (x * x - 1.0);
      return 1.0 / ((1.0 - x * x) * derivative * derivative);
    }

    double radiance_at_index(std::size_t angle_index, double phi) const {
      double result = 0.0;
      const auto& fourier = components_[angle_index];
      for (std::size_t order = 0; order < fourier.size(); ++order)
        result += fourier[order] *
                  std::cos(static_cast<double>(order) * phi);
      return result;
    }

    double azimuth_integral_at_index(std::size_t angle_index,
                                     double phi1, double phi2) const {
      const auto& fourier = components_[angle_index];
      double result = fourier[0] * (phi2 - phi1);
      for (std::size_t order = 1; order < fourier.size(); ++order) {
        const double m = static_cast<double>(order);
        result += fourier[order] *
                  (std::sin(m * phi2) - std::sin(m * phi1)) / m;
      }
      return result;
    }

    double azimuth_integral(double mu, double phi1, double phi2) const {
      const auto upper_iterator = std::lower_bound(
          quadrature_cosines_.begin(), quadrature_cosines_.end(), mu);
      if (upper_iterator == quadrature_cosines_.begin())
        return azimuth_integral_at_index(0, phi1, phi2);
      if (upper_iterator == quadrature_cosines_.end())
        return azimuth_integral_at_index(
            quadrature_cosines_.size() - 1, phi1, phi2);

      const auto upper = static_cast<std::size_t>(
          upper_iterator - quadrature_cosines_.begin());
      const auto lower = upper - 1;
      const double fraction =
          (mu - quadrature_cosines_[lower]) /
          (quadrature_cosines_[upper] - quadrature_cosines_[lower]);
      return (1.0 - fraction) *
                 azimuth_integral_at_index(lower, phi1, phi2) +
             fraction * azimuth_integral_at_index(upper, phi1, phi2);
    }

    double radiance_unrotated(double mu, double phi) const {
      const auto upper_iterator = std::lower_bound(
          quadrature_cosines_.begin(), quadrature_cosines_.end(), mu);
      if (upper_iterator == quadrature_cosines_.begin())
        return radiance_at_index(0, phi);
      if (upper_iterator == quadrature_cosines_.end())
        return radiance_at_index(components_.size() - 1, phi);
      const auto upper = static_cast<std::size_t>(
          upper_iterator - quadrature_cosines_.begin());
      const auto lower = upper - 1;
      const double fraction =
          (mu - quadrature_cosines_[lower]) /
          (quadrature_cosines_[upper] - quadrature_cosines_[lower]);
      return (1.0 - fraction) * radiance_at_index(lower, phi) +
             fraction * radiance_at_index(upper, phi);
    }

    double rotated_conical_integral(double theta1, double theta2,
                                    double phi1, double phi2,
                                    bool cosine_weighted) const {
      const double mu_low = std::cos(theta2);
      const double mu_high = std::cos(theta1);
      const std::size_t half = quadrature_cosines_.size() / 2;
      const std::size_t azimuth_count = 2 * components_.front().size() + 1;
      const double delta_phi = (phi2 - phi1) / azimuth_count;
      double result = 0.0;
      for (std::size_t i = 0; i < half; ++i) {
        const double x = 2.0 * quadrature_cosines_[half + i] - 1.0;
        const double mu = 0.5 * ((mu_high - mu_low) * x +
                                 mu_high + mu_low);
        const double mu_weight =
            (mu_high - mu_low) * quadrature_weights_[half + i];
        const double theta = std::acos(std::clamp(mu, -1.0, 1.0));
        for (std::size_t j = 0; j < azimuth_count; ++j) {
          const double phi = phi1 + (j + 0.5) * delta_phi;
          double value = radiance(theta, phi);
          if (cosine_weighted)
            value *= std::max(0.0, mu);
          result += mu_weight * delta_phi * value;
        }
      }
      return result;
    }

    double conical_integral(double theta1, double theta2,
                            double phi1, double phi2,
                            bool cosine_weighted) const {
      const double pi = std::acos(-1.0);
      if (!std::isfinite(theta1) || !std::isfinite(theta2) ||
          theta1 < 0.0 || theta2 > pi || theta1 > theta2)
        throw std::invalid_argument(
            "Polar-angle bounds must satisfy 0 <= theta1 <= theta2 <= pi");
      if (!std::isfinite(phi1) || !std::isfinite(phi2) || phi1 > phi2)
        throw std::invalid_argument(
            "Azimuth bounds must be finite and satisfy phi1 <= phi2");

      if (!frame_is_unrotated())
        return rotated_conical_integral(
            theta1, theta2, phi1, phi2, cosine_weighted);

      const double tolerance =
          100.0 * std::numeric_limits<double>::epsilon();
      const bool upper_hemisphere =
          std::abs(theta1) <= tolerance &&
          std::abs(theta2 - pi / 2.0) <= tolerance;
      const bool lower_hemisphere =
          std::abs(theta1 - pi / 2.0) <= tolerance &&
          std::abs(theta2 - pi) <= tolerance;
      if (upper_hemisphere || lower_hemisphere) {
        double result = 0.0;
        for (std::size_t i = 0; i < quadrature_cosines_.size(); ++i) {
          const double mu = quadrature_cosines_[i];
          if ((upper_hemisphere && mu <= 0.0) ||
              (lower_hemisphere && mu >= 0.0))
            continue;
          double value = azimuth_integral_at_index(i, phi1, phi2);
          if (cosine_weighted)
            value *= std::max(0.0, mu);
          result += quadrature_weights_[i] * value;
        }
        return result;
      }

      const double mu_low = std::cos(theta2);
      const double mu_high = std::cos(theta1);
      std::vector<double> boundaries = {mu_low};
      if (cosine_weighted && mu_low < 0.0 && mu_high > 0.0)
        boundaries.push_back(0.0);
      for (double mu : quadrature_cosines_) {
        if (mu > mu_low && mu < mu_high)
          boundaries.push_back(mu);
      }
      boundaries.push_back(mu_high);
      std::sort(boundaries.begin(), boundaries.end());

      double result = 0.0;
      for (std::size_t i = 1; i < boundaries.size(); ++i) {
        const double a = boundaries[i - 1];
        const double b = boundaries[i];
        const double fa = azimuth_integral(a, phi1, phi2);
        const double fb = azimuth_integral(b, phi1, phi2);
        if (!cosine_weighted) {
          result += 0.5 * (fa + fb) * (b - a);
        } else {
          const double slope = (fb - fa) / (b - a);
          if (b > 0.0) {
            result += 0.5 * fa * (b * b - a * a) +
                      slope * ((b * b * b - a * a * a) / 3.0 -
                               0.5 * a * (b * b - a * a));
          }
        }
      }
      return result;
    }

  public:
    using vector_matrix = std::vector<std::vector<double>>;
    fourier_field(std::vector<std::vector<double>> components,
                  std::vector<double> quadrature_cosines,
		  std::vector<dirac_delta_beam> dirac_delta_beams = {})
        : components_(std::move(components)),
          quadrature_cosines_(std::move(quadrature_cosines)),
          dirac_delta_beams_(std::move(dirac_delta_beams))
    {
      if (components_.empty() || components_.front().empty())
        throw std::invalid_argument("Fourier component field cannot be empty");
      if (quadrature_cosines_.empty() ||
          quadrature_cosines_.size() > components_.size())
        throw std::invalid_argument(
            "Quadrature grid must match the Fourier component rows");
      components_.resize(quadrature_cosines_.size());
      const auto order_count = components_.front().size();
      if (std::any_of(components_.begin(), components_.end(),
                      [order_count](const auto& row) {
                        return row.size() != order_count;
                      }))
        throw std::invalid_argument(
            "Fourier component rows must have equal sizes");
      if (!std::is_sorted(quadrature_cosines_.begin(),
                          quadrature_cosines_.end()))
        throw std::invalid_argument(
            "Quadrature direction cosines must be sorted");
      quadrature_weights_.reserve(quadrature_cosines_.size());
      for (double mu : quadrature_cosines_)
        quadrature_weights_.push_back(
            gauss_legendre_weight(mu, quadrature_cosines_.size()));
    }

    const std::vector<std::vector<double>>& components() const noexcept {
      return components_;
    }

    double radiance(double theta, double phi) const {
      const double pi = std::acos(-1.0);
      if (!std::isfinite(theta) || theta < 0.0 || theta > pi)
        throw std::invalid_argument("Polar angle must be in [0, pi]");
      if (!std::isfinite(phi))
        throw std::invalid_argument("Azimuth angle must be finite");

      const double sine = std::sin(theta);
      const std::array<double, 3> direction =
          {sine * std::cos(phi), sine * std::sin(phi), std::cos(theta)};
      std::array<double, 3> original{};
      for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column)
          original[row] +=
              frame_rotation_[3 * row + column] * direction[column];
      }
      const double original_phi = std::atan2(original[1], original[0]);
      return radiance_unrotated(
          std::clamp(original[2], -1.0, 1.0), original_phi);
    }

    void rotate_frame_about_x_axis(double angle) {
      if (!std::isfinite(angle))
        throw std::invalid_argument("Rotation angle must be finite");
      const double cosine = std::cos(angle);
      const double sine = std::sin(angle);
      rotate_frame({1.0, 0.0, 0.0,
                    0.0, cosine, -sine,
                    0.0, sine, cosine});
    }

    void rotate_frame_about_y_axis(double angle) {
      if (!std::isfinite(angle))
        throw std::invalid_argument("Rotation angle must be finite");
      const double cosine = std::cos(angle);
      const double sine = std::sin(angle);
      rotate_frame({cosine, 0.0, sine,
                    0.0, 1.0, 0.0,
                    -sine, 0.0, cosine});
    }

    void rotate_frame_about_z_axis(double angle) {
      if (!std::isfinite(angle))
        throw std::invalid_argument("Rotation angle must be finite");
      const double cosine = std::cos(angle);
      const double sine = std::sin(angle);
      rotate_frame({cosine, -sine, 0.0,
                    sine, cosine, 0.0,
                    0.0, 0.0, 1.0});
    }

    double conical_scalar_irradiance(double theta1, double theta2,
                                     double phi1, double phi2) const {
      return conical_integral(theta1, theta2, phi1, phi2, false) +
             dirac_delta_integral(theta1, theta2, phi1, phi2, false);
    }

    double conical_plane_irradiance(double theta1, double theta2,
                                    double phi1, double phi2) const {
      return conical_integral(theta1, theta2, phi1, phi2, true) +
             dirac_delta_integral(theta1, theta2, phi1, phi2, true);
    }

    double scalar_irradiance_upper_heimisphere() const {
      return conical_scalar_irradiance(0, pi/2, 0, 2*pi);
    }
    double scalar_irradiance_lower_heimisphere() const {
      return conical_scalar_irradiance(pi/2, pi, 0, 2*pi);
    }
    double plane_irradiance_upper_heimisphere() const {
      return conical_plane_irradiance(0, pi/2, 0, 2*pi);
    }
    double plane_irradiance_lower_heimisphere() const {
      return conical_plane_irradiance(pi/2, pi, 0, 2*pi);
    }
    double beam_plane_irradiance_upper_heimisphere() const {
      return beam_hemisphere_integral(true, true);
    }

    double beam_plane_irradiance_lower_heimisphere() const {
      return beam_hemisphere_integral(false, true);
    }

    double beam_scalar_irradiance_upper_heimisphere() const {
      return beam_hemisphere_integral(true, false);
    }

    double beam_scalar_irradiance_lower_heimisphere() const {
      return beam_hemisphere_integral(false, false);
    }
  };
}
#endif
