#ifndef flick_wall_interactor
#define flick_wall_interactor

#include "../geometry/volume.hpp"
#include "../component/content.hpp"
#include "../component/emitter.hpp"

namespace flick {
  class wall_interactor {
    geometry::navigator<content>& nav_;
    std::optional<pose> next_wall_intersection_;
    geometry::volume<content>* next_volume_;
    const coating::base* coating_;
    radiation_package& rp_;
    uniform_random& rnd_;
    coating::lambert_angle_generator& lag_;
    bool is_reflected_{false};
    bool is_transmitted_{true};
  public:
    wall_interactor(geometry::navigator<content>& nav,
		    radiation_package& rp,
		    uniform_random& ur,
		    coating::lambert_angle_generator& lag)
      : nav_{nav}, rp_{rp}, rnd_{ur}, lag_{lag} {   
      next_wall_intersection_ = nav_.next_intersection(rp_.pose());
      next_volume_ = &nav_.next_volume(rp_.pose());
      geometry::volume<content>* v = next_volume_;
      if (is_moving_outward())
      	v = &nav_.current_volume();
      if (v->content().has_coating() && next_wall_intersection_.has_value()) {
	coating_ = &v->content().coating();
	lag_ = coating::lambert_angle_generator{rnd_(0,1),rnd_(0,1)};
	double r = rnd_(0,1);
	is_reflected_ = (r < coating_->reflectivity());
	is_transmitted_ = (1-r < coating_->transmissivity());
      }
    }
    void move_to_wall() {
      if (next_wall_intersection_.has_value())
	rp_.move_to((*next_wall_intersection_).position());
    }
    void move_close_to_wall() {
      move_to_wall();
      rp_.move(-nav_.current_volume().small_step()); //change to moving allong surf norm
    }
    void move_through_wall() {
      move_to_wall();
      rp_.move(+nav_.current_volume().small_step()); //change to moving allong surf norm
    }    
    void interact_with_wall() {
      unit_vector n = facing_surface_normal();
      if (is_reflected_) {
	move_close_to_wall();
	reorient_traveling_direction(n);
	//reshape_polarization(n);
	likelihood_scale_intensity();
      } else if (is_transmitted_) {
	//std::cout << "transmitted" << std::endl;
	if (nav_.current_volume().content().has_coating() &&
	    next_wall_intersection_.has_value()) {
	  reorient_traveling_direction(-n);
	  reshape_polarization(-n);
	  likelihood_scale_intensity();
	}
	move_through_wall();
	//std::cout << "moved through wall" << std::endl;
	if (next_wall_intersection_.has_value())
	  increment_activated_receiver();
	//std::cout << "will go next volume" << std::endl;
	nav_.go_to(*next_volume_);
	//std::cout << "in next volume" << std::endl;
      } else {
	absorb_radiation_package();
      }
    }
    void increment_activated_receiver() {
      if (is_moving_inward())
	next_volume_->content().inward_receiver().receive(rp_);
      else if (is_moving_outward())
	nav_.current_volume().content().outward_receiver().receive(rp_);
    }
    auto& current_volume() const {
      return nav_.current_volume();
    }  
    unit_vector facing_surface_normal() const {
      unit_vector n = (*next_wall_intersection_).z_direction();
      if (is_moving_outward())
	n = -n;
      return n;
    }
    void absorb_radiation_package() {
      rp_.scale_intensity(0);
    }
    bool is_moving_inward() const {
      if (!next_wall_intersection_.has_value())
	return false;
      if (nav_.is_moving_inward(*next_wall_intersection_,rp_.pose()))
	return true;
      return false;
    }    
    bool is_moving_outward() const {
      return !is_moving_inward();
    }
    void set_x_axis_parallel_with_plane_of_incidence(const unit_vector&
						     surface_normal) {
      pose p = rp_.pose();
      double ang = acos(dot(surface_normal, p.y_direction()));
      rp_.rotate_about_local_z(ang);
    }
    void reorient_traveling_direction(const unit_vector& facing_surface_normal) {      
      rp_.traveling_direction(facing_surface_normal);
      rp_.rotate_about_local_z(lag_.azimuth_angle());
      rp_.rotate_about_local_x(lag_.polar_angle());
    }
    void reshape_polarization(const unit_vector& facing_surface_normal) {
      set_x_axis_parallel_with_plane_of_incidence(facing_surface_normal);
      if (is_reflected_)
	rp_.reshape_polarization(coating_->reflection_mueller_matrix());
      else if (is_transmitted_)
	rp_.reshape_polarization(coating_->transmission_mueller_matrix());	
    }
    void likelihood_scale_intensity() {
      double f = 1;
      if (is_reflected_) {
	double brdf = coating_->reflection_mueller_matrix().value(0,0);
	f = brdf/lag_.brdf();
      } else if (is_transmitted_) {
	double brdf = coating_->transmission_mueller_matrix().value(0,0);
	f = brdf/lag_.brdf();
      }
      rp_.scale_intensity(f);
    }
    double distance_to_wall() {
      if (next_wall_intersection_.has_value())
	return norm((*next_wall_intersection_).position()-rp_.pose().position());
      return std::numeric_limits<double>::max(); 
    }
  };
}

#endif