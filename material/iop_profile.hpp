#ifndef flick_iop_profile
#define flick_iop_profile

#include "../numeric/constants.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/pose.hpp"

namespace flick {

  class constant_iop {
    double value_;
    const double epsilon_ = std::numeric_limits<double>::epsilon()*10;
  public:
    constant_iop(double value) : value_{value}{}
    double optical_depth(double distance) const {
      return value_ * distance;
    }
    double distance(double optical_depth) const {
      if (value_ >  epsilon_)
	return optical_depth / value_;
      return std::numeric_limits<double>::max(); 
    }
  };
  
  class basic_iop_profile {
  protected:
    pe_function profile_;
    const double epsilon_ = std::numeric_limits<double>::epsilon()*10;
  public:
    basic_iop_profile() = default;
    basic_iop_profile(const pe_function& vertical_profile)
      : profile_{vertical_profile} {}
    //    void basic_iop_profile(const std::vector<double>& heights,
    //			   const std::vector<double>& iop_values) {
    //profile_ = pe_function{heights, iop_values};
    ///}
    virtual double optical_depth(const pose& start, double distance) const = 0;
    virtual double distance(const pose& start, double optical_depth) const = 0;
    //void append(double height, double iop_value) {
    //  profile_.append(point{height,iop_value});
    //}
  };
    
  class iop_z_profile : public basic_iop_profile {
    const unit_vector z_direction_{0,0,1};
    public:
    using basic_iop_profile::basic_iop_profile;
    double optical_depth(const pose& start, double distance) const {
      double mu = zk(start);
      if (fabs(mu) < epsilon_) {
	double v = profile_.value(start.position().z());
	return constant_iop(v).optical_depth(distance);
      };
      double tau = profile_.integral(next_z(start,0),next_z(start,distance))/mu;
      return tau;
    }
    double distance(const pose& start, double optical_depth) const {
      double mu = zk(start);
      if (fabs(mu) < epsilon_) {
	double v = profile_.value(start.position().z());
	return constant_iop(v).distance(optical_depth);
      };
      std::optional<double> z =
	profile_.integral_limit_b(next_z(start,0),optical_depth * mu);
      if (z.has_value()) {
	double zr0 = dot(z_direction_, start.position());
	return (*z - zr0) / mu;
      }
      return std::numeric_limits<double>::max(); 
    }
  private:
    vector next_position(const pose& start, double distance) const {
      return start.position() + start.direction() * distance;
    }
    double next_z(const pose& start, double distance) const {
      return dot(next_position(start, distance), z_direction_);
    }
    double zk(const pose& start) const {
      return dot(z_direction_, start.direction());
    }
  };
    
    /*
    class atmospheric : public base {
      double planet_radius_;
    public:
      planetary(const pose& start_pose, pe_function vertical_profile, double planet_radius)
	: base{start_pose,vertical_profile}, planet_radius{planet_radius} {
	for (size_t i=0; i < sl.size(); ++i) {
	  slant_profile_.push_back(vertical_profile_.y[i]/r*(dot(r0,uv)+sl[i]));
	}
      }
      double optical_depth(double distance) const {
	return slant_profile(p.direction()).integral(0,distance);
      }
      double distance(double optical_depth) const {
	return slant_profile(p.direction()).integral_limit_b(0, optical_depth);
      }
    private:      
      compute_slant_profile() {
	std::vector<double> tl  = trajectory_length();
	const std::vector<double>& h = vertical_profile_.x();
	const std::vector<double>& coefficient = vertical_profile_.y();
	for (size_t i=0; i < tl.size(); ++i) {
	  double r = planet_radius_ + h[i];
	  double c = coefficient[i] / r * (dot(start_position_,uv)+tl[i]);
	  slant_profile_.insert({tl[i],c});
	}
      }
      std::vector<double> trajectory_length() {
	const std::vector<double>& h = vertical_profile_.x();
	std::vector<double> tl;
	for (size_t i=0; i < h.size(); ++i) {
	  const vector& r0v = start_pose.position();
	  const unit_vector& uv = start_pose.direction();
	  double r = planet_radius_ + h[i];
	  double r0uv = dot(r0,uv);
	  double radicand = pow(r0uv,2)+pow(r,2)-dot(r0,r0),2);
	  if (radicand > 0) {
	    double value;	    
	    if (rd > 0)
	      value = -r0uv + sqrt(radicand);
	    else
	      value = -r0uv - sqrt(radicand);
	    assert(value > 0);
	    tl.push_back(value);
	  }
	}
      assert(tl.size()>0);
      return tl;
      }
    };
    */

}

#endif