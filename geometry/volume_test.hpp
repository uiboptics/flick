#include "volume.hpp"

namespace flick {
namespace geometry {
  begin_test_case(volume_test_A) {
    class content
    // User defined volume content. May include material, coating,
    // etc.
    {      
    public:
      double mass{9};
    };
    cube<content> c{1};       
    auto ui = c.get_uniform_intersections(10000);
    check_close(ui.boundary_area(),6,9);
    check_close(ui.enclosed_volume(),1,9);
    check_small(rms(ui.center_of_gravity(),{0,0,0}),0.1);

    cube<content> c2{1};
    check_small(c2().mass-9,1e-9);
  
    sphere<content> small_sphere{2};
    small_sphere.insert(c2);
     
    sphere<content> large_sphere{3};
    large_sphere.insert(small_sphere);
    navigator<content> nav(large_sphere);
    nav.go_inward();
    nav.go_inward();
    nav.go_outward();
    nav.go_outward();
    nav.go_to("cube");
    nav.go_outward();
    nav.go_outward(); 
    pose observer{{0,0,-10},{0,0}};
    observer.move_to({0,0,-2.99});
    vector pos = (*nav.next_intersection(observer)).position();
    check_small(rms(pos,{0,0,-2}));
    check(nav.next_volume(observer).name()=="sphere");
    nav.go_inward();
    observer.move_to({0,0,-1.99});
    pos = (*nav.next_intersection(observer)).position();
    check_small(rms(pos,{0,0,-0.5}));
    check(nav.next_volume(observer).name()=="cube");
    nav.go_inward();
    check_throw(nav.go_inward());
    nav.go_outward();
    nav.go_outward();
    check_throw(nav.go_outward());
    
    sphere sph_copy = large_sphere;
    vector new_pos = {1,1,1};
    sph_copy.move_by(new_pos);
    navigator nav2(sph_copy);
    nav2.go_inward();
    nav2.go_inward();
    vector sph_pos = nav2.go_outward().placement().position();
    check_small(rms(sph_pos,new_pos));

    cube<content> c3{0.1};
    c3.move_by({1,0,0});
    c3.rotate_by(rotation_about_z(constants::pi/2),{-1,0,0});
    check_small(rms(c3.placement().position(),{-1,2,0}));
    
  } end_test_case()

  begin_test_case(volume_test_B) {
    class content {      
    public:
      double mass{9};
    };
    using semi_infinite_box = semi_infinite_box<content>;
    semi_infinite_box space;
    semi_infinite_box atmosphere;
    semi_infinite_box bottom;

    space.name("space");
    atmosphere.name("atmosphere");
    bottom.name("bottom");
    
    space.move_by({0,0,10});
    atmosphere.move_by({0,0,1});
    atmosphere.insert(bottom);
    space.insert(atmosphere);

    auto nav = navigator<content>(space);
    pose observer1{{0,0,5},{0,0,-1}};
    nav.go_to(nav.next_volume(observer1));
    pose observer2{{0,0,2},{0,0,1}};
    nav.go_to(nav.next_volume(observer2));

    nav.go_to("bottom");
    check(nav.find("atmosphere").name()=="atmosphere","find error");
    check_throw(nav.find("wrongname"));
  } end_test_case()
}
}
