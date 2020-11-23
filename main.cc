#include "rtweekend.h"
#include "camera.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "material.h"

#include <iostream>
color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;
   //iff we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <=0)
       return color(0,0,0); 
    if (world.hit(r, 0.001, infinity, rec)) {
       ray scattered;
       color attenuation;
       if(rec.mat_ptr->scatter(r,rec,attenuation, scattered))
          return attenuation *ray_color(scattered,world,depth-1);
       return color(0,0,0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}


hittable_list random_scene(int n) {
    hittable_list world;
    shared_ptr<material> sphere_material;
      double r =500;
    auto ground_material = make_shared<lambertian>(color(1, 0.7, 0.5));
    world.add(make_shared<sphere>(point3(0,-r,0), r, ground_material));

    auto glass = make_shared<dielectric>(1.5);
    auto metall = make_shared<metal>(color(0.6, 0.6, 0.6), 0.0);

     for(int y =0;y<n; y++){
        for (int z = 0; z > -n; z--) {
          for (int x = 0; x < n; x++) {

                point3 bola(x,y+0.5,z-3);
          
              //  auto rand = random_double();
              //  sphere_material = make_shared<lambertian>(color(1,0.7,random_double() ));
             //   if(rand <=1.0/3.0){
             //   world.add(make_shared<sphere>(bola, 0.5,  glass ));
             //   }
            //    if (rand>0.5){
                world.add(make_shared<sphere>(bola, 0.5,metall));
            //    }
            //    else {
            //        world.add(make_shared<sphere>(bola, 0.5,sphere_material));
            //    }

        }
    }
    }
    
    return world;
}

int main() {

    // Image

    const auto aspect_ratio = 2.0 / 3.0;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel=500; //100
    const int max_depth=50; //max number of bounces of a ray

    // World
    int n =4;
    auto world = random_scene(n);

    // Camera
    point3 lookfrom(n-1+3 + 4*n,6+n,4*n);
    point3 lookat(n-1,n-0.5,-3);
    vec3 vup(0,1,0);
    auto dist_to_focus = 12.0;
    auto aperture = 0.0001; //small -> landscape , big ->small and close object       
       camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
     // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
           color pixel_color(0,0,0);
           for (int s =0; s< samples_per_pixel;++s){
               auto u = (i+random_double ())/ (image_width-1);
               auto v = (j+random_double ())/ (image_height-1);
               ray r=cam.get_ray(u,v);
               pixel_color +=ray_color(r,world,max_depth);

           }
            write_color(std::cout, pixel_color,samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
