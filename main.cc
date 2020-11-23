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


hittable_list random_scene() {
    hittable_list world;
    shared_ptr<material> sphere_material;
      double r =10000;
    auto ground_material = make_shared<lambertian>(color(1, 0.2, 0.2));
    world.add(make_shared<sphere>(point3(0,-r,0), r, ground_material));

    point3 bola(1,1+sqrt(2),-2);
    vec3 vbola1(-1,-sqrt(2),+1);
    vec3 vbola2(1,-sqrt(2),+1);
    vec3 vbola3(1,-sqrt(2),-1);
    vec3 vbola4(-1,-sqrt(2),-1);

    point3 bola1= bola+ vbola1;
    point3 bola2= bola+ vbola2;
    point3 bola3= bola+ vbola3;
    point3 bola4= bola+ vbola4;

    auto glass = make_shared<dielectric>(1.5);
    auto metall = make_shared<metal>(color(0.7, 0.7, 0.5), 0.0);

    for (int a = -2; a < 2; a++) {
        for (int b = -2; b < 2; b++) {

                //point3 bola(9*a + 0.9*random_double(), 1+sqrt(2) -(r-sqrt(abs(r^2-z^2))),z);
                point3 bola(9*a + 0.9*random_double(), 1+sqrt(2) ,9*b+ 0.9*random_double());
                  point3 bola1= bola+ vbola1;
                  point3 bola2= bola+ vbola2;
                  point3 bola3= bola+ vbola3;
                  point3 bola4= bola+ vbola4;
          
                auto rand = random_double();
                sphere_material = make_shared<lambertian>(color(1,0.7,random_double() ));
                if(rand <=1.0/3.0){
                world.add(make_shared<sphere>(bola1, 1.0,  glass ));
                }
                else if (rand>1.0/3.0 && rand<=2.0/3.0){
                world.add(make_shared<sphere>(bola1, 1.0,metall));
                }
                else {
                    world.add(make_shared<sphere>(bola1, 1.0,sphere_material));
                }

                rand = random_double();
                if(rand <=1.0/3.0){
                world.add(make_shared<sphere>(bola2, 1.0,  glass ));
                }
                else if (rand>1.0/3.0 && rand<=2.0/3.0){
                world.add(make_shared<sphere>(bola2, 1.0,metall));
                }
                else {
                    world.add(make_shared<sphere>(bola2, 1.0,sphere_material));
                }
 
                rand = random_double();
                sphere_material = make_shared<lambertian>(color(1,0.7,random_double() ));
                if(rand <=1.0/3.0){
                world.add(make_shared<sphere>(bola3, 1.0,  glass ));
                }
                else if (rand>1.0/3.0 && rand<=2.0/3.0){
                world.add(make_shared<sphere>(bola3, 1.0,metall));
                }
                else {
                    world.add(make_shared<sphere>(bola3, 1.0,sphere_material));
                }
                rand = random_double();
                if(rand <=1.0/3.0){
                world.add(make_shared<sphere>(bola4, 1.0,  glass ));
                }
                else if (rand>1.0/3.0 && rand<=2.0/3.0){
                world.add(make_shared<sphere>(bola4, 1.0,metall));
                }
                else {
                    world.add(make_shared<sphere>(bola4, 1.0,sphere_material));
                }
                rand = random_double();
                if(rand <=1.0/3.0){
                world.add(make_shared<sphere>(bola, 1.0,  glass ));
                }
                else if (rand>1.0/3.0 && rand<=2.0/3.0){
                world.add(make_shared<sphere>(bola, 1.0,metall));
                }
                else {
                    world.add(make_shared<sphere>(bola, 1.0,sphere_material));
                }

 
 
 

                    }
    }
    
    return world;
}

int main() {

    // Image

    const auto aspect_ratio = 3.0 / 2.0;
    const int image_width = 800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel=300; //100
    const int max_depth=50; //max number of bounces of a ray

    // World
    auto world = random_scene();

    // Camera
    point3 lookfrom(14,4,8);
    point3 lookat(1,1,-2);
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
