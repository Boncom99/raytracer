//#include "rtweekend.h"
#include "camera.h"

#include "scenes.h"

#include <iostream>
#include <thread>
#include <fstream>
#include <string.h>

color ray_color(const ray &r, const color &background, const hittable &world, int depth)
{
    hit_record rec;
    //iff we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);
    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
}

int main()
{

    // Image

    const int image_width = 200;
    double aspect_ratio;
    int samples_per_pixel;    //100
    const int max_depth = 50; //max number of bounces of a ray
    //
    //World
    hittable_list world;

    point3 lookfrom;
    point3 lookat;
    auto vfov = 40.0; //Zoom
    auto aperture = 0.0;
    color background(0, 0, 0);

    int n = 5;
    switch (0)
    {
        //default:

    case 1:
        world = random_scene();
        samples_per_pixel = 400;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        aperture = 0.1;
        break;
        //default:

    case 2:
        world = two_perlin_spheres();
        samples_per_pixel = 400;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        break;
    //default:
    case 3:
        world = simple_light();
        samples_per_pixel = 200;
        background = color(0, 0, 0);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(26, 3, 6);
        lookat = point3(0, 2, 0);
        vfov = 20.0;
        break;
    //default:
    case 4:
        world = earth_sun();
        samples_per_pixel = 200;
        background = color(0, 0, 0);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(0, 2, 6);
        lookat = point3(0, 0, 0);
        vfov = 60.0;
        break;
    //default:
    case 5:

        world = square_broken(n);
        samples_per_pixel = 100;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 3.0 / 2.0;
        lookfrom = point3(n + 2 + 3 + 4 * n, 6 + n, 4 * n);
        lookat = point3(n - 1, n - 0.75, -3);
        vfov = 20.0;
        break;
    //default:
    case 6:

        world = piramides();
        samples_per_pixel = 100;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(14, 4, 8);
        lookat = point3(1, 1, -2);
        vfov = 20.0;
        break;
    default:
    case 7:

        world = fractal();
        samples_per_pixel = 50;
        background = color(0, 0, 0);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(4, 0, 10);
        lookat = point3(4, 0, 0);
        vfov = 20.0;
        break;
    }

    // Camera

    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);
    // Render

    //multithreading

    auto f = [](const camera &cam, const int image_height, const int longitud, const int finish, const int image_width, const hittable &world, const color &background, const int samples_per_pixel, const int l) {
        std::ofstream myfile;
        myfile.open("../Desktop/" + to_string(l) + ".ppm");
        if (l == 0)
        {
            myfile << "P3\n"
                   << image_width << ' ' << image_height << "\n255\n";
        }
        for (int j = longitud; j > finish; j--)

        {
            for (int i = 0; i < image_width; ++i)
            {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s)
                {
                    auto u = (i + random_double()) / (image_width - 1);
                    auto v = (j + random_double()) / (image_height - 1);
                    ray r = cam.get_ray(u, v);
                    pixel_color += ray_color(r, background, world, max_depth);
                }
                myfile << write_color(pixel_color, samples_per_pixel);
            }
        }
        myfile.close();
    };
    int core = 1;

    if (core == 1)
    {
        f(std::ref(cam), image_height, image_height - 1, -1, image_width, std::ref(world), std::ref(background), samples_per_pixel, 0);
    }
    else
    {
        int m = image_height / core;
        int mm[core + 1];
        for (int i = 0; i < core; i++)
        {
            mm[i] = image_height - m * (i)-1;
        }
        mm[core] = -1;
        std::vector<std::thread> threads;

        for (int t = 0; t < core; t++)
            threads.emplace_back(f, std::ref(cam), image_height, mm[t], mm[t + 1], image_width, std::ref(world), std::ref(background), samples_per_pixel, t);

        for (auto &th : threads)
            th.join();

        std::ifstream file1("../Descktop/0.ppm");
        std::ifstream file2("../Descktop/1.ppm");
        std::ifstream file3("../Descktop/2.ppm");
        std::ifstream file4("../Descktop/3.ppm");
        std::ifstream file5("../Descktop/4.ppm");
        std::ifstream file6("../Descktop/5.ppm");
        std::ifstream file7("../Descktop/6.ppm");
        std::ifstream file8("../Descktop/7.ppm");
        std::ofstream FINAL("../Descktop/fractalFINAL2.ppm");
        FINAL << file1.rdbuf() << file2.rdbuf() << file3.rdbuf() << file4.rdbuf() << file5.rdbuf() << file6.rdbuf() << file7.rdbuf() << file8.rdbuf();

        remove("../Desktop/0.ppm");
        remove("../Desktop/1.ppm");
        remove("../Desktop/2.ppm");
        remove("../Desktop/3.ppm");
        remove("../Desktop/4.ppm");
        remove("../Desktop/5.ppm");
        remove("../Desktop/6.ppm");
        remove("../Desktop/7.ppm");
    }
}
