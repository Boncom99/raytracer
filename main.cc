#include "rtweekend.h"
#include "camera.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "material.h"
#include "moving_sphere.h"

#include <iostream>
color ray_color(const ray &r, const hittable &world, int depth)
{
    hit_record rec;
    //iff we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);
    if (world.hit(r, 0.001, infinity, rec))
    {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0, 0, 0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

hittable_list random_scene(int n)
{
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(1, 0.2, 0.2));
    world.add(make_shared<sphere>(point3(n - 1, -800, -4), 800, ground_material));

    shared_ptr<material> sphere_material;
    sphere_material = make_shared<metal>(color(0.6, 0.6, 0.6), 0.0);
    vec3 v(static_cast<int>(-n / 2), 0.5, static_cast<int>(n / 2));

    for (int y = 0; y < n; y++)
    {
        for (int z = 0; z > -n; z--)
        {
            for (int x = 0; x < n; x++)
            {
                if (y < n - 2)
                {
                    point3 bola(x, y + 0.5, z - 3);
                    world.add(make_shared<sphere>(bola + v, 0.5, sphere_material));
                }
                else if (y == n - 2 && (x < n - 1 || z > -n + 1))
                {
                    point3 bola(x, y + 0.5, z - 3);
                    world.add(make_shared<sphere>(bola + v, 0.5, sphere_material));
                }
                else if (y == n - 1 && ((x < n - 2 || z > -n + 2)))
                {
                    point3 bola(x, y + 0.5, z - 3);
                    world.add(make_shared<sphere>(bola + v, 0.5, sphere_material));
                }
            }
        }
    }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            auto rand = random_double();
            point3 center(1.2 * i + 0.2 * rand + n - 1, 0.5, -1.2 * j - 3 - n + 1 - 0.2 * rand);
            world.add(make_shared<sphere>(center + v, 0.5, sphere_material));
        }
    }

    return world;
}

int main()
{

    // Image

    const auto aspect_ratio = 3.0 / 2.0;
    const int image_width = 400;
    const int samples_per_pixel = 20; //100
    const int max_depth = 50;         //max number of bounces of a ray

    // World
    int n = 5;
    vec3 v(static_cast<int>(-n / 2), 0.5, static_cast<int>(n / 2));
    auto world = random_scene(n);

    // Camera
    point3 lookfrom(n + 2 + 3 + 4 * n, 6 + n, 4 * n);
    point3 lookat(n - 1, n - 0.75, -3);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.01;
    int image_height = static_cast<int>(image_width / aspect_ratio);
    camera cam(lookfrom + v, lookat + v, vup, 20, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);
    // Render

    std::cout << "P3\n"
              << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
