#include <string>
#include <float.h>
#include <fstream>
#include <chrono>
#include <thread>
#include <vector>
#include <mutex>
#include <atomic>
#include "camera.h"

#include "scenes.h"

color ray_color(const ray &r, const color &background, const hittable &world, int depth)
{
    hit_record rec;
    //iff we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);
    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.01, infinity, rec))
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

    double aspect_ratio;
    int samples_per_pixel;    //100
    const int max_depth = 50; //max number of bounces of a ray

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
    default:
    case 1:
        world = random_scene();
        samples_per_pixel = 10;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 4.0 / 3.0;
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
    //default:
    case 7:

        world = fractal();
        samples_per_pixel = 50;
        background = color(0, 0, 0);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(4, 0, 10);
        lookat = point3(4, 0, 0);
        vfov = 20.0;
        break;
    //default:
    case 8:
        world = torus();
        samples_per_pixel = 20;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        aperture = 0.1;
        break;
    }

    // Camera
    const int image_width = 800;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    int pixelCount = image_height * image_width;

    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);
    //
    vec3 *image = new vec3[pixelCount];
    memset(&image[0], 0, pixelCount * sizeof(vec3));

    // Render
    auto start = std::chrono::high_resolution_clock::now();
    // std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; j--)
    //for (int j = 0; j < image_height; j++)

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
                pixel_color += ray_color(r, background, world, max_depth);
            }
            pixel_color /= float(samples_per_pixel);
            pixel_color = vec3(sqrt(pixel_color[0]), sqrt(pixel_color[1]), sqrt(pixel_color[2]));
            const unsigned int index = (image_height - j - 1) * image_width + i;
            image[index] = pixel_color;

            //write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    int frameTimeMs = static_cast<int>(diff.count());
    std::cout << " - time " << frameTimeMs << " ms \n";
    //Open file

    std::string filename =
        "../Desktop/out-x" + std::to_string(image_width) + "-y" + std::to_string(image_height) + "-s" + std::to_string(samples_per_pixel) + "-ms" + std::to_string(frameTimeMs) + ".ppm";

    std::ofstream fileHandler;
    fileHandler.open(filename, std::ios::out | std::ios::binary);
    if (!fileHandler.is_open())
    {
        std::cout << "error opening file" << std::endl;
        return -1;
    }
    else
    {
        std::cout << "files opened" << std::endl;
    }

    fileHandler << "P3\n"
                << image_width << " " << image_height << "\n255\n";

    //write color
    for (unsigned int i = 0; i < pixelCount; ++i)
    {
        fileHandler
            << static_cast<int>(255.99f * image[i].e[2]) << " "
            << static_cast<int>(255.99f * image[i].e[1]) << " "
            << static_cast<int>(255.99f * image[i].e[0]) << "\n";
    }
    std::cout << "File Saved" << std::endl;
    fileHandler.close();
    delete[] image;
}
