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
std::mutex writeM;

struct BlockJob
{
    int rowStart;
    int rowEnd;
    int colSize;
    int spp;
    std::vector<int> indices;
    std::vector<vec3> colors;
};
void CalculateColor(BlockJob job, std::vector<BlockJob> &imageBlocks, int ny, camera cam, hittable &world,
                    std::mutex &mutex, std::condition_variable &cv, std::atomic<int> &completedThreads, const color background)
{
    {
        for (int j = job.rowStart; j > job.rowEnd; j--)
        {

            for (int i = 0; i < job.colSize; i++)
            {
                vec3 col(0, 0, 0);
                for (int s = 0; s < job.spp; ++s)
                {
                    float u = float(i + random_double()) / float(job.colSize);
                    float v = float(j + random_double()) / float(ny);
                    ray r = cam.get_ray(u, v);
                    col += ray_color(r, background, world, 50); //max depth
                }
                col /= float(job.spp);
                col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

                const unsigned int index = j * job.colSize + i;
                job.indices.push_back(index);
                job.colors.push_back(col);
            }
        }

        {
            std::lock_guard<std::mutex> lock(mutex);
            imageBlocks.push_back(job);
            completedThreads++;
            cv.notify_one();
        }
    }
}
int main()
{

    // Image

    double aspect_ratio;
    int samples_per_pixel; //100
    //const int max_depth = 50; //max number of bounces of a ray

    //World
    hittable_list world;

    point3 lookfrom;
    point3 lookat;
    auto aperture = 0.0; //if diferent than 0 it blur(desenfocar) the image.  in 0, there is no blur
    auto vfov = 40.0;    //Zoom
    color background(0, 0, 0);
    int n = 5;

    switch (0)
    {
    //default:
    case 1:
        world = random_scene();
        samples_per_pixel = 10;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 4.0 / 3.0;
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        aperture = 0;
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
    default:
    case 8:
        world = torus();
        samples_per_pixel = 5;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        aperture = 0.1;
        break;
    }

    // Camera
    const int image_width = 50;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    int pixelCount = image_height * image_width;

    vec3 vup(0, 1, 0);
    auto dist_to_focus = (lookfrom - lookat).length(); //10.0;

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);

    vec3 *image = new vec3[pixelCount];
    memset(&image[0], 0, pixelCount * sizeof(vec3));

    // Render
    std::cout << "start raytracing... \n"
              << std::flush;

    auto start = std::chrono::high_resolution_clock::now();

    const int nThreads = std::thread::hardware_concurrency();
    int rowsPerThread = image_height / nThreads;
    int leftOver = image_height % nThreads;
    std::mutex mutex;
    std::condition_variable cvResults;
    std::vector<BlockJob> imageBlocks;
    std::atomic<int> completedThreads = {0};
    std::vector<std::thread> threads;
    for (int i = nThreads - 1; i >= 0; i--)
    {
        BlockJob job;
        job.rowStart = (i + 1) * rowsPerThread - 1;
        job.rowEnd = job.rowStart - rowsPerThread;
        if (i == 0)
        {
            job.rowEnd = job.rowStart - rowsPerThread - leftOver;
        }
        job.colSize = image_width;
        job.spp = samples_per_pixel;

        std::thread t([job, &imageBlocks, image_height, &cam, &world, &mutex, &cvResults, &completedThreads, &background]() {
            CalculateColor(job, imageBlocks, image_height, cam, world, mutex, cvResults, completedThreads, background);
        });
        threads.push_back(std::move(t));
    }

    // launched jobs. need to build image.
    // wait for number of jobs = pixel count
    {
        std::unique_lock<std::mutex> lock(mutex);
        cvResults.wait(lock, [&completedThreads, &nThreads] {
            return completedThreads == nThreads;
        });
    }

    for (std::thread &t : threads)
    {
        t.join();
    }

    for (BlockJob job : imageBlocks)
    {
        // int index = job.rowStart;
        int colorIndex = 0;
        for (vec3 &col : job.colors)
        {
            int colIndex = job.indices[colorIndex];
            image[colIndex] = col;
            ++colorIndex;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    int frameTimeMs = static_cast<int>(diff.count());
    std::cout << " - time " << frameTimeMs << " ms \n";
    //Open file

    std::string filename =
        "../Desktop/oout-x" + std::to_string(image_width) + "-y" + std::to_string(image_height) + "-s" + std::to_string(samples_per_pixel) + "-ms" + std::to_string(frameTimeMs) + ".ppm";

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

    for (int j = image_height - 1; j >= 0; j--)
        for (int i = 0; i < image_width; i++)
        {
            fileHandler
                << static_cast<int>(255.99f * image[j * image_width + i].e[0]) << " "
                << static_cast<int>(255.99f * image[j * image_width + i].e[1]) << " "
                << static_cast<int>(255.99f * image[j * image_width + i].e[2]) << "\n";
        }
    std::cout << "File Saved" << std::endl;
    fileHandler.close();
    delete[] image;
    return 0;
}
