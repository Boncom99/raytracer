#include "rtweekend.h"
#include "camera.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "material.h"
#include "moving_sphere.h"
#include "aarect.h"

#include <iostream>
#include <thread>
#include <fstream>
void triangle(double a, double R, double &x, double &y, int i)
{
    double A = 0, b = 0, c = 0;
    //pow(1.0 / a, i)
    A = R * (1.0 + a * a);
    b = R * a * (1 + a);
    c = R * (1.0 + a);
    y = 0.5 * sqrt((A + b + c) * (-A + b + c) * (A - b + c) * (A + b - c)) / b;
    x = sqrt(c * c - y * y);
}
void newpoint(double a, double R, int i, double &p0, double &p1, vec3 u, int ident)
{
    double x;
    double y;
    triangle(a, R, x, y, i);
    vec3 v = perpen(u);

    if (ident == 1)
    {
        p0 = p0 - x * u[0] + y * v[0];
        p1 = p1 - x * u[1] + y * v[1];
    }
    else if (ident == 2)
    {
        p0 = p0 - x * u[0] - y * v[0];
        p1 = p1 - x * u[1] - y * v[1];
    }
    else
    { //canviar. en comptes d'allargar el vector u, hem de fer la ¿suma? dels dos
        p0 = p0 + (R + R * a) * u[0];
        p1 = p1 + (R + R * a) * u[1];
    }
}

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
hittable_list fractal()
{
    hittable_list world;

    double r = 3.0;
    double a = 1.0 / 3.0;

    auto material1 = make_shared<lambertian>(color(1, 0.2, 0.6));
    auto metall = make_shared<metal>(color(0.7, 0.7, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(0, 0, 0), r, material1));
    world.add(make_shared<sphere>(point3(r * (1 + a), 0, 0), r * a, material1));
    auto difflight = make_shared<diffuse_light>(color(10, 10, 10));
    // world.add(make_shared<sphere>(point3(17, 3, -2.5), 2.5, difflight));
    world.add(make_shared<sphere>(point3(13, 0, 1), 2, difflight));

    static int n = 4;
    static int L = 2000;
    double R = r * a;

    int j1 = 1;
    int j2 = 2;
    int j3 = 3;

    double p[L][4];

    vec3 u;

    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            p[i][j] = 0;
        }
    }
    p[1][0] = r + a * r;
    p[2][0] = r + a * r;
    p[3][0] = r + a * r;
    p[0][0] = r + a * r;

    int q1 = 0, q2 = 0, q3 = 0;

    for (int i = 0; i < n; i++)
    {

        for (int k = 0; k < pow(3, i); k++)
        {

            u = vec3(p[j1][0] - p[j1][2], p[j1][1] - p[j1][3], 0);

            u = unit_vector(u);

            p[j2][0] = p[j1][0];
            p[j2][1] = p[j1][1];
            p[j2][2] = p[j1][2];
            p[j2][3] = p[j1][3];
            p[j3][0] = p[j1][0];
            p[j3][1] = p[j1][1];
            p[j3][2] = p[j1][2];
            p[j3][3] = p[j1][3];
            double w1 = p[j1][0];
            double w2 = p[j1][1];
            newpoint(a, R, i, p[j1][0], p[j1][1], u, 1);
            newpoint(a, R, i, p[j2][0], p[j2][1], u, 2);
            newpoint(a, R, i, p[j3][0], p[j3][1], u, 3);

            if (i == 0)
            {
                q1 = 4;
                q2 = 7;
                q3 = 10;
            }
            else
            {

                q1 = j1 + (pow(3, i) - k - 1) * 3 + k * 9 + 3;
                q2 = j1 + (pow(3, i) - k - 1) * 3 + k * 9 + 3 + 3 * 1;
                q3 = j1 + (pow(3, i) - k - 1) * 3 + k * 9 + 3 + 3 * 2;
            }
            p[q1][0] = p[j1][0];
            p[q1][1] = p[j1][1];
            p[q1][2] = w1;
            p[q1][3] = w2;
            p[q2][0] = p[j2][0];
            p[q2][1] = p[j2][1];
            p[q2][2] = w1;
            p[q2][3] = w2;
            p[q3][0] = p[j3][0];
            p[q3][1] = p[j3][1];
            p[q3][2] = w1;
            p[q3][3] = w2;

            world.add(make_shared<sphere>(point3(p[j1][0], p[j1][1], 0), R * a, material1));
            world.add(make_shared<sphere>(point3(p[j2][0], p[j2][1], 0), R * a, material1));
            world.add(make_shared<sphere>(point3(p[j3][0], p[j3][1], 0), R * a, material1));
            j1 += 3; //numero de cada punt
            j2 += 3;
            j3 += 3;
        }
        R *= a;
    }

    return world;
}
hittable_list piramides()
{
    hittable_list world;
    shared_ptr<material> sphere_material;
    double r = 10000;
    auto ground_material = make_shared<lambertian>(color(1, 0.2, 0.2));
    world.add(make_shared<sphere>(point3(0, -r, 0), r, ground_material));

    vec3 vbola1(-1, -sqrt(2), +1);
    vec3 vbola2(1, -sqrt(2), +1);
    vec3 vbola3(1, -sqrt(2), -1);
    vec3 vbola4(-1, -sqrt(2), -1);

    auto glass = make_shared<dielectric>(1.5);
    auto metall = make_shared<metal>(color(0.7, 0.7, 0.5), 0.0);

    for (int a = -2; a < 2; a++)
    {
        for (int b = -2; b < 2; b++)
        {

            //point3 bola(9*a + 0.9*random_double(), 1+sqrt(2) -(r-sqrt(abs(r^2-z^2))),z);
            point3 bola(9 * a + 0.9 * random_double(), 1 + sqrt(2), 9 * b + 0.9 * random_double());
            point3 bola1 = bola + vbola1;
            point3 bola2 = bola + vbola2;
            point3 bola3 = bola + vbola3;
            point3 bola4 = bola + vbola4;

            auto rand = random_double();
            sphere_material = make_shared<lambertian>(color(1, 0.7, random_double()));
            if (rand <= 1.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola1, 1.0, glass));
            }
            else if (rand > 1.0 / 3.0 && rand <= 2.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola1, 1.0, metall));
            }
            else
            {
                world.add(make_shared<sphere>(bola1, 1.0, sphere_material));
            }

            rand = random_double();
            if (rand <= 1.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola2, 1.0, glass));
            }
            else if (rand > 1.0 / 3.0 && rand <= 2.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola2, 1.0, metall));
            }
            else
            {
                world.add(make_shared<sphere>(bola2, 1.0, sphere_material));
            }

            rand = random_double();
            sphere_material = make_shared<lambertian>(color(1, 0.7, random_double()));
            if (rand <= 1.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola3, 1.0, glass));
            }
            else if (rand > 1.0 / 3.0 && rand <= 2.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola3, 1.0, metall));
            }
            else
            {
                world.add(make_shared<sphere>(bola3, 1.0, sphere_material));
            }
            rand = random_double();
            if (rand <= 1.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola4, 1.0, glass));
            }
            else if (rand > 1.0 / 3.0 && rand <= 2.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola4, 1.0, metall));
            }
            else
            {
                world.add(make_shared<sphere>(bola4, 1.0, sphere_material));
            }
            rand = random_double();
            if (rand <= 1.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola, 1.0, glass));
            }
            else if (rand > 1.0 / 3.0 && rand <= 2.0 / 3.0)
            {
                world.add(make_shared<sphere>(bola, 1.0, metall));
            }
            else
            {
                world.add(make_shared<sphere>(bola, 1.0, sphere_material));
            }
        }
    }

    return world;
}
hittable_list square_broken(int n)
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
hittable_list earth_sun()
{
    hittable_list objects;
    //hearth
    auto earth_texture = make_shared<image_texture>("images/earth2.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    objects.add(make_shared<sphere>(point3(0, 0, 0), 2, earth_surface));
    //lights

    //moon

    //     auto difflight2 = make_shared<diffuse_light>(color(0.5,0.5,0.0));
    //    objects.add(make_shared<sphere>(point3(4,2.3,0), 2*0.3, difflight2));
    auto pertext = make_shared<moon_texture>(4);
    objects.add(make_shared<sphere>(point3(4, 2.3, -2), 2 * 0.30, make_shared<lambertian>(pertext)));
    //sun
    auto difflight = make_shared<diffuse_light>(color(30, 30, 20));
    objects.add(make_shared<sphere>(point3(-15, 14.89, 3), 3, difflight));

    return objects;
}
hittable_list simple_light()
{
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    auto difflight = make_shared<diffuse_light>(color(16, 16, 15));
    objects.add(make_shared<sphere>(point3(0, 5, 0), 0.5, difflight));

    return objects;
}
hittable_list random_scene()
{
    hittable_list world;

    auto checker = make_shared<checker_texture>(color(0.1, 0.1, 0.1), color(1, 0.2, 0.2));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

    for (int a = -9; a < 9; a++)
    {
        for (int b = -9; b < 9; b++)
        {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9)
            {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.2)
                {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.6)
                {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else
                {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(1, 0.6, 0.2));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list two_perlin_spheres()
{
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    return objects;
}

int main()
{

    // Image

    const int image_width = 500;
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
        samples_per_pixel = 200;
        background = color(0.70, 0.80, 1.00);
        aspect_ratio = 16.0 / 9.0;
        lookfrom = point3(14, 4, 8);
        lookat = point3(1, 1, -2);
        vfov = 20.0;
        break;
    default:
    case 7:

        world = fractal();
        samples_per_pixel = 500;
        //background = color(0, 0, 0);
        background = color(0.3, 0.7, 0.8);
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
    int core = 8;
    int m = image_height / core;
    std::cerr << "image_height = " << image_height << std::endl;
    std::cerr << "m =  " << m << std::endl;
    std::cerr << "m*8=" << m * 8 << std::endl;
    int mm[core + 1];
    for (int i = 0; i < core; i++)
    {
        mm[i] = image_height - m * (i)-1;
        std::cerr << "mm[" << i << "]=" << mm[i] << std::endl;
    }
    mm[core] = -1;
    for (int i = 0; i < core; i++)
    {
        std::cerr << "mm[" << i << "]-m[" << i + 1 << "]=" << mm[i] - mm[i + 1] << std::endl;
    }

    auto f = [](camera &cam, int image_height, int longitud, int finish, int image_width, const hittable &world, const color &background, int samples_per_pixel, const std::string A, bool l) {
        std::ofstream myfile;
        myfile.open("../" + A + ".ppm");
        if (l)
        {
            myfile << "P3\n"
                   << image_width << ' ' << image_height << "\n255\n";
        }
        // myfile << "document " << A;
        for (int j = longitud; j > finish; j--)

        {
            // std::cerr << "j=" << j << std::endl;
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
                //std::cerr << pixel_color << std::endl;
                myfile << write_color(pixel_color, samples_per_pixel);
            }
        }
        //myfile.close();
    };

    std::thread A(f, std::ref(cam), image_height, mm[0], mm[1], image_width, std::ref(world), std::ref(background), samples_per_pixel, "A", 1);
    std::thread B(f, std::ref(cam), image_height, mm[1], mm[2], image_width, std::ref(world), std::ref(background), samples_per_pixel, "B", 0);
    std::thread C(f, std::ref(cam), image_height, mm[2], mm[3], image_width, std::ref(world), std::ref(background), samples_per_pixel, "C", 0);
    std::thread D(f, std::ref(cam), image_height, mm[3], mm[4], image_width, std::ref(world), std::ref(background), samples_per_pixel, "D", 0);
    std::thread E(f, std::ref(cam), image_height, mm[4], mm[5], image_width, std::ref(world), std::ref(background), samples_per_pixel, "E", 0);
    std::thread F(f, std::ref(cam), image_height, mm[5], mm[6], image_width, std::ref(world), std::ref(background), samples_per_pixel, "F", 0);
    std::thread G(f, std::ref(cam), image_height, mm[6], mm[7], image_width, std::ref(world), std::ref(background), samples_per_pixel, "G", 0);
    std::thread H(f, std::ref(cam), image_height, mm[7], mm[8], image_width, std::ref(world), std::ref(background), samples_per_pixel, "H", 0);
    A.join();
    B.join();
    C.join();
    D.join();
    E.join();
    F.join();
    G.join();
    H.join();
    std::ifstream file1("../A.ppm");
    std::ifstream file2("../B.ppm");
    std::ifstream file3("../C.ppm");
    std::ifstream file4("../D.ppm");
    std::ifstream file5("../E.ppm");
    std::ifstream file6("../F.ppm");
    std::ifstream file7("../G.ppm");
    std::ifstream file8("../H.ppm");
    std::ofstream FINAL("../fractalFINAL.ppm");
    FINAL << file1.rdbuf() << file2.rdbuf() << file3.rdbuf() << file4.rdbuf() << file5.rdbuf() << file6.rdbuf() << file7.rdbuf() << file8.rdbuf();
}
