
#ifndef SCENES_H
#define SCENES_H
#include "rtweekend.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "tor.h"
#include "material.h"
//#include "moving_sphere.h"
#include "aarect.h"
#include "fractal.h"

#include <iostream>
//#include <thread>
//#include <fstream>
#include <string.h>
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
    world.add(make_shared<sphere>(point3(13, 0, 0), 2, difflight));

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
hittable_list torus()
{
    {
        hittable_list world;

        auto checker = make_shared<checker_texture>(color(0.1, 0.1, 0.1), color(1, 0.2, 0.2));
        world.add(make_shared<tor>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

        for (int a = -9; a < 9; a++)
        {
            for (int b = -9; b < 9; b++)
            {
                auto choose_mat = random_double();
                point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

                if ((center - point3(4, 0.2, 0)).length() > 0.9)
                {
                    shared_ptr<material> tor_material;

                    if (choose_mat < 0.2)
                    {
                        // diffuse
                        auto albedo = color::random() * color::random();
                        tor_material = make_shared<lambertian>(albedo);
                        world.add(make_shared<tor>(center, 0.2, tor_material));
                    }
                    else if (choose_mat < 0.6)
                    {
                        // metal
                        auto albedo = color::random(0.5, 1);
                        auto fuzz = random_double(0, 0.5);
                        tor_material = make_shared<metal>(albedo, fuzz);
                        world.add(make_shared<tor>(center, 0.2, tor_material));
                    }
                    else
                    {
                        // glass
                        tor_material = make_shared<dielectric>(1.5);
                        world.add(make_shared<tor>(center, 0.2, tor_material));
                    }
                }
            }
        }

        auto material1 = make_shared<dielectric>(1.5);
        world.add(make_shared<tor>(point3(0, 1, 0), 1.0, material1));

        auto material2 = make_shared<lambertian>(color(1, 0.6, 0.2));
        world.add(make_shared<tor>(point3(-4, 1, 0), 1.0, material2));

        auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
        world.add(make_shared<tor>(point3(4, 1, 0), 1.0, material3));

        return world;
    }
}

#endif