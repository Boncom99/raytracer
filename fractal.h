#ifndef FRACTAL_H
#define FRACTAL_H
#include "rtweekend.h"
#include <iostream>

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
#endif