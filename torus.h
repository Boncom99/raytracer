#ifndef TORUS_H
#define TORUS_H

#include "hittable.h"
#include "vec3.h"
#include "aabb.h"
#include <complex>

using namespace std;

class torus : public hittable
{
public:
    torus() {}
    torus(double radius1, double radius2, shared_ptr<material> m)
        : R(radius1), r1(radius2), mat_ptr(m){};

    virtual bool hit(
        const ray &r, double tmin, double tmax, hit_record &rec) const override;
    virtual bool bounding_box(double time0, double time1, aabb &output_box) const override;

public:
    double r1;
    double R;
    shared_ptr<material> mat_ptr;
};

bool torus::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    auto root = 0;
    vec3 d = r.direction();
    vec3 o = r.origin();

    auto sum_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    auto od = o[0] * d[0] + o[1] * d[1] + o[2] * d[2];
    auto o_sq = o[0] * o[0] + o[1] * o[1] + o[2] * o[2];
    auto o_minus = o_sq - (r1 * r1 + R * R);
    auto a = pow(sum_sq, 2);
    auto b = 4 * sum_sq * od;
    auto c = 2 * sum_sq * (o_minus) + 4 * pow(od, 2) + 4 * R * R * d[1] * d[1];
    auto d = 4 * (o_minus)*od + 8 * R * R * o[1] * d[1];
    auto e = pow(o_minus, 2) - 4 * R * R * (r1 * r1 - o[1] * o[1]);

    /*
    auto a = r.direction().length_squared();
    auto half_b = dot(o, r.direction());
    auto c = oc.length_squared() - radius * radius;
   
    auto discriminant = half_b * half_b - a * c; */
    if (root < 0)
        return false;

    //find the nearest root that lies in the acceptable range.
    if (root < t_min || t_max < root)
    {
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    //get_torus_uv(outward_normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;

    return true;
}

#endif
