#ifndef TOR_H
#define TOR_H

#include "hittable.h"
#include "vec3.h"
#include "aabb.h"

class tor : public hittable
{
public:
    tor() {}
    tor(point3 cen, double r1, /*double r2,*/ shared_ptr<material> m)
        : center(cen), radius(r1), /*radius2(r2),*/ mat_ptr(m){};

    virtual bool hit(
        const ray &r, double tmin, double tmax, hit_record &rec) const override;
    virtual bool bounding_box(double time0, double time1, aabb &output_box) const override;

public:
    point3 center;
    double radius;
    //double radius2;
    shared_ptr<material> mat_ptr;

private:
    static void get_tor_uv(const point3 &p, double &u, double &v)
    {
        // p: a given point on the tor of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        auto theta = acos(-p.y());
        auto phi = atan2(-p.z(), p.x()) + pi;

        u = phi / (2 * pi);
        v = theta / pi;
    }
};

bool tor::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{

    double t = t_min;
    double dist = 0;
    //std::cerr << "t_min= " << t_min << std::endl;
    // std::cerr << "t_max= " << t_max << std::endl;
    // std::cerr << "...............t= " << t << std::endl;

    while (t >= t_min && t < t_max) //TODO
    {
        dist = distance(center, radius, r, t);

        if (dist <= 0.00001)
        {
            // std::cerr << ".............................DENTROOOO................" << std::endl;

            rec.t = t;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            rec.set_face_normal(r, outward_normal);
            get_tor_uv(outward_normal, rec.u, rec.v);
            rec.mat_ptr = mat_ptr;

            return true;
        }
        t += dist;
        // std::cerr << "t= " << t << std::endl;
    }

    return false;
}
bool tor::bounding_box(double time0, double time1, aabb &output_box) const //TODO
{
    output_box = aabb(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius));
    return true;
}
#endif
