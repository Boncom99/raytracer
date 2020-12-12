#ifndef RAY_H
#define RAY_H
#include "vec3.h"

class ray
{
public:
	ray() {}															//constructor if there is no argument given
	ray(const point3 &origin, const vec3 &direction, double time = 0.0) //constructor if there IS argument ;
		: orig(origin), dir(direction), tm(time)						//initialization on variables ray.orig = origin, ray.dir= direction,...
	{
	}

	point3 origin() const { return orig; }
	vec3 direction() const { return dir; }
	double time() const { return tm; }

	point3 at(double t) const
	{
		return orig + t * dir;
	}
	point3 at_unit(double t) const
	{
		return orig + t * unit_vector(dir);
	}

public:
	point3 orig;
	vec3 dir;
	double tm;
};
#endif
