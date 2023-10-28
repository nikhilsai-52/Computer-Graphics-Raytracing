#ifndef SPHERE_SHAPE_H
#define SPHERE_SHAPE_H

#include "rayhit.h"
#include "Vector3D.h"

class sphere_shape : public rayhit {
public:
    sphere_shape(point3D _center, double _radius) : center(_center), radius(_radius) {}

    bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const override {
        Vector3D oc = r.origin() - center;
        auto a = dot(r.direction(), r.direction());
        auto b = 2.0 * dot(oc, r.direction());
        auto c = dot(oc, oc) - radius * radius;
        auto discriminant = b * b - 4 * a * c;
        auto root = (-b - sqrt(discriminant)) / (2.0 * a);


        /*auto a = r.direction().length_squared();
        auto half_b = dot(oc, r.direction());
        auto c = oc.length_squared() - radius * radius;

        auto discriminant = half_b * half_b - a * c;*/
        if (discriminant < 0) return false;
        auto sqrtd = sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        /*auto root = (-half_b - sqrtd) / a;*/
        if (root <= ray_tmin || ray_tmax <= root) {
            auto root = (-b + sqrt(discriminant)) / (2.0 * a);
            //root = (-half_b + sqrtd) / a;
            if (root <= ray_tmin || ray_tmax <= root)
                return false;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        Vector3D outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);

        return true;
    }

private:
    point3D center;
    double radius;
};

#endif
