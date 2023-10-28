#ifndef RAYHIT_H
#define RAYHIT_H
#include <memory>
#include "Vector3D.h"

class material;

class ray {
public:
    ray() {}

    ray(const point3D& origin, const Vector3D& direction) : orig(origin), dir(direction) {}

    point3D origin() const { return orig; }
    Vector3D direction() const { return dir; }

    point3D at(double t) const {
        return orig + t * dir;
    }
//previously private
public: 
    point3D orig;
    Vector3D dir;
};

class hit_record {
public:
    point3D p;
    Vector3D normal;
    std::shared_ptr<material> mat_ptr;
    double t;
    
    bool front_face;

    inline void set_face_normal(const ray& r, const Vector3D& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0; // If ray and outward normal are less than 90 degrees apart
        normal = front_face ? outward_normal : -outward_normal; // If its the front face then leave the normal outward, else have it point inward
    }
};

class  rayhit{
public:
    virtual ~rayhit() = default;

    virtual bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const = 0;
};

#endif
