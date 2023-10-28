#include <math.h>
#include <Vector3D.h>
#include <rayhit.h>

Vector3D loc = *new Vector3D(-1,2, -1);
using color = Vector3D;

class Light {


public:
    Light() {}
    
    color static illuminate(const ray& r, color& surface_col, Vector3D interection, Vector3D normal);
};

color Light::illuminate(const ray& r, color& surface_col, Vector3D interection, Vector3D normal)
{
    color setcolor;

    //ambient
    double x = surface_col.x() * 0.5;
    double y = surface_col.y() * 0.5;
    double z = surface_col.z() * 0.5;

    // diffuse
    Vector3D L = loc - interection;
    unit_vector(L);
    double dotp = dot(L,normal);
    if (dotp > 0) {
        x += dotp * surface_col.x();
        y += dotp * surface_col.y();
        z += dotp * surface_col.z();

        // specular
        Vector3D R = L - normal * (2. * dot(normal,L));
        unit_vector(R);
        double dot2 = dot(R,r.dir);
        if (dot2 > 0) {
            x += pow(dot2, 200);
            y += pow(dot2, 200);
            z += pow(dot2, 200);
        }
    }

    return color(x, y, z);
    
}



