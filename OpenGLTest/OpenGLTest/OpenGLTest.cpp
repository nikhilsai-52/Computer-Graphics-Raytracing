// Based on templates from learnopengl.com
#include <glew.h>
#include <glfw3.h>
#include <Vector3D.h>
#include<vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <random>
#include<string>
#include <rayhit.h>
#include<algorithm>
#include<Light.h>
#define GLFW_INCLUDE_NONE
#include <GL/gl.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


#include <iostream>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 800;

using namespace std;


const char* vertexShaderSource = "#version 330 core\n"
"layout (location = 0) in vec3 aPos;\n"
"layout (location = 1) in vec3 aColor;\n"
"layout (location = 2) in vec2 aTexCoord;\n"
"out vec3 ourColor;\n"
"out vec2 TexCoord;\n"
"void main()\n"
"{\n"
"gl_Position = vec4(aPos, 1.0);\n"
"ourColor = aColor;\n"
"TexCoord = vec2(aTexCoord.x, aTexCoord.y);\n"
"}\0";

const char* fragmentShaderSource = "#version 330 core\n"
"out vec4 FragColor;\n"
"in vec3 ourColor;\n"
"in vec2 TexCoord;\n"
"uniform sampler2D texture1;\n"
"void main()\n"
"{\n"
"   FragColor = texture(texture1, TexCoord);\n"
"}\n\0";

using color = Vector3D;
//const double EPSILON = 0.0001; del
constexpr double pi = 3.14159265358979323846;


double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}


double write_colorR(std::ostream& out, color pixel_color) {
    // Write the translated [0,255] value of each color component.
    return static_cast<int>(255.999 * pixel_color.x());

}

double write_colorG(std::ostream& out, color pixel_color) {
    // Write the translated [0,255] value of each color component.
    return static_cast<int>(255.999 * pixel_color.y());

}

double write_colorB(std::ostream& out, color pixel_color) {
    // Write the translated [0,255] value of each color component.
    return static_cast<int>(255.999 * pixel_color.z());

}

class Plane {
public:
    point3D point;  // A point on the plane
    Vector3D normal; // Normal of the plane

    Plane(const point3D& p, const Vector3D& n) : point(p), normal(n) {}

    // Check if a ray intersects with this plane
    bool intersects(const ray& r, float& t) const {
        float denom = dot(normal, r.direction());
        if (abs(denom) > 1e-6) {
            Vector3D p0l0 = point - r.origin();
            t = dot(p0l0, normal) / denom;
            return (t >= 0);
        }
        return false;
    }
};

class Camera {
public:
    Camera(point3D lookfrom, point3D lookat, Vector3D vup, double vfov, double aspect_ratio) {
        double theta = degrees_to_radians(vfov);
        double h = tan(theta / 2);
        double viewport_height = 2.0 * h;
        double viewport_width = aspect_ratio * viewport_height;

        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        origin = lookfrom;
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - w;
    }

    ray get_ray(double s, double t) const {
        return ray(origin, lower_left_corner + s * horizontal + t * vertical - origin);
    }

private:
    point3D origin;
    point3D lower_left_corner;
    Vector3D horizontal;
    Vector3D vertical;
    Vector3D u, v, w;
};

double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

void saveImage(const char* filename, const unsigned char* image, int width, int height) {
    // Save the image in PNG format (you can also save in BMP, TGA, or HDR with corresponding functions)
    stbi_write_png(filename, width, height, 3, image, width * 3);
}





double hit_sphere(const point3D& center, double radius, const ray& r) {
    Vector3D oc = center - r.origin();
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius * radius;
    auto discriminant = b * b - 4 * a * c;
    auto root = (-b - sqrt(discriminant)) / (2.0 * a);
    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return root;
    }

}

bool hit_plane(const point3D& pointOnPlane, const Vector3D& normal, const ray& r, double& t) {
    auto denominator = dot(normal, r.direction());
    if (denominator > 1e-6) { // Check for non-parallel ray and plane
        auto p0l0 = pointOnPlane - r.origin();
        t = dot(p0l0, normal) / denominator;
        return (t >= 0);
    }
    return false;
}



double length(const Vector3D& v) {
    return sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}

color ray_color(const ray& r, const Plane& groundPlane, int depth) {

    if (depth <= 0) return color(0, 0, 0);

    // two spheres

    //Sphere 1 dimensions
    point3D sphere_center1(-0.5, 1, 2);
    double sphere_radius1 = 0.5;
    
    // Sphere 2 dimensions
    point3D sphere_center3(-0.8, 0, 5);
    double sphere_radius3 = 1;

    //hit Sphere 1
    auto t1 = hit_sphere(sphere_center1, sphere_radius1, r);
    //hit Sphere 2
    auto t3 = hit_sphere(sphere_center3, sphere_radius3, r);

    point3D lightPosition(-0.6, 2.5, 4);


    if (t1 > 0) {


        point3D cameracenter(0, 0, 0);

        Vector3D N = unit_vector(r.at(t1) - sphere_center1);
        Vector3D lightDir = unit_vector(lightPosition - r.at(t1));


        Vector3D call = *new Vector3D;


        // Calculate ambient, diffuse, and specular components
        color ambientColor = color(1, 0, 0); // Ambient color 
        color diffuseColor = color(1, 0, 0); // Diffuse color 
        color specularColor = color(1, 1, 1); // Specular color 

        double ambientIntensity = 0.1; 
        double diffuseIntensity = 0.4;
        double specularIntensity = 0.5;

        double diffuseFactor = std::max(dot(N, lightDir), 0.0);
        color diffuseContribution = diffuseIntensity * diffuseFactor * diffuseColor;

        Vector3D viewDir = unit_vector(cameracenter - r.at(t1));
        Vector3D reflectDir = reflect(-lightDir, N);
        double specularFactor = std::pow(std::max(dot(viewDir, reflectDir), 0.0), 100); // Shininess factor 
        color specularContribution = specularIntensity * specularFactor * specularColor;

        color shading = ambientColor * ambientIntensity + diffuseContribution + specularContribution;
        return shading;


    }


    if (t3 > 0) {

        point3D cameracenter(0, 0, 0);

        Vector3D N = unit_vector(r.at(t3) - sphere_center3);
        Vector3D lightDir = unit_vector(lightPosition - r.at(t3));


        Vector3D call = *new Vector3D;


        // Calculate ambient, diffuse, and specular components
        color ambientColor = color(0, 1, 0); // Ambient color 
        color diffuseColor = color(0, 1, 0); // Diffuse color 
        color specularColor = color(1, 1, 1); // Specular color 

        double ambientIntensity = 0.1; 
        double diffuseIntensity = 0.4;
        double specularIntensity = 0.8;

        double diffuseFactor = std::max(dot(N, lightDir), 0.0);
        color diffuseContribution = diffuseIntensity * diffuseFactor * diffuseColor;

        Vector3D viewDir = unit_vector(cameracenter - r.at(t3));
        Vector3D reflectDir = reflect(-lightDir, N);
        double specularFactor = std::pow(std::max(dot(viewDir, reflectDir), 0.0), 64); // Shininess factor 
        color specularContribution = specularIntensity * specularFactor * specularColor;

        color shading = ambientColor * ambientIntensity + diffuseContribution + specularContribution;
        return shading;


    }

    //Check for intersection with glazed plane
    float t_plane;
    if (groundPlane.intersects(r, t_plane)) {
        point3D intersectPoint = r.origin() + r.direction() * t_plane;
        Vector3D normal = unit_vector(groundPlane.normal);  

        // Calculate reflection ray

        Vector3D incident_dir = unit_vector(r.direction());
        Vector3D reflected_dir = reflect((r.direction()), normal);
        ray reflected_ray(intersectPoint * 0.001 * reflected_dir, reflected_dir); // to avoid self-intersection

        // Reflectivity of the plane
        float reflectivity = 0.8; 

        // Recursively cast the reflection ray
        color reflected_color = ray_color(reflected_ray, groundPlane, depth - 1);

        // Compute plane's inherent color ( checkerboard pattern)
        color plane_color;
        if (static_cast<int>(intersectPoint.y() + intersectPoint.z()) % 2 == 0) {
            plane_color = color(1, 1, 1);  // white
        }
        else {
            plane_color = color(0, 0, 0);  // black
        }

        return reflected_color * reflectivity + plane_color * (1.0 - reflectivity);
    }

    //std::cout << "No Sphere Hit.\n";
    return color(0, 0, 0); // White color for no hit


}


int main()
{


    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef APPLE
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif


    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Display RGB Array", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);



    // // GLEW: load all OpenGL function pointers
    glewInit();

    // build and compile the shaders
    // ------------------------------------
    // vertex shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    // check for shader compile errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // fragment shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    // check for shader compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // link shaders
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    // check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);


    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    float vertices[] = {
        // positions          // colors           // texture coords
         0.5f,  0.5f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // top right
         0.5f, -0.5f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // bottom right
        -0.5f, -0.5f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // bottom left
        -0.5f,  0.5f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // top left 
    };
    unsigned int indices[] = {
        0, 1, 3, // first triangle
        1, 2, 3  // second triangle
    };
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);


    // load and create a texture 
    // -------------------------
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    // Create the image (RGB Array) to be displayed
    const int width = 512; // keep it in powers of 2!
    const int height = 512; // keep it in powers of 2!

    

    int numFrames = 300; // Noumber of frames for the animation
    double stepSize = 1.2; //  amt camera will move on each frame
    
    // User input for perspective camera
    bool usePerspectiveGeometry = false;
    double distToImage = 2;



    std::cout << "Press 'p' to switch to perspective geometry or any other key for orthographic: ";
    char input;
    std::cin >> input;
    if (input == 'p' || input == 'P') {
        usePerspectiveGeometry = true;
    }

    
    
    //light
    Light light = *new Light;


    //camera
    double focal_length = 1.0; //change the variable
    double viewh = 2.0; // view auto height
    double vieww = viewh * ((double)(width) / height); //viewport width
    Vector3D camera_center(0.0, 0.0, 0.0);

    //vectors across view
    auto vu = Vector3D(vieww, 0, 0);
    auto vv = Vector3D(0, -viewh, 0);

    //vector diff
    auto diff_u = vu / width;
    auto diff_v = vv / height;

    Plane groundPlane(point3D(-5, 0, 0), Vector3D(1, 0, 0));  // plane lying on the XZ plane, 1 unit below the origin


    auto v_upperleft = camera_center - Vector3D(0, 0, focal_length) - vu / 2 - vv / 2;
    auto zerozero_pixel = v_upperleft + 0.5 * (diff_u + diff_v);


    

    Vector3D centerOfRotation(-1, 0, 2);
    double radius = (camera_center - centerOfRotation).length();

    // Image Plane Setup -perspective camera-
    Vector3D image_plane_center = camera_center + distToImage * unit_vector(Vector3D(0, 0, -1)); 
    Vector3D image_plane_upperleft = image_plane_center - (vu / 2.0) - (vv / 2.0);



    unsigned char image[width * height * 3];
    int depth = 5;
    for (int i = 0; i < height; i++)
    {
        std::clog << "\r scanlines remaining: " << (height - i) << std::flush;

        for (int j = 0; j < width; j++)
        {

            auto pix_center = zerozero_pixel + (i * diff_u) + (j * diff_v);
            //PerspectiveCamera condition
            Vector3D ray_dir;
            Vector3D w = unit_vector(camera_center - pix_center); // normalized direction towards the scene
            Vector3D u = diff_u * i;
            Vector3D v = diff_v * j;

            if (usePerspectiveGeometry)  // If perspective view is chosen
            {
                //ray_dir = -distToImage * w + u + v;
                //ray_dir = unit_vector(ray_dir); // Normalize the ray direction
                Vector3D position_on_image_plane = image_plane_upperleft + (i * diff_u) + (j * diff_v);
                ray_dir = position_on_image_plane - camera_center;

            }

            //ray r;
            //old camera settings for ray dir
            else
            {
                ray_dir = pix_center - camera_center;
            }

            //auto ray_dir = pix_center - camera_center;
            ray_dir = unit_vector(ray_dir);  // normalize//perspective
            ray r(camera_center, ray_dir);




            color pix_color = ray_color(r, groundPlane, depth);

            int idx = (i * width + j) * 3;
            image[idx] = write_colorR(std::cout, pix_color);
            image[idx + 1] = write_colorG(std::cout, pix_color);
            image[idx + 2] = write_colorB(std::cout, pix_color);
        }
    }
    
    
    
    //Generate a movie variables
    bool userMovieInput = false;

    
    
    
    
    cout << "Do you want to generate a movie ? (It is recommended you generate a movie only after you reviewed orthogonal and the perspective view of the scene as movie generation is time consuming) (Y/N)" << endl;
    char movieinput;
    cin >> movieinput;

    if (movieinput == 'Y' || movieinput == 'y') {
        userMovieInput = true;
    }



  
    if(userMovieInput == true)
    {

        for (int frame = 0; frame < numFrames; frame++) {

            // Calculate the current angle
            double currentAngle = (2.0 * pi / numFrames) * frame;

            // Move the camera
            camera_center.setX(centerOfRotation.x() + radius * cos(currentAngle));
            camera_center.setZ(centerOfRotation.z() + radius * sin(currentAngle));

            // Re-render the scene with the new camera position
            v_upperleft = camera_center - Vector3D(0, 0, focal_length) - vu / 2 - vv / 2;
            zerozero_pixel = v_upperleft + 0.5 * (diff_u + diff_v);


            int depth = 5;
            for (int i = 0; i < height; i++)
            {
                std::clog << "\r scanlines remaining: " << (height - i) << std::flush;

                for (int j = 0; j < width; j++)
                {

                    auto pix_center = zerozero_pixel + (i * diff_u) + (j * diff_v);

                    auto ray_dir = pix_center - camera_center;
                    //std::cout << pix_center << ray_dir;
                    ray r(camera_center, ray_dir);

                    color pix_color = ray_color(r, groundPlane, depth);

                    int idx = ((height - 1 - i) * width + j) * 3;//(i * width + j) * 3;
                    image[idx] = write_colorR(std::cout, pix_color);
                    image[idx + 1] = write_colorG(std::cout, pix_color);
                    image[idx + 2] = write_colorB(std::cout, pix_color);
                }
            }
            //char filename[256];
            std::ostringstream oss;
            oss << "output_" << std::setw(3) << std::setfill('0') << frame << ".png";
            std::string filename = oss.str();

            // Save the image to a file
            saveImage(filename.c_str(), image, width, height);
        }

    }
    


    // std::clog << "\r Done.    \n";

    unsigned char* data = &image[0];

    if (data)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    else
    {
        std::cout << "Failed to load texture" << std::endl;
    }





    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // bind Texture
        glBindTexture(GL_TEXTURE_2D, texture);

        // render container
        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}