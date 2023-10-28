
# CAP 5705 Assignment 1: First-Hit Ray Tracer

A simple ray tracer to render a scene composed of simple objects (Two spheres of different colors and shading along with a patterned glazed surface plane)


## Authors

- Nikhil Sai Siddam Shetty
- Ganti Balasai Srikanth


## FAQ

#### How to run this program?

Extract the file which contains OpenGLTest.sln in your local Visual Studio Community application.

Run the project.

The program will ask the user to provide 2 inputs:

    1: To view the scene in either orthogonal or perspective camera view:    
        press 'p' to proceed to perspective view or select any other key to 
        continue with orthogonal view.
    2: To generate a movie(Y/N): 
        It is recommended to generate a movie only after viewing the 
        orthoganal and perspective camera views of the scene as it could take 
        some time to generate all 300 Frames of the movie. 

To genearate a movie go to the directory where the images are stored that is the project dir and run the below command

        ffmpeg -framerate 30 -i output_%03d.png -c:v libx264 -pix_fmt yuv420p output.mp4




## Demo

![Scene1](https://imgur.com/skxgnuF)


## Features

- Spheres with ambient, diffusion and specular lighting
- Bonus feature: A patterened plane to add some playfulness to the scene
- Glazed Plane reflecting the objects above it
- Orthogonal and perspective camera views based on user input
- Movie Generated based on user input


## Acknowledgements

 - [OpenGL Lighting tutorial](https://learnopengl.com/Lighting/Basic-Lighting)
 - [Raytracing Basics and code inspiration](https://raytracing.github.io/)
 - [How to write a Good readme](https://bulldogjob.com/news/449-how-to-write-a-good-readme-for-your-github-project)

