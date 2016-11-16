# Project 4 - Recursive Raytracing

This project build on C language on unix environment's gcc compiler is intended to implement a recursive ray tracing algorithm to provide reflection and refraction and color objects based on the reflection of light. The reflective and refractive indexes are calculated for mathematical primitives based on an input file into a pixel buffer in positive z- axis direction. Also the program writes the image output in P6 format to a ppm file. In this process of conversion the program will be able to handle undefined errors with a prefixed error message.

There are two input json files: one is for sphere and plane types(Example specified by Dr. Palmer), the other uses the quadric values to run the test for quadric type.

Steps to run the program :

To run this program - use gcc compiler on unix machine and follow the instructions specified in MakeFile.txt.
Step1: Execute all command from MakeFile to compile the program.
Step2: Use run-sphere command from MakeFile to calculate the reflection of light for sphere and plane types.
Step3: Use run-quadric command from MakeFile to calculate the reflection of light for quadric types.

Usage Steps with Example : 
./a.out <width> <height> <input.json> <output.ppm>

width: The width of the image
height: The height of the image
input.json :  The input scene file consisting all the mathematical primitives and their attribute values in JSON format.
output.ppm : The location/name of the output file with the file type extension.

Example : ./a.out 1000 1000 sphere_input.json sphere_output.ppm