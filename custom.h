/*
  File Name    : custom.h
  Assignment   : Project 4 - Raytracing
  Created by   : Harika Hari (hh453). */

#ifndef custom_headers
#define custom_headers

//All Library Includes that are required
#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Pre Defined Values 
#define ObjectsCount 1024  //maximum number of JSON objects that can be parsed
#define MAX_COLOR_VAL 255 //max value of RGB color 
#define CAM 1
#define SPH 2
#define PLN 3
#define QUAD 4
#define LITE 5
#define REFLECTION 7
#define MAXREC 7

//structures based on the JSON sample
typedef struct CAMERA {
    double width;
    double height;
} CAMERA;

typedef struct SPHERE {
    double *diffuse_color;
	double *specular_color;
    double *position;
    double radius;
	double reflection;
    double refraction;
    double ior;
} SPHERE;

typedef struct PLANE {
    double *diffuse_color;
	double *specular_color;
    double *position;
    double *normal;
    double reflection;
    double refraction;
    double ior;
} PLANE;

typedef struct QUADRIC {
    double *diffuse_color;
	double *specular_color;
    double *position;
	double *coefficients;
	double reflection;
    double refraction;
    double ior;
} QUADRIC;

//Structure for light object
typedef struct LIGHT{
    double *color;
    double *position;
    double *direction;
    double radial_a0;
    double radial_a1;
    double radial_a2;
	double theta;
    double angular_a0;
}LIGHT;

//structure to store data as an object after parsing JSON data
typedef struct OBJECT{
    int type;
    union {
        CAMERA camera;
        SPHERE sphere;
        PLANE plane;
        QUADRIC quadric;
    } data;
} OBJECT;

//structure of image to store image data
typedef struct Image
{
    int width;
    int height;
    char *data; 
    int depth;
    char *tupltype;
    int maxval;
}Image;

// structure for ray equation
typedef struct ray{
    double origin[3];
    double direction[3];
} ray;

//custom defined type Vector
typedef double Vector[3]; 

//member functions to read JSON data and write image
void read_scene(const char* filename);
int ImageWrite(Image *image, const char *filename,int format);

//member functions to raycast and color image intersections 
int getCameraPosition(OBJECT *objects);
void colorPixel(double *color, int row, int col,Image *image);
double planeIntersection(double *Ro, double *Rd, double *Pos, double *Norm);
double sphereIntersection(double *Ro, double *Rd, double *C, double r);
double quadricIntersection (double *Ro, double *Rd, double *pos, double *coefficient);
void raycast(Image *image, double cam_width, double cam_height, OBJECT *objects, LIGHT *lights);

//member functions for vector calculations
double sqr(double v);
double vectorLength(Vector a);
double vectorDistance(Vector a, double *b);
void VectorAddition(Vector a, Vector b, Vector c);
void VectorSubstraction(Vector a, Vector b, Vector c);
void VectorCopy(Vector a, Vector b);
double VectorDotProduct(Vector a, Vector b);
void VectorReflection(Vector a, Vector b, Vector c);
void normalize(double *v);

#endif /* custom_h */
