/*
  File Name    : vector.c
  Assignment   : Project 4 - Raytracing
  Created by   : Harika Hari (hh453). */

#include "custom.h"

double sqr(double v){
    return v*v;
}

double vectorLength(Vector a){
    return sqrt(sqr(a[0])+sqr(a[1])+sqr(a[2]));
}

void vectorUnit(Vector a, Vector b){
    double len = vectorLength(a);
    b[0] = a[0]/len;
    b[1] = a[1]/len;
    b[2] = a[2]/len;
    
}

double vectorDistance(Vector a, double *b) {
	 double x =  sqr(a[0] - b[0]);
	 double y =	 sqr(a[1] - b[1]);
	 double z =  sqr(a[2] - b[2]);
	 return sqrt(x+y+z);
}

void VectorAddition(Vector a, Vector b, Vector c){
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

void VectorSubstraction(Vector a, Vector b, Vector c){
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

void VectorCopy(Vector a, Vector b){
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
}

double VectorDotProduct(Vector a, Vector b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


void VectorScale(Vector a, double s, Vector b){
    b[0] = s * a[0];
    b[1] = s * a[1];
    b[2] = s * a[2];
}

void VectorReflection(Vector a, Vector b, Vector c){
    normalize(b);
    double p = 2.0 * VectorDotProduct(a, b);
    Vector new;
    VectorScale(b, p, new);
    VectorSubstraction(a, new, c);
}

void normalize(double *v) {
    double len = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    len = sqrt(len);
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
}

void vectorCorssProduct(Vector a, Vector b, Vector c){
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

double calcFresenel(double a, double b, double mix)
{
    return b * mix + a * (1 - mix);
}

