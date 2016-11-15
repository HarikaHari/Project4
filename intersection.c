/*
  File Name    : intersection.c
  Assignment   : Project 4 - Raytracing
  Created by   : Harika Hari (hh453). */

#include "custom.h"
#include "vector.c"
#include "ppm.c"

Vector diffuseColor = {0,0,0},specularColor = {0,0,0}, normal = {0,0,0};
Vector pixelColor = {0,0,0};
int hitIndex; double hitDistance;


//Returns the position of camera object from the JSON file
int getCameraPosition(OBJECT *objects){
    int i=0;
    while(i < ObjectsCount && objects[i].type !=0){
        if (objects[i].type == CAM) {
            return i;
        }
        i++;
    }
	 if (i == ObjectsCount) {
        fprintf(stderr, "Error: No Camera object is found to raycast \n");
        exit(1);
    }
    return -1;
}


// fill in pixel color to our image - need to flip the y axis due to the placement of the viewplane

void colorPixel(double *color, int row, int col,Image *image){
 //store results in image->data which store pixels
    image->data[row * image->width*4 + col*4] = (char) (MAX_COLOR_VAL * color[0]);
    image->data[row * image->width*4 + col*4+1] = (char)(MAX_COLOR_VAL * color[1]);
    image->data[row * image->width*4 + col*4+2]= (char)(MAX_COLOR_VAL * color[2]);
    image->data[row * image->width*4 + col*4+3]= 255;

}


//Ro means Ray origin - Rd means Ray direction (other end of the ray)
 
double planeIntersection(double *Ro, double *Rd, double *Pos, double *Norm){

    double a,d;
    normalize(Norm);
    a = VectorDotProduct(Norm, Rd);
    if (fabs(a) <0.0001) { 
	// checks with absolute value and verifies if the plane is parallel to ray 
        return -1;
    }
    Vector incident;
    VectorSubstraction(Pos, Ro, incident);
    d = VectorDotProduct(incident, Norm);
    double t = d/a; 
    if (t<0.0) {
	 // no intersection of ray with plane
        return -1;
    }
    return t; 
	}

double sphereIntersection(double *Ro, double *Rd, double *pos, double r) {

    double a, b, c;
    double t0,t1; 
	//calculating the coefficients of the sphere equation. 
    a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
    b = 2 * (Rd[0]*(Ro[0]-pos[0]) + Rd[1]*(Ro[1]-pos[1]) + Rd[2]*(Ro[2]-pos[2]));
    c = sqr(Ro[0]-pos[0]) + sqr(Ro[1]-pos[1]) + sqr(Ro[2]-pos[2]) - sqr(r);
	//calculate the determinant value to get the point of intersection
    double d = sqr(b) - 4*a*c;
    
    if (d < 0) 
	  return -1; //no intersection
    
    else if (d == 0) {
        t0 = -1*(b / (2*a));
        t1 = -1*(b / (2*a));
       }
    else {  
        t0 = (-1*b - sqrt(d))/2;
        t1 = (-1*b + sqrt(d))/2;
    }        
    if (t0 < 0 && t1 < 0) 
        return -1;

    else if (t0 < 0 && t1 > 0) 
        return t1;
    
    else if (t0 > 0 && t1 < 0) 
        return t0;
    
    else { 
        if (t0 <= t1)
            return t0;
        else
            return t1;
    }
    
}
//function to calculate quadric types
double quadricIntersection (double *Ro, double *Rd, double *pos, double *coefficient) {

    double a, b, c;
    double t; 
		//calculating the coefficients of the equation. 
a = (coefficient[0]) * sqr(Rd[0]) + (coefficient[1]) * sqr(Rd[1]) + (coefficient[2]) * sqr(Rd[2]) + (coefficient[3]) * (Rd[0]) * (Rd[1]) + (coefficient[4]) * (Rd[0]) * (Rd[2]) + (coefficient[5]) * (Rd[1]) * (Rd[2]);
    
    
    b = 2*(coefficient[0]) * (Ro[0] - pos[0]) * (Rd[0]) + 2*(coefficient[1]) * (Ro[1] - pos[1]) * (Rd[1]) + 2*(coefficient[2]) * (Ro[2] - pos[2]) * (Rd[2]) + (coefficient[3]) * ((Ro[0] - pos[0]) * (Rd[1]) + (Ro[1] - pos[1]) * (Rd[0])) + (coefficient[4]) * (Ro[0] - pos[0]) * (Rd[2]) + (coefficient[5]) * ((Ro[1] - pos[1]) * (Rd[2]) + (Rd[1]) * (Ro[2] - pos[2])) + (coefficient[6]) * (Rd[0]) + (coefficient[7]) * (Rd[1]) + (coefficient[8]) * (Rd[2]);
    
   // printf("\n a : %lf ::: b : %lf ::: c : %lf ", a,b,c);
    
    c = (coefficient[0]) * sqr(Ro[0] - pos[0]) + (coefficient[1]) * sqr(Ro[1] - pos[1]) + (coefficient[2]) * sqr(Ro[2] - pos[2]) + (coefficient[3]) * (Ro[0] - pos[0]) * (Ro[1] - pos[1]) + (coefficient[4]) * (Ro[0] - pos[0]) * (Ro[2] - pos[2]) + (coefficient[5]) * (Ro[1] - pos[1]) * (Ro[2] - pos[2]) + (coefficient[6]) * (Ro[0] - pos[0]) + (coefficient[7]) * (Ro[1] - pos[1]) + (coefficient[8]) * (Ro[2] - pos[2]) + (coefficient[9]);
		
		
    if(a == 0)
	 return (-c/b); 
	//calculate the determinant value to get the point of intersection
    double d = (sqr(b) - 4*a*c);
    
    if (d < 0) 
	  return -1;  //no intersection since the det is imaginary 
    
    else if (d == 0) 
        t = -1*(b / (2*a));
       
       
    else {  
    d = sqrt(d);
    // Tests which is the closest one and in front of camera.
    double t = (-b - d)/(2 * a);
    if (t > 0) return t;
    t = (-b + d)/(2 * a);
    if (t > 0) return t;
    
    return -1;
	}
    }	
	



double getRadialAttenuation(LIGHT *light, double distance) {
	//check if the light distance is infinity
	if(distance == INFINITY) {
		return 1.0;
	}
	//check if denominator is zero to remove divide by zero exception
    if (light->radial_a0 == 0 && light->radial_a1 == 0 && light->radial_a2 == 0) {
        return 1.0;
    }
    return 1.0 / (light->radial_a2 * sqr(distance) + light->radial_a1 * distance + light->radial_a1);
}

double getAngularAttenuation(LIGHT *light, Vector objDir) {

	//check if the light type is spot lights
	if (!(light->theta || light->theta == 0)) {
	if (light->direction == NULL) {
        fprintf(stderr, "Error: direction vector cant be NULL for spotlight\n");
        exit(1);
    }
    //use pi/180 to calculate the cosine angle value
    double cosineOfTheta = cos(light->theta * (M_PI / 180));
	//Vobject and Vlight dot product
    double VobjdotVlight = VectorDotProduct(objDir,light->direction);
    if (VobjdotVlight < cosineOfTheta){
        return 0.0;
    }
    return pow(VobjdotVlight, light->angular_a0);
	}
	else return 1.0; //return for point light
}

void calculateDiffuseColor(double *N, double *L, double *IL, double *KD) {
    
    double ndotl = VectorDotProduct(N, L);
    if (ndotl > 0) {
        double temp[3];
        temp[0] = KD[0] * IL[0];
        temp[1] = KD[1] * IL[1];
        temp[2] = KD[2] * IL[2];
        VectorScale(temp, ndotl, diffuseColor);
    }
}


void calculateSpecularColor(double *R, double *V, double *KS, double *IL, double ns) {
    
   double vdotr = VectorDotProduct(V, R);
    if (vdotr > 0) {
        double vrpowns = pow(vdotr, ns);
        Vector temp;
        temp[0] = KS[0] * IL[0];
        temp[1] = KS[1] * IL[1];
        temp[2] = KS[2] * IL[2];
        VectorScale( temp, vrpowns,specularColor );
    }
}
void initializePixelColors() {
   specularColor[0]  = 0;
   diffuseColor[0] = 0;
   pixelColor[0] = 0;
   specularColor[1]  = 0;
   diffuseColor[1] = 0;
   pixelColor[1] = 0;
   specularColor[2]  = 0;
   diffuseColor[2] = 0;
   pixelColor[2] = 0;
}

void computeIntersection(Vector Ro, Vector Rd, int selfObjectIndex, double maxDist) {
   int objIndex = -1, counter = 0;
            double objDistance = INFINITY;
            for (counter=0; objects[counter].type!=0; counter++) {
			if (selfObjectIndex == counter) continue;
                double t =0;
                switch (objects[counter].type) {
                    case 0:
                        printf("no object found\n");
                        break;
                        
                    case CAM:
                        break;
                        
                    case SPH:
                        t = sphereIntersection(Ro, Rd, objects[counter].data.sphere.position,objects[counter].data.sphere.radius);
						break;
                        
                    case PLN:
                        t = planeIntersection(Ro, Rd, objects[counter].data.plane.position, objects[counter].data.plane.normal);
                        break;
						
					case QUAD:
                        t = quadricIntersection(Ro, Rd, objects[counter].data.quadric.position, objects[counter].data.quadric.coefficients);
                        break;
                        
                    default:
                        exit(1);
                }
                if (maxDist != INFINITY && t > maxDist)
                  continue;
				if (t > 0 && t < objDistance) {
		            //printf("\n objects[counter].type %d %lf",objects[counter].type,t);
                    objDistance = t;
                    objIndex = counter;
					//printf("selfObjectIndex %d : hitDistance %ld",hitIndex,hitDistance);
                }
                
            }
			
    hitIndex = objIndex;
    hitDistance = objDistance;
}

void computeNormal(int objIndex, Vector newRo, Vector Rd) {
   if(objects[objIndex].type != CAM) {
   if (objects[objIndex].type == PLN) {
      VectorCopy(objects[objIndex].data.plane.normal, normal);
    }
    else if (objects[objIndex].type == SPH) {
       VectorSubstraction(newRo, objects[objIndex].data.sphere.position, normal);
    }
    else if (objects[objIndex].type == QUAD) {
        Vector temp ={0,0,0};
		double *quad_coef = objects[objIndex].data.quadric.coefficients;
		double *pos = objects[objIndex].data.quadric.position;
        normal[0] =-2*quad_coef[0]*(pos[0])-quad_coef[3]*(pos[1])-quad_coef[4]*(pos[2])-quad_coef[6];
	    normal[1] =-2*quad_coef[1]*(pos[1])-quad_coef[3]*(pos[0])-quad_coef[5]*(pos[2])-quad_coef[7];
	    normal[2] =-2*quad_coef[2]*(pos[2])-quad_coef[4]*(pos[0])-quad_coef[5]*(pos[1])-quad_coef[8];
	    normalize(normal);
        }
		else {
        fprintf(stderr, "Error: Can't get normal vector for this Object Type\n");
        }
		}
}

//calculate the color of objects with illumination
void computeIlluminationColor(Vector light_Ro, Vector light_Rd, Vector Rd, int objIndex, double lightDistance,LIGHT *light, int ref) {
         
	double diffuseTemp[3] = {0,0,0}, specularTemp[3] = {0,0,0};	

			if(objects[objIndex].type != CAM) {
            if (objects[objIndex].type == PLN) {
			    VectorCopy(objects[objIndex].data.plane.normal, normal);
                VectorCopy(objects[objIndex].data.plane.diffuse_color, diffuseTemp);
				if(!(objects[objIndex].data.plane.specular_color)) 
                VectorCopy(objects[objIndex].data.plane.specular_color, specularTemp);	
            }
			else if (objects[objIndex].type == SPH) {
			    VectorSubstraction(light_Ro, objects[objIndex].data.sphere.position, normal);
                VectorCopy(objects[objIndex].data.sphere.diffuse_color, diffuseTemp);
                VectorCopy(objects[objIndex].data.sphere.specular_color, specularTemp);
            } 
			else if (objects[objIndex].type == QUAD) {
			 Vector temp ={0,0,0};
				double *quad_coef = objects[objIndex].data.quadric.coefficients;
				double *pos = objects[objIndex].data.quadric.position;
				normal[0] =-2*quad_coef[0]*(pos[0])-quad_coef[3]*(pos[1])-quad_coef[4]*(pos[2])-quad_coef[6];
				normal[1] =-2*quad_coef[1]*(pos[1])-quad_coef[3]*(pos[0])-quad_coef[5]*(pos[2])-quad_coef[7];
				normal[2] =-2*quad_coef[2]*(pos[2])-quad_coef[4]*(pos[0])-quad_coef[5]*(pos[1])-quad_coef[8];
			
            VectorCopy(objects[objIndex].data.quadric.diffuse_color,diffuseTemp);
            VectorCopy(objects[objIndex].data.quadric.specular_color, specularTemp);
            }
	        else {
                fprintf(stderr, "Error: not a valid object\n");
                exit(1);
            }
						
			//variables to start and read angular and radial attenuation 
			double fang = 1,frad = 1;
			Vector L,R,V;
			VectorCopy(light_Rd, L);
			normalize(normal);	
            normalize(L);
            VectorReflection(L, normal, R);
            VectorCopy(Rd, V);
			Vector objectToLightVector;
			//VectorSubstraction(intersection,lights[i].position,objectToLightVector);
			VectorCopy(L, objectToLightVector);
            VectorScale(objectToLightVector, -1, objectToLightVector);
			normalize(objectToLightVector);
			if(ref) {
			frad = getRadialAttenuation(light, lightDistance);
			fang = getAngularAttenuation(light, objectToLightVector);
			}			
            //using vectors specified by DR Palmer to calculate diffuse and specular colors
			
            calculateDiffuseColor(normal, L, light->color, diffuseTemp);
	         //To caluculate specular color set ns to 20
            calculateSpecularColor(R, V, specularTemp, light->color, 20);

            pixelColor[0] += frad*fang*(specularColor[0] + diffuseColor[0]);
            pixelColor[1] += frad*fang*(specularColor[1] + diffuseColor[1]);
            pixelColor[2] += frad*fang*(specularColor[2] + diffuseColor[2]);
	
    }
}

void computePixelColor(double *pixelColor) {
	int i;
	for(i = 0; i<3; i++) {
	
	if (pixelColor[i] < 0) 
		pixelColor[i] = 0;		
				
	else if (pixelColor[i] > 1) 
		pixelColor[i] = 1;
	}        
}

double getReflectionCoefficient(int objIndex){
    if(objects[objIndex].type != CAM) {
	if (objects[objIndex].type == PLN) 
        return objects[objIndex].data.plane.reflection;
    
    else if (objects[objIndex].type == SPH) 
        return objects[objIndex].data.sphere.reflection;
    
    else if (objects[objIndex].type == QUAD)
        return objects[objIndex].data.quadric.reflection;
    
    else {
        fprintf(stderr, "Error: Can't get reflectivity for this Object Type\n");
        return -1;
    }
	}
}

double getRefractiveCoeffiecient(int objIndex) {
    if(objects[objIndex].type != CAM) {
	if (objects[objIndex].type == PLN) 
        return objects[objIndex].data.plane.refraction;
    
    else if (objects[objIndex].type == SPH) 
        return objects[objIndex].data.sphere.refraction;
    
    else if (objects[objIndex].type == QUAD) 
        return objects[objIndex].data.quadric.refraction;
    
    else {
        fprintf(stderr, "Error: Can't get refractivity for this Object Type\n");
        return -1;
    }
	}
}

double getIor(int objectIndex){
    double ior = 0.0;
	if(objects[objectIndex].type != CAM) {
    if (objects[objectIndex].type == PLN) 
        ior = objects[objectIndex].data.plane.ior;
    
    else if (objects[objectIndex].type == SPH) 
        ior = objects[objectIndex].data.sphere.ior;
    
    else if(objects[objectIndex].type == QUAD)
        ior = objects[objectIndex].data.quadric.ior;
    
    else {
        fprintf(stderr, "Error: Can't get ior for this Object Type\n");
        exit(1);
    }
	}
    if (fabs(ior) < 0.0001)
        return 1;
    else
        return ior;
}

void computeRefractionValue(Vector Rd, Vector newRo, int objectIndex, double outIor, Vector *refractedVector) {
    // initializations and variables setup
    Vector dir, pos;
    VectorCopy(Rd, dir);
    VectorCopy(newRo, pos);
    normalize(dir);
    normalize(pos);
    
    double inIor = getIor(objectIndex);
    
    if (inIor == outIor) {
        inIor =1;
    }
    Vector a, b;
    
    // find normal vector of current object
	computeNormal(objectIndex, newRo, Rd);
    normalize(normal);
    
    // create coordinate frame with a and b, where b is tangent to the object intersection
    vectorCorssProduct(normal, dir, a);
    normalize(a);
    vectorCorssProduct(a, normal, b);
    
    // find transmission vector angle and direction
    double sin_theta = VectorDotProduct(dir, b);
    double sin_phi = (outIor / inIor) * sin_theta;
    double cos_phi = sqrt(1 - sqr(sin_phi));
    
    VectorScale(normal, -1*cos_phi, normal);
    VectorScale(b, sin_phi, b);
    VectorAddition(normal , b, *refractedVector);

}

void computeIlluminationAndReflectionColor(Vector Ro, Vector Rd, int objIndex, double objDistance, int level, double ior) {

    if (level > MAXREC) {
		VectorScale(pixelColor, 0, pixelColor);
        return;
    }
	
	if (objIndex == -1) {  
        return;
    }
    Vector newRo = {0,0,0};
    Vector newRd = {0,0,0};
	
	Vector reflectedRo;
	Vector reflectedRd;
	Vector refractedRo;
	Vector refractedRd;
		
    int reflection_index, refraction_index;
    double reflection_dist, refraction_dist;
    
	
    // finding new rays origin Ro and dir Rd vectors
    VectorScale(Rd, objDistance, newRo);
    VectorAddition(newRo, Ro, newRo);
    
    Vector reflection = {0,0,0};
	Vector refraction = {0,0,0};
    normalize(Rd);

    computeNormal(objIndex, newRo, Rd);
    normalize(normal);
    VectorReflection(Rd, normal, reflection);
    
	computeRefractionValue(Rd, newRo, objIndex, ior, &refraction);
	
	VectorCopy(newRo, reflectedRo);
	VectorCopy(reflection, reflectedRd);
	
	VectorCopy(newRo, refractedRo);
	VectorCopy(refraction, refractedRd);
	
	normalize(reflectedRd);
	normalize(refractedRd);
    
	computeIntersection(reflectedRo, reflectedRd ,objIndex, INFINITY);
	reflection_index =  hitIndex ;
    reflection_dist = hitDistance;
		
	if (objects[objIndex].type == PLN)
        computeIntersection(refractedRo, refractedRo ,objIndex, INFINITY);
    else
        computeIntersection(refractedRo, refractedRd ,-1, INFINITY);
		
	refraction_index =  hitIndex ;
    refraction_dist = hitDistance;
    
    if (reflection_index == -1 && refraction_index == -1) 
        VectorScale(pixelColor, 0, pixelColor);
        
    else {
        // recursive ray tracing 
        Vector reflectionColor ={0,0,0};
		Vector refractionColor ={0,0,0};
		
		double reflectIor =1;
        double refractIor =1;
		
        double reflectivity = getReflectionCoefficient(objIndex);
		double refractivity = getRefractiveCoeffiecient(objIndex);
        
        LIGHT reflectionLight;
        reflectionLight.direction = malloc(sizeof(Vector));
        reflectionLight.color = malloc(sizeof(Vector));
		
	if (reflection_index >= 0) {
            reflectIor = getIor(reflection_index);

            computeIlluminationAndReflectionColor(reflectedRo, reflectedRd, reflection_index, reflection_dist, level+1,reflectIor);
            VectorScale(reflectionColor, reflectivity, reflectionColor);
            VectorScale(reflection, -1, reflectionLight.direction);
            VectorCopy(reflectionColor, reflectionLight.color);
            
            VectorScale(reflectedRd, reflection_dist, reflectedRd);
            VectorSubstraction(reflectedRd, newRo, newRd);
            normalize(newRd);
            computeIlluminationColor(newRo, newRd, Rd, objIndex, INFINITY, &reflectionLight, 0);
        }
		
	LIGHT refractionLight;
        refractionLight.direction = malloc(sizeof(Vector));
        refractionLight.color = malloc(sizeof(Vector));
		
        if (refraction_index >= 0) {
            refractIor = getIor(refraction_index);
            VectorScale(refractedRd, 0.01, refractedRd);
            computeIlluminationAndReflectionColor(refractedRo, refractedRd, refraction_index, refraction_dist, level+1,refractIor);
            VectorScale(refractionColor, refractivity, refractionColor);
            VectorScale(refraction, -1, refractionLight.direction);
            VectorCopy(refractionColor, refractionLight.color);
            
            VectorScale(refractedRd, refraction_dist,refractedRd);
            VectorSubstraction(refractedRd, newRo, newRd);
            normalize(newRd);
            computeIlluminationColor(newRo, newRd, Rd, objIndex, INFINITY, &refractionLight, 0);
        }
        if (reflectivity == -1) {
            reflectivity = 0;
        }
        if (refractivity == -1) {
            refractivity = 0;
        }
		
         if (fabs(reflectivity) > 0.00001 && fabs(refractivity) > 0.00001) {
            
            double color_coefficient = 1.0 - reflectivity - refractivity;
            if (fabs(color_coefficient) < 0.0001) 
                color_coefficient = 0;
            Vector object_color = {0, 0, 0};
            VectorCopy(objects[objIndex].data.plane.diffuse_color, object_color);
            VectorScale(object_color, color_coefficient, object_color);
            pixelColor[0] += object_color[0];
            pixelColor[1] += object_color[1];
            pixelColor[2] += object_color[2];
        }
        //removing the allocated memory for light sources through malloc
        free(reflectionLight.direction);
        free(reflectionLight.color);
        free(refractionLight.direction);
        free(refractionLight.color);
    }
	int i;
    for (i=0; lights[i].color != NULL; i++) {
        // find new ray direction
	
		int lightIndex;
        double light_dist;
		VectorScale(newRd , 0 , newRd);
        VectorSubstraction(lights[i].position, newRo, newRd);
        double lightDistance = vectorLength(newRd);
        normalize(newRd);
        
        // new check new ray for intersections with other objects
        computeIntersection(newRo, newRd, objIndex, lightDistance);
        light_dist =  hitDistance;
		lightIndex =  hitIndex;
			
	if (lightIndex == -1) // this means there was no object in the way between the current one and the light
            computeIlluminationColor(newRo, newRd, Rd, objIndex, lightDistance, &lights[i], 1);
        
    }
}

void raycast(Image *image, double cameraWidth, double cameraHeight) {
   
    int x, y, counter; 

    Vector viewPlanePosition= {0,0,1}; //view plane position
    Vector Ro = {0,0,0}; //Camera position
	Vector Rd = {0,0,0}; // ray direction
    Vector point = {0,0,0}; //initial point on view plane
    
    double pixheight = cameraHeight / image->height;
    double pixwidth = cameraWidth / image->width;
    
   
     // set viewplane to Z direction
    
    for(x=0;x<image->height;x++){
        
        for (y=0; y<image->width; y++) {
		     
             point[0] = viewPlanePosition[0] - cameraWidth/2.0 + pixwidth*(y + 0.5);
			 point[1] = -(viewPlanePosition[1] - cameraHeight/2.0 + pixheight*(x + 0.5)); 
			 point[2] = viewPlanePosition[2];
            normalize(point); 		
            Rd[0] = point[0];
            Rd[1] = point[1];
            Rd[2] = point[2];
            
              // use this to compute the final color of the pixel	
	       initializePixelColors();

			computeIntersection(Ro, Rd, -1, INFINITY);
			double objDistance =  hitDistance;
			int objIndex =  hitIndex;
							 
            if (objIndex != -1 && objDistance > 0 && objDistance != INFINITY) {		
			   int level = 0;
			   double ior = 0.0;
			   computeIlluminationAndReflectionColor(Ro, Rd, objIndex, objDistance, level, ior);
			   computePixelColor(pixelColor);
			   colorPixel(pixelColor, x, y, image); 
            }
            else{
			//colouring the pixel to default color since there was no intersection
			Vector defaultColor ={0,0,0};
			colorPixel(defaultColor,x,y,image);
            }
            
        }
    }
    
}

