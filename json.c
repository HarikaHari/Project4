/*
  File Name    : json.c
  Assignment   : Project 4 - Raytracing
  Created by   : Harika Hari (hh453). */

#include "custom.h"

int line = 1;  // global variable, it will tells us which line is not correct

OBJECT objects[ObjectsCount]; // Array of All geometric Objects from JSON File
LIGHT lights[ObjectsCount]; //Array of all light sources from JSON file

//JSON PARSER CODE BY DR.PALMER START

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
    int c = fgetc(json);
#ifdef DEBUG
    printf("next_c: '%c'\n", c);
#endif
    if (c == '\n') {
        line++;;
    }
    if (c == EOF) {
        fprintf(stderr, "Error: next_c: Unexpected EOF: %d\n", line);
        exit(1);
    }
    return c;
}


// expect_c() checks that the next character is d.
// It is not d, give us an error.
void expect_c(FILE* json, int d) {
    int c = next_c(json);
    if (c == d) return;
    fprintf(stderr, "Error: Expected '%c': %d\n", d, line);
    exit(1);
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
    int c = next_c(json);
    while (isspace(c)) {
        c = next_c(json);
    }
    ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
    char buffer[129];
    int c = next_c(json);
    if (c != '"') {
        fprintf(stderr, "Error: Expected string on line %d.\n", line);
        exit(1);
    }
    c = next_c(json);
    int i = 0;
    while (c != '"') {
        if (i >= 128) {
            fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
            exit(1);
        }
        if (c == '\\') {
            fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
            exit(1);
        }
        if (c < 32 || c > 126) {
            fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
            exit(1);
        }
        buffer[i] = c;
        i += 1;
        c = next_c(json);
    }
    buffer[i] = 0;
    return strdup(buffer);
}

double next_number(FILE* json) {
    double value;
    fscanf(json, "%lf", &value);
    if (value == NAN) {
        fprintf(stderr, "Error: Expected a number but found NAN: %d\n", line);
        exit(1);
    }
    return value;
}


// Parse rgb colors and verify if the colors are valid
double* parse_color(FILE* json) {
    double* v = malloc(sizeof(double)*3);
    skip_ws(json);
    expect_c(json, '[');
    skip_ws(json);
    v[0] =  next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[1] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[2] = next_number(json);
    skip_ws(json);
    expect_c(json, ']');
    // check that all values are valid
    if ((v[0] < 0.0 || v[0] > 1.0) ||
        (v[1] < 0.0 || v[1] > 1.0) ||
        (v[2] < 0.0 || v[2] > 1.0)) {
        fprintf(stderr, "Error: rgb value out of range: %d\n", line);
        exit(1);
    }
    return v;
}

//parse light color 
double* parse_light_color(FILE* json) {
    double* v = malloc(sizeof(double)*3);
    skip_ws(json);
    expect_c(json, '[');
    skip_ws(json);
    v[0] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[1] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[2] = next_number(json);
    skip_ws(json);
    expect_c(json, ']');
    // check that all values are valid
    if ((v[0] < 0.0) ||
        (v[1] < 0.0) ||
        (v[2] < 0.0)) {
        fprintf(stderr, "Error: rgb value out of range: %d\n", line);
        exit(1);
    }
    return v;
}



// This function at here is for quadric
double* nextCoefficient(FILE* json){
    
	int i;
	double* v = malloc(10*sizeof(double));
	
    expect_c(json, '[');
    skip_ws(json);
    v[0] = next_number(json);
    
	for(i=1;i<10;i++){
        skip_ws(json);
        expect_c(json, ',');
        skip_ws(json);
        v[i] = next_number(json);
    }
    skip_ws(json);
    expect_c(json, ']');
    return v;
}


double* next_vector(FILE* json) {
    double* v = malloc(3*sizeof(double));
    expect_c(json, '[');
    skip_ws(json);
    v[0] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[1] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[2] = next_number(json);
    skip_ws(json);
    expect_c(json, ']');
    return v;
}

void read_scene(const char* filename) {
    FILE* json = fopen(filename, "r");
    int counter = 0;
	int lightCounter = 0;
    if (json == NULL) {
        fprintf(stderr, "Error: Could not open file\n");
        exit(1);
    }
    skip_ws(json);
    // Find the beginning of the list
    int c  = next_c(json);
	if (c != '[') {
        fprintf(stderr, "Error: JSON file must begin with [\n");
        exit(1);
    }
    skip_ws(json);
	c = next_c(json);
    //Conditional check to verify JSON file is not empty
    if (c == ']' || c == EOF) {
        fprintf(stderr, "Error: Empty JSON file\n");
        exit(1);
    }
	skip_ws(json);
    
    //Find the objects
    while (1) {
        if (counter > ObjectsCount) {
            fprintf(stderr, "Error: Number of objects is too large: %d\n", line);
            exit(1);
        }
        if (c == ']') {
            fprintf(stderr, "Error:  Unexpected ']': %d\n", line);
            fclose(json);
            exit(1);
        }
        if (c == '{') {     
            skip_ws(json);
            char *key = next_string(json);
            if (strcmp(key, "type") != 0) {
                fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
                exit(1);
            }
            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            char *type = next_string(json);
            int object_type;
            if (strcmp(type, "camera") == 0) {
                object_type = CAM;
                objects[counter].type = CAM;
			}
            else if (strcmp(type, "sphere") == 0) {
                object_type = SPH;
                objects[counter].type = SPH;
			}
            else if (strcmp(type, "plane") == 0) {
                object_type = PLN;
                objects[counter].type = PLN;
			}
            else if (strcmp(type, "quadric") == 0) {
                object_type = QUAD;
                objects[counter].type = QUAD;
			}
			else if (strcmp(type, "light") == 0)
                object_type = LITE;
            
            else {
                exit(1);
            }
            
            skip_ws(json);
            
            while (1) {
                //  , }
                c = next_c(json);
                if (c == '}') {
                    // stop parsing this object
                    break;
                }
                else if (c == ',') {
                    // read another field
                    skip_ws(json);
                    char* key = next_string(json);
                    skip_ws(json);
                    expect_c(json, ':');
                    skip_ws(json);
					double value;
                    if (strcmp(key, "width") == 0) {
						if(object_type == CAM) {
						value = next_number(json);
						if(value >= 0) {
                        objects[counter].data.camera.width = value;
						}
						else {
						fprintf(stderr, "Error: width cannot be negative %d\n", line);
                         exit(1);
						}
						}
						else {
						fprintf(stderr, "Error: width cannot be applied to other object types than CAMERA %d\n", line);
                         exit(1);
						}
                    }
                    else if (strcmp(key, "height") == 0) {
						value = next_number(json);
						if(object_type == CAM) {
						if(value >= 0) {
                        objects[counter].data.camera.height = value;
						}
						else {
						fprintf(stderr, "Error: height cannot be negative %d\n", line);
                         exit(1);
						}
						}
						else {
						fprintf(stderr, "Error: height cannot be applied to other object types than CAMERA %d\n", line);
                         exit(1);
						}
                    }
                    else if (strcmp(key, "radius") == 0) {
						if (object_type == SPH) {
                        value = next_number(json);
						if(value >= 0) {
                        objects[counter].data.sphere.radius = value;
						}
						else {
						fprintf(stderr, "Error: radius cannot be negative %d\n", line);
                         exit(1);
						}
						}
						 
						 else {
						 fprintf(stderr, "Error: radius can't be applied to other object types than sphere: %d\n", line);
                            exit(1);
						 }
                    }
					 else if (strcmp(key, "radial-a0")==0){
                        value = next_number(json);
						if (object_type == LITE) {
                        if(value >= 0) {
                            lights[lightCounter].radial_a0 = value;
                        }
                        else {
						fprintf(stderr, "Error: radial-a0 cannot be negative: %d\n", line);
                            exit(1);
						}
						}
						else {
						 fprintf(stderr, "Error: radial-a0 can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
						 }
                    }
                    else if (strcmp(key, "radial-a1")==0){
                        value = next_number(json);
						if (object_type == LITE) {
						if(value >= 0) {
                             lights[lightCounter].radial_a1 = value;
                        }
                        else {
						fprintf(stderr, "Error: radial-a1 cannot be negative: %d\n", line);
                            exit(1);
						}
						}
						else {
						 fprintf(stderr, "Error: radial-a1 can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
						 }
                    }
                    else if (strcmp(key, "radial-a2")==0){
                        value = next_number(json);
						if (object_type == LITE) {
						if(value >= 0) {
                             lights[lightCounter].radial_a2 = value;
                        }
                        else {
						fprintf(stderr, "Error: radial-a2 cannot be negative: %d\n", line);
                            exit(1);
						}
						}
						else {
						 fprintf(stderr, "Error: radial-a2 can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
						 }
						
                    }
                    else if (strcmp(key, "angular-a0")==0){
                        value = next_number(json);
						if (object_type == LITE) {
						if(value >= 0) {
                             lights[lightCounter].angular_a0 = value;
                        }
                        else {
						fprintf(stderr, "Error: angular_a0 cannot be negative: %d\n", line);
                            exit(1);
						}
						}
						else {
						 fprintf(stderr, "Error: angular-a0 can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
						 }
                    }
					
					else if (strcmp(key, "theta")==0){
                        value = next_number(json);
						if (object_type == LITE) {
						//if(value >= 0) {
                             lights[lightCounter].theta = value;
                       // }
                       // else {
						//fprintf(stderr, "Error: theta cannot be negative: %d\n", line);
                       //    exit(1);
						//}
						}
						else {
						 fprintf(stderr, "Error: theta can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
						 }
                    }
					
					 else if (strcmp(key, "color") == 0) {
                        if (object_type == LITE)
                            lights[lightCounter].color = parse_light_color(json);
                        else {
                            fprintf(stderr, "Error: color vector can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
                        }
                    }
					
					else if (strcmp(key, "direction") == 0) {
                        if (object_type == LITE)
                            lights[lightCounter].direction = next_vector(json);
                        else {
                            fprintf(stderr, "Error: direction vector can't be applied to other object types than LIGHT: %d\n", line);
                            exit(1);
                        }
                    }
					
                    else if (strcmp(key, "specular_color") == 0) {
                        if (object_type == SPH)
                            objects[counter].data.sphere.specular_color = parse_color(json);
                        else if (object_type == PLN)
                            objects[counter].data.plane.specular_color = parse_color(json);
                        else if (object_type == QUAD)
                            objects[counter].data.quadric.specular_color = parse_color(json);
                        else {
                            fprintf(stderr, "Error: specular_color vector can't be applied here: %d\n", line);
                            exit(1);
                        }
                    }
					else if (strcmp(key, "diffuse_color") == 0) {
                        if (object_type == SPH)
                            objects[counter].data.sphere.diffuse_color = parse_color(json);
                        else if (object_type == PLN)
                            objects[counter].data.plane.diffuse_color = parse_color(json);
                        else if (object_type == QUAD)
                            objects[counter].data.quadric.diffuse_color = parse_color(json);
                        else {
                            fprintf(stderr, "Error: diffuse_color vector can't be applied here: %d\n", line);
                            exit(1);
                        }
                    }
					
                    else if (strcmp(key, "position") == 0) {
                        if (object_type == SPH)
                            objects[counter].data.sphere.position = next_vector(json);
                        else if (object_type == PLN)
                            objects[counter].data.plane.position = next_vector(json);
                        else if (object_type == QUAD)
                            objects[counter].data.quadric.position = next_vector(json);
                        else if (object_type == LITE)
                            lights[lightCounter].position = next_vector(json);
                        
						else {
                            fprintf(stderr, "Error: Position vector can't be applied here: %d\n", line);
                            exit(1);
                        }
                        
                    }
					
					 else if (strcmp(key, "reflectivity") == 0) {
                        
						if(object_type == PLN)
                            objects[counter].data.plane.reflection = next_number(json);
                        
                        else if (object_type == SPH)
                            objects[counter].data.sphere.reflection = next_number(json);
                        
                        else if (object_type == QUAD)
                         objects[counter].data.quadric.reflection = next_number(json);
                        
                        else{
                            fprintf(stderr, "Error: reflectivity can't be applied here: %d\n", line);
                            exit(1);
                        }
                    }
                    else if (strcmp(key, "refractivity") == 0) {
                       if(object_type == PLN)
                            objects[counter].data.plane.refraction = next_number(json);
                        
                        else  if (object_type == SPH)
                             objects[counter].data.sphere.refraction = next_number(json);
                        
                        else if (object_type == QUAD)
                            objects[counter].data.quadric.refraction = next_number(json);
                        
                        else{
                            fprintf(stderr, "Error: Refractivity can't be applied here: %d\n", line);
                            exit(1);
                        }
                    }
                    else if (strcmp(key, "ior") == 0) {
                        if (object_type == PLN)
                            objects[counter].data.plane.ior = next_number(json);
                        
                        else if(object_type == SPH)
                            objects[counter].data.sphere.ior = next_number(json);
                        
                        else if (object_type == QUAD)
                            objects[counter].data.quadric.ior = next_number(json);
                        
                        else{
                            fprintf(stderr, "Error: ior can't be applied here: %d\n", line);
                            exit(1);
                        }
                    }
                    else if (strcmp(key, "normal") == 0) {
                        if (object_type != PLN) {
                            fprintf(stderr, "Error: Normal vector can't be applied to other object types than Plane: %d\n", line);
                            exit(1);
                        }
                        else
                            objects[counter].data.plane.normal = next_vector(json);
                    }
                    
                    else if (strcmp(key, "coefficient") == 0) {
                        if (object_type != QUAD) {
                            fprintf(stderr, "Error:  Normal vector can't be applied to other object types than Quadrics: %d\n", line);
                            exit(1);
                        }
                        else
                             objects[counter].data.quadric.coefficients = nextCoefficient(json);
                        
                    }
                    else {
                        fprintf(stderr, "Error: '%s' not a valid object: %d\n", key, line);
                        exit(1);
                    }
				skip_ws(json);
                }
                else {
                    fprintf(stderr, "Error: Unexpected value '%c': %d\n", c, line);
                    exit(1);
                }
            }
			if(type != NULL)
				printf("Parsing %s Completed...\n",type);
            skip_ws(json);
            c = next_c(json);
            if (c == ',') {
                // noop
                skip_ws(json);
            }
            else if (c == ']') {
                printf("Completed Parsing JSON file\n");
                fclose(json);
                return;
            }
            else {
                fprintf(stderr, "Error:Expecting comma or ]: %d\n", line);
                exit(1);
            }
        if (object_type == LITE) 
            lightCounter++;
        
        else
            counter++;
		}
		 
        
        c = next_c(json);
    }
	
//JSON PARSER CODE BY DR.PALMER END	
}
