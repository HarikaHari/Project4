all: gcc raycast.c -lm

run-sphere: ./a.out 1000 1000 sphere_input.json sphere_output.ppm

run-quadric: ./a.out 1000 1000 quad_input.json quad_output.ppm