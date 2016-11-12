/*
  File Name    : ppm.c
  Assignment   : Project 4 - Raytracing
  Created by   : Harika Hari (hh453). */


#include "custom.h"

 int ImageWrite(Image *buffer, const char *filename,int format) {
    size_t num2;
    int size = buffer->width * buffer->height * 4;
	int i;
    FILE *fp = fopen(filename, "w");
    if (!fp) { 
        fprintf(stderr,"cannot open file for writing \n");
		exit(0);
	}
    
	//write the header file
    //image format
	if(format  && (format == 3 || format == 6))
    fprintf(fp, "P%d\n",format);
	else {
	fprintf(stderr,"invalid image format \n");
	}
    //comments
    fprintf(fp, "# Created by %s\n","Harika");

    //image size
	if(buffer->width && buffer->height)
    fprintf(fp, "%d %d\n",buffer->width, buffer->height);

	else {
	fprintf(stderr,"invalid height or width of image\n");
	}
    // rgb component depth
    fprintf(fp, "%d\n",buffer->maxval);

    // writing pixel data to file
	for(i=1; i<size+1;i++){    
            char ch=buffer->data[i-1];
            if (i%4 !=0) {
               fwrite(&ch, 1, 1, fp);  //writing the image pixel to file
            }
   	}
	fclose(fp);
	return 1;
}