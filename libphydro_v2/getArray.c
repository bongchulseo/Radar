#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "libphydro.h"

/** Dynamic allocation **/

/* 1-dimensional array */
char* get1dCharArray(int dim_x)
{
    int i;
    char *array;

    if ((array=(char*)malloc(dim_x*sizeof(char))) == NULL) {
         exit(-1);
    }

    for (i=0; i<dim_x; i++) {
         array[i]=INITIAL_VALUE;
    }

    return array;
}
short* get1dShortArray(int dim_x)
{
    int i;
    short *array;

    if ((array=(short*)malloc(dim_x*sizeof(short))) == NULL) {
         exit(-1);
    }

    for (i=0; i<dim_x; i++) {
         array[i]=INITIAL_VALUE;
    }

    return array;
}
int* get1dIntArray(int dim_x)
{
    int i;
    int *array;

    if ((array=(int*)malloc(dim_x*sizeof(int))) == NULL) {
         exit(-1);
    }

    for (i=0; i<dim_x; i++) {
         array[i]=INITIAL_VALUE;
    }

    return array;
}
float* get1dFloatArray(int dim_x)
{
    int i;
    float *array;

    if ((array=(float*)malloc(dim_x*sizeof(float))) == NULL) {
         exit(-1);
    }

    for (i=0; i<dim_x; i++) {
         array[i]=INITIAL_VALUE;
    }

    return array;
}
double* get1dDoubleArray(int dim_x)
{
    int i;
    double *array;

    if ((array=(double*)malloc(dim_x*sizeof(double))) == NULL) {
         exit(-1);
    }

    for (i=0; i<dim_x; i++) {
         array[i]=INITIAL_VALUE;
    }

    return array;
}

/* 2-dimensional array */
short** get2dShortArray(int dim_x,int dim_y)
{
    int i, j;
    short **array;

    if ((array=(short**)malloc(dim_x*sizeof(short*)))!=NULL) {
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(short*)malloc(dim_y*sizeof(short)))==NULL) {
                exit(-1);
            }
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
             array[i][j]=INITIAL_VALUE;
        }
    }

    return array;
}
int** get2dIntArray(int dim_x,int dim_y)
{
    int i, j;
    int **array;

    if ((array=(int**)malloc(dim_x*sizeof(int*)))!=NULL) {
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(int*)malloc(dim_y*sizeof(int)))==NULL) {
                exit(-1);
            }
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
             array[i][j]=INITIAL_VALUE;
        }
    }

    return array;
}
float** get2dFloatArray(int dim_x,int dim_y)
{
    int i, j;
	float **array;

    if ((array=(float**)malloc(dim_x*sizeof(float*)))!=NULL) {
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(float*)malloc(dim_y*sizeof(float)))==NULL) {
                exit(-1);
            }
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
             array[i][j]=INITIAL_VALUE;
        }
    }

    return array;
}
double** get2dDoubleArray(int dim_x,int dim_y)
{
    int i, j;
    double **array;

    if ((array=(double**)malloc(dim_x*sizeof(double*)))!=NULL) {
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(double*)malloc(dim_y*sizeof(double)))==NULL) {
                exit(-1);
            }
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
             array[i][j]=INITIAL_VALUE;
        }
    }

    return array;
}

/* 3-dimensinal array */
short*** get3dShortArray(int dim_x,int dim_y,int dim_z)
{
    int i, j, k;
    short ***array;

    if ((array=(short***)malloc(dim_x*sizeof(short**)))!=NULL){
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(short**)malloc(dim_y*sizeof(short*)))!=NULL){
                 for (j=0; j<dim_y; j++) {
					 if ((array[i][j]=(short*)malloc(dim_z*sizeof(short)))==NULL){
                          exit(-1);
                     }
                 }
            }
            else
                exit(-1);
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            for (k=0; k<dim_z; k++) {
                 array[i][j][k]=INITIAL_VALUE;
            }
        }
    }

    return array;
}
int*** get3dIntArray(int dim_x,int dim_y,int dim_z)
{
    int i, j, k;
    int ***array;

    if ((array=(int***)malloc(dim_x*sizeof(int**)))!=NULL){
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(int**)malloc(dim_y*sizeof(int*)))!=NULL){
                 for (j=0; j<dim_y; j++) {
					 if ((array[i][j]=(int*)malloc(dim_z*sizeof(int)))==NULL){
                          exit(-1);
                     }
                 }
            }
            else
                exit(-1);
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            for (k=0; k<dim_z; k++) {
                 array[i][j][k]=INITIAL_VALUE;
            }
        }
    }

    return array;
}
float*** get3dFloatArray(int dim_x,int dim_y,int dim_z)
{
    int i, j, k;
    float ***array;

    if ((array=(float***)malloc(dim_x*sizeof(float**)))!=NULL){
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(float**)malloc(dim_y*sizeof(float*)))!=NULL){
                 for (j=0; j<dim_y; j++) {
					 if ((array[i][j]=(float*)malloc(dim_z*sizeof(float)))==NULL){
                          exit(-1);
                     }
                 }
            }
            else
                exit(-1);
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            for (k=0; k<dim_z; k++) {
                 array[i][j][k]=INITIAL_VALUE;
            }
        }
    }

    return array;
}
double*** get3dDoubleArray(int dim_x,int dim_y,int dim_z)
{
    int i, j, k;
    double ***array;

    if ((array=(double***)malloc(dim_x*sizeof(double**)))!=NULL){
        for (i=0; i<dim_x; i++) {
            if ((array[i]=(double**)malloc(dim_y*sizeof(double*)))!=NULL){
                 for (j=0; j<dim_y; j++) {
					 if ((array[i][j]=(double*)malloc(dim_z*sizeof(double)))==NULL){
                          exit(-1);
                     }
                 }
            }
            else
                exit(-1);
        }
    }
    else
        exit(-1);

    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            for (k=0; k<dim_z; k++) {
                 array[i][j][k]=INITIAL_VALUE;
            }
        }
    }

    return array;
}

/** Release allocated memory **/

/* 1-dimensional array */
int free1dArray(void *array)
{
    int status;

    status=FALSE;
    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

/* 2-dimensional array */
int free2dArray(void **array,int dim_x)
{
    int status, i;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        if (array[i] != NULL) {
            free(array[i]);
        }
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free2dShortArray(short **array,int dim_x)
{
    int status, i;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        if (array[i] != NULL) {
            free(array[i]);
        }
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free2dIntArray(int **array,int dim_x)
{
    int status, i;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        if (array[i] != NULL) {
            free(array[i]);
        }
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free2dFloatArray(float **array,int dim_x)
{
    int status, i;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        if (array[i] != NULL) {
            free(array[i]);
        }
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free2dDoubleArray(double **array,int dim_x)
{
    int status, i;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        if (array[i] != NULL) {
            free(array[i]);
        }
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

/* 3-dimensional array */
int free3dArray(void ***array,int dim_x,int dim_y)
{
    int status, i, j;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            if (array[i][j] != NULL)
                free(array[i][j]);
        }
        if (array[i] != NULL)
            free(array[i]);
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free3dIntArray(int ***array,int dim_x,int dim_y)
{
    int status, i, j;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            if (array[i][j] != NULL)
                free(array[i][j]);
        }
        if (array[i] != NULL)
            free(array[i]);
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free3dFloatArray(float ***array,int dim_x,int dim_y)
{
    int status, i, j;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            if (array[i][j] != NULL)
                free(array[i][j]);
        }
        if (array[i] != NULL)
            free(array[i]);
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

int free3dDoubleArray(double ***array,int dim_x,int dim_y)
{
    int status, i, j;

    status=FALSE;
    for (i=0; i<dim_x; i++) {
        for (j=0; j<dim_y; j++) {
            if (array[i][j] != NULL)
                free(array[i][j]);
        }
        if (array[i] != NULL)
            free(array[i]);
    }

    if (array != NULL) {
        free(array);
        status=TRUE;
    }

    return status;
}

void initialize3dShortArray(short ***array,
                            int dim_x,
                            int dim_y,
                            int dim_z,
                            short data_value)
{
    int i, j, k;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              for (k=0; k<dim_z; k++) {
                   array[i][j][k]=data_value;
              }
         }
    }
}
void initialize2dShortArray(short **array,
                            int dim_x,
                            int dim_y,
                            short data_value)
{
    int i, j;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              array[i][j]=data_value;
         }
    }
}
void initialize1dShortArray(short *array,
                            int dim_x,
                            short data_value)
{
    int i;

    for (i=0; i<dim_x; i++) {
         array[i]=data_value;
    }
}
void initialize3dIntArray(int ***array,
                          int dim_x,
                          int dim_y,
                          int dim_z,
                          int data_value)
{
    int i, j, k;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              for (k=0; k<dim_z; k++) {
                   array[i][j][k]=data_value;
              }
         }
    }
}
void initialize2dIntArray(int **array,
                          int dim_x,
                          int dim_y,
                          int data_value)
{
    int i, j;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              array[i][j]=data_value;
         }
    }
}
void initialize1dIntArray(int *array,
                          int dim_x,
                          int data_value)
{
    int i;

    for (i=0; i<dim_x; i++) {
         array[i]=data_value;
    }
}
void initialize3dFloatArray(float ***array,
                            int dim_x,
                            int dim_y,
                            int dim_z,
                            float data_value)
{
    int i, j, k;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              for (k=0; k<dim_z; k++) {
                   array[i][j][k]=data_value;
              }
         }
    }
}
void initialize2dFloatArray(float **array,
                            int dim_x,
                            int dim_y,
                            float data_value)
{
    int i, j;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              array[i][j]=data_value;
         }
    }
}
void initialize1dFloatArray(float *array,
                            int dim_x,
                            float data_value)
{
    int i;

    for (i=0; i<dim_x; i++) {
         array[i]=data_value;
    }
}
void initialize3dDoubleArray(double ***array,
                             int dim_x,
                             int dim_y,
                             int dim_z,
                             double data_value)
{
    int i, j, k;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              for (k=0; k<dim_z; k++) {
                   array[i][j][k]=data_value;
              }
         }
    }
}
void initialize2dDoubleArray(double **array,
                             int dim_x,
                             int dim_y,
                             double data_value)
{
    int i, j;

    for (i=0; i<dim_x; i++) {
         for (j=0; j<dim_y; j++) {
              array[i][j]=data_value;
         }
    }
}
void initialize1dDoubleArray(double *array,
                             int dim_x,
                             double data_value)
{
    int i;

    for (i=0; i<dim_x; i++) {
         array[i]=data_value;
    }
}
