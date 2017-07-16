#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "../headers/functions.h"
#include "functions.h"
#include <stdbool.h>

void vectNorm(double * v, int size, double * norm)
{
	double sum = 0.;
	double res;
	int i;
	
	for (i = 0; i < size; i++)
		sum += (v[i] * v[i]);
	
	res =  sqrt(sum);
	*norm = res;
	//printf("normeeee: %lf\n",*norm);
}

void symMatInit(double ** A, int n)
{
	int i,j;
	
	for(i=0;i<n;i++)
		for(j=i;j<n;j++)
		{
			A[i][j]=i*j+2*n;
			A[j][i]=A[i][j];
		}
}
/*
void matPrint(const char * string, double * A, int n, int m)
{
	int i,j;
	printf("\nMatrix name %s\n",string);
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
			printf("%lf ",A[(i*n)+j]);
	
		printf("\n");
	}
}
*/
void vectPrint(const char * string, double * v, int n)
{
	int i;
	printf("\nVector name %s\n",string);
	for( i = 0; i < n; i++)
		printf("%lf ",v[i]);
	printf("\n");
}

double * matToVect(double ** A, int n)
{
	double * res = (double *)malloc(sizeof(double)*n*n);
	int i,j;

	for( i = 0; i < n; i++)
		for( j = 0; j < n; j++)
			res[(i*n)+j] = A[i][j];
	
	return res;
}

double dotProd(double *v, double *w, int t)
{
	double res = 0.;
	int i;

	for( i = 0; i < t; i++)
		res += v[i] * w[i];
	
	return res;
}

double * vectProd(double *v, double *w, int t)
{
	double * res = (double *)malloc(sizeof(double)*t);
	int i;

	for( i = 0; i < t; i++)
		res[i] = v[i] * w[i];
	
	return res;
}

double * matVectProd( double * A, int nbl, int nbc, double * v, int t)
{
	double * res = calloc ( nbl, sizeof( double));
	double somme = 0.;
	int i,j;

	for( i = 0; i < nbl; i++)
	{
		somme = 0.;
		
		for( j = 0; j < nbc; j++)
			somme += v[j] * A[(i*nbc)+j];
			
		res[i] = somme;
	}
	
	return res;
}

double * matMatProd( double * A, double * B, int nbl, int nbc)
{
	
}

double * vectSum(double *v, double *w, int size)
{
	int i;
	double * res=(double *)malloc(size*sizeof(double));
	
	for(i=0;i<size;i++)
		res[i]=v[i]+w[i];
		
	return res;
}

double * vectSub(double *v, double *w, int size)
{
	int i;
	double * res=(double *)malloc(size*sizeof(double));
	
	for(i=0;i<size;i++)
		res[i]=v[i]-w[i];
		
	return res;
}

void scalarVectprod(double * v, double scalar, int size, double * res)
{
	int i;
	
	for(i=0;i<size;i++)
		res[i] = scalar * v[i];
}

void twoD_mat(double ** A, int nbl, int nbc)
{
	int i;
	A = (double **) malloc (nbl*sizeof(double*));
	
	for(i=0;i<nbl;i++)
		A[i] = (double *) malloc (nbc*sizeof(double));
}

void twoD_oneD_Row(double ** A, double * B, int lines, int columns)
{
	int i,j,k;
	k=0;
	
	for(i=0;i<lines;i++)
		for(j=0;j<columns;j++)
		{
			B[k]=A[i][j];
			k++;
		}
}

void twoD_oneD_Col(double ** A, double * B, int lines, int columns)
{
	int i,j,k;
	k=0;
	
	for(i=0;i<lines;i++)
		for(j=0;j<columns;j++)
		{	
			B[k]=A[j][i];
			k++;
		}
}

void minimumVal(double * array, int size, int * out)
{
	double minimum;
	int c;
	*(out) = 0;
	
	minimum = array[0];
 
    for ( c = 0 ; c < size ; c++ ) 
    {
        if ( array[c] < minimum ) 
        {
           minimum = array[c];
           *(out) = c;
        }
    } 
}

void init_2D_Matrix(double ** Matrix, int size)
{
	int i;
	
	Matrix = (double **) malloc (sizeof(double*)*size);
	
	for(i=0;i<size;i++)
		Matrix[i] = (double *) malloc (sizeof(double)*size);
}

void insertVect(double * Mat, double * v, int location, int size, int stride)
{
	int i;
	
	for(i=0;i<size;i++)
	{
		Mat[location] = v[i];
		location += stride;
	}
}

void extractVect(double * Mat, double * v, int location, int size, int stride)
{
	int i;

	for(i=0;i<size;i++)
	{
		v[i]=Mat[location];
		location+=stride;
	}
}
