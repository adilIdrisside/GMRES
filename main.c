#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
//#include <lapacke.h>
//#include "../headers/functions.h"
#include "functions.h"
#include "gmres.h"
#include "matrix.h"


extern void mat_init(double * A, int size);
extern double * init_b(double * b,int size);
extern void poissonMat(double * A, int size);
extern void matPrint(const char * string, double * A, int n, int m);
extern void gmres(double * A, double * x, double * b, double * r, double * q, double * Hessenberg, double * Krylov, int size, int m);




int main(int argc, char ** argv)
{
	int i, m, size;
	
	if(argc != 3)
	{
		perror("USAGE: ./gmres <m> <Matrix size>");
		exit(0);
	}else{
		m = atoi(argv[1]);
		size = atoi(argv[2]);
	}
	
	//m = 2; size = 4;
	/* Matrix */
	double * A = (double *)malloc(sizeof(double)*size*size);
	double * B = (double *)malloc(sizeof(double)*size*size);
	double * Krylov = (double *)malloc(sizeof(double)*size*m);
	double * Hessenberg = (double *)malloc(sizeof(double)*(m+1)*m);
	
	/* Vectors */
	double * x = (double *)malloc(sizeof(double)*size);
	double * b = (double *)malloc(sizeof(double)*size);
	double * r = (double *)malloc(sizeof(double)*size);
	double * q = (double *)malloc(sizeof(double)*size);
	
	printf("m:%d, N:%d\n",m,size);
	printf("A Matrix size: (%dx%d)\n",size,size);
	printf("Hessenberg Matrix size: (%dx%d)\n",m+1,m);
	printf("Krylov sub-space size: (%dx%d)\n",size,m);
	
	//mat_init(A, size);
	nonsymSparseMatInit(A, size);
	matPrint("A", A, size, size);
	
	//for(i=0;i<size;i++)
		//x[i] = 1.0;
	
	x[0] = 1.0;
	b = init_b(b, size);
	vectPrint("b",b,size);
	printf("Initial guess\n");
	vectPrint("x",x,size);
	gmres(A, x, b, r, q, Hessenberg, Krylov, size, m);
	
	return 0;
}
