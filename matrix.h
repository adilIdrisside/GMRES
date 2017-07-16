#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <lapacke.h>
#define RAND_INIT 100

void mat_init(double * A, int size);

void mat_assign(double * A, double * B, int size);

void matPrint(const char * string, double * A, int n, int m);

void matPrint_2D(const char * string, double ** A, int n, int m);

void nonsymSparseMatInit(double * A, int size);

void grabMatFromFile(const char * path, double * A, int size);

#endif //MATRIX_H_INCLUDED

