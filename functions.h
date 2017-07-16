#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <lapacke.h>
#include <stdbool.h>

void vectNorm(double * v, int size, double * norm); /* Norme du vecteur */

void symMatInit(double ** A, int n); /* Matrix initialisation */

void vectPrint(const char * string, double * v, int n); /* Vector printing */

double * matToVect(double ** A, int n); /* Matrix linearisation */

double dotProd(double *v, double *w, int t); /* Produit scalaire */

void scalarVectprod(double * v, double scalar, int size, double * res); /* produit scalaire-vecteur */

void scalarVectdiv(double * v, double scalar, int size, double * res); /* division scalaire-vecteur */

double * vectProd(double *v, double *w, int t); /* Produit vectoriel */

double * vectSum(double *v, double *w, int size);

double * vectSub(double *v, double *w, int size); /* Soustraction de deux vecteurs */

double * matVectProd( double * A, int nbl, int nbc, double * v, int t); /* Produit matrice-vecteur */

double * matMatProd( double * A, double * B, int nbl, int nbc); /* Produit matrice-matrice */

void twoD_mat(double ** A, int nbl, int nbc);

void twoD_oneD_Row(double ** A, double * B, int lines, int columns);

void twoD_oneD_Col(double ** A, double * B, int lines, int columns);

void minimumVal(double * array, int size, int * out);

void init_2D_Matrix(double ** Matrix, int size);

void insertVect(double * Mat, double * v, int location, int size, int stride);

void extractVect(double * Mat, double * v, int location, int size, int stride);

#endif //FUNCTIONS_H_INCLUDED
