#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "functions.h"
#include <stdbool.h>
#define EPSILON 1e-3F

double * init_b(double * b,int size); /* initialisation de b */

double min(double a, double b); /* retourne a si a<b, retourne b sinon */

void Hessenberg(double * hessenberg, int size, double h); /* construction de la matrice triangulaire supérieure de Hessenberg */

void minResidualSolver(double beta, double * e1, double * Hessenberg, double * y, int size, double * y_min); /* retourne le résidu avec la plus petite valeur */

void earlyConvergenceTest(bool * convergence, double * x, double * Krylov, double * y, double ** hess, double * Hessenberg, double beta, double * res, int size, int m, double s);

void secondConvergenceTest(double beta, double * e1, double * Hessenberg, double * y, int m, int size, double s);

void solverY(double * e1, double beta, double * Hessenberg, double * y, int m, int size, double s, int * location); /* trouver le plus petit vecteur y */

void gmres(double * A, double * x, double * b, double * r, double * q, double * Hessenberg, double * Krylov, int size, int m); /* Itération de GMRES avec redémarrage implicite */

#endif //FUNCTIONS_H_INCLUDED
