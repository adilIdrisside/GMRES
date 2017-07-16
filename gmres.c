#include "gmres.h"
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

extern double dotProd(double *v, double *w, int t);
extern double * vectSum(double *v, double *w, int size);
extern double * vectSub(double *v, double *w, int size);
extern void vectNorm(double * v, int size, double * norm);
extern void minimumVal(double * array, int size, int * out);
extern void vectPrint(const char * string, double * v, int size);
extern void matPrint(const char * string, double * A, int n, int m);
extern void matPrint_2D(const char * string, double ** A, int n, int m);
extern void twoD_oneD_Row(double ** A, double * B, int lines, int columns);
extern void twoD_oneD_Col(double ** A, double * B, int lines, int columns);
extern void scalarVectdiv(double * v, double scalar, int size, double * res);
extern void scalarVectprod(double * v, double scalar, int size, double * res);
extern double * matVectProd( double * A, int nbl, int nbc, double * v, int t);
extern void insertVect(double * Mat, double * v, int location, int size, int m);
extern void extractVect(double * Mat, double * v, int location, int size, int stride);


/*
*	Initialisation du vecteur résultat b
*
*	Return: double * 
*	
*	@Param:
*		double * b : vecteur résultat
*		int size   : taille du vecteur
*
*/

double * init_b(double * b, int size)
{
	int i;
	
	b[0]=1;
	b[size-1]=1;
	for(i=1;i<size-1;i++)
		b[i]=2;
		
	return b;
}


double min(double a, double b)
{
	return (a<b?a:b) ;
}

void solverY(double * e1, double beta, double * Hessenberg, double * y, int m, int size, double s, int * location)
{
	int i,j,k,loc,stride;
	double inter[m+1];
	double val[m+1];
	
	e1[0]=beta;
	//matPrint("Hessenberg",Hessenberg,m+1,m);
	//vectPrint("e1",e1,size);
	j=0; 
	for(k=0;k<m+1;k++)
	{
		extractVect(Hessenberg, inter, k*m, m, 1);
		vectNorm(vectSub(e1,matVectProd(Hessenberg, m+1, m, inter, m), m+1), m+1, &s);
		val[k]=s;
	}
	//vectPrint("val",val,m+1);
	minimumVal(val, m+1, &loc);
	//printf("loc: %d\n",loc);
	
	*(location) = loc;
}

void earlyConvergenceTest(bool * convergence, double * x, double * Krylov, double * y, double ** hess, double * Hessenberg, double beta, double * minResidual, int size, int m, double s)
{
	int i,j; 
	int location;
	y = (double *) malloc (sizeof(double)*m);
	double * eigVect = (double *) malloc (sizeof(double)*m*m);
	double * w = (double *) malloc (sizeof(double)*m);
	double * e1 = (double *) malloc (sizeof(double)*m);
	
	e1[0] = beta;
	
	//matPrint("Hessenberg",Hessenberg,(m+1),m);
	twoD_oneD_Row(hess, Hessenberg, m+1, m);
	//matPrint("Hessenberg",Hessenberg,(m+1),m);
	
	if(s<EPSILON)
	{
		solverY(e1, beta, Hessenberg, y, m, size, s, &location);
		
		for(i=0;i<m;i++)
			y[i]=Hessenberg[i+location];
		
		x = vectSum(x,matVectProd(Krylov,size,m,y,m), size);
		*convergence = true;
	}
}

/*
*	Trouver le vecteur (ligne) de Hessenberg avec la plus petite norme vectorielle
*
*	Return : void
*
*	@Param:
*		double   beta 	    : Norme du résidue r
*		double * e1 	    : Vecteur cannonique (e1=[1,0,...,0])
*		double * Hessenberg : Matrice triangulaire supérieure de Hessenberg
*		double * y	    : vecteur avec norme minimale
*		int 	 m	    : Taille du sous-espace de Krylov
*		int 	 size	    : Taille de la matrice A (size x size)
*		double   s	    : Norme minimale des vecteurs ligne de Hessenberg
*/

void secondConvergenceTest(double beta, double * e1, double * Hessenberg, double * y, int m, int size, double s)
{
	int i,j,k,location;
	double * v = (double *)malloc(m*sizeof(double));

	solverY(e1, beta, Hessenberg, y, m, size, s, &location);

	for(i=0;i<m;i++)
	{
		y[i]=Hessenberg[location];
		location++;
	}
}

/*
*	Algorithme de GMRES(m) avec redémarrage implicite, toutes les matrices sont en 1D
*
*	Return: void
*
*	@Param:
*		double * A	    : Matrice du système d'équations linéaires non-symmétrique 
*		double * x	    : Vecteur solution du système (initialisé avec x:=[1,0,0,0,...,0])
*		double * b	    : Vecteur résultat du système Ax = b
*		double * r	    : Vecteur résiduelle
*		double * q	    : Vecteur utilisé pour construire la matrice accompagnante du sous-espace de Krylov
*		double * Hessenberg : Matrice de Hessenberg
*		double * Krylov	    : Matrice de Krylov
*		int 	 size       : (size x size) taille de la matrice A
*		int	 m	    : (m << size)
*/

void gmres(double * A, double * x, double * b, double * r, double * q, double * Hessenberg, double * Krylov, int size, int m)
{
	int i, j, count, iter, location, restart, krylov_loc, Q_loc;
	double beta, scalar, h, s;
	bool convergence;
	
	double * w = (double *) malloc (sizeof(double)*size);
	double * inter = (double *) malloc (sizeof(double)*size);
	double * y = (double *) malloc (sizeof(double)*m);
	double * e1 = (double *) malloc (sizeof(double)*(m+1));
	double ** hess = (double **) malloc (sizeof(double*)*(m+1));
	double Q[(m+1)*size];
	
	for(i=0;i<m+1;i++)
		hess[i] = (double *) malloc (sizeof(double)*m);
		
	//printf("m+1:%d, m:%d\n",m+1,m);
	convergence = false;
	e1[0]=1;
	restart=0;
	
	while ((!convergence) && (restart <20000))
	{
		iter = 0;
		location = 0;
		krylov_loc=0;
		Q_loc=0;
		r = vectSub(b, matVectProd( A, size, size, x, size), size);
		vectNorm(r, size, &beta);
		scalar = 1/(beta);
		scalarVectprod(r, scalar, size, q);
		insertVect(Q, q, Q_loc, size, m+1);
		Q_loc++;
		insertVect(Krylov, q, krylov_loc, size, m);
		
		for(i=0;i<m;i++)
		{
			w = matVectProd( A, size, size, q, size);
			//vectPrint("w",w,size);
			/*
			if(i < (m-1))
			{
				krylov_loc++;
				insertVect(Krylov, w, krylov_loc, size, m);
			}
			*/
			for(j=0;j<=i;j++)
			{
				double _q[size];
				extractVect(Q, _q, j, size, m+1);
				//vectPrint("_q",_q,size);
				h = dotProd(w, _q, size);
				//printf("h(%d,%d): %lf\n",j,i,h);
				hess[j][i] = h;
				scalarVectprod(_q, h, size, inter);
				w = vectSub(w, inter, size);
				//vectPrint("w",w,size);
			}
			
			vectNorm(w, size, &h);
			hess[i+1][i] = h;
			//printf("h(%d,%d): %lf\n",i+1,i,h);
			scalar = 1/h;
			scalarVectprod(w, scalar, size, q);
			
			if(i < (m-1))
			{
				krylov_loc++;
				//double * _x = (double *)malloc(sizeof(double)*size);
				//_x = matVectProd(A, size, size, q, size);
				insertVect(Krylov, q, krylov_loc, size, m);
				//free(_x);
			}
			
			//vectPrint("q",q,size);
			insertVect(Q, q, Q_loc, size, m+1);
			Q_loc++;
			//vectPrint("q",q,size);
			vectNorm(vectSub(b, matVectProd( A, size, size, x, size), size), size, &s);
			earlyConvergenceTest(&convergence, x, Krylov, y, hess, Hessenberg, beta, y, size, m, s);
			
			//matPrint("Krylov",Krylov,size,m);
			//matPrint_2D("Hessenberg",hess,m+1,m);
		}
		if(convergence)
			return;
		//matPrint("Krylov",Krylov,size,m);
		//matPrint("Q",Q,size,m+1);
		//vectPrint("Krylov", Krylov, size*m);
		twoD_oneD_Row(hess, Hessenberg, m+1, m);
		//matPrint("Hessenberg",Hessenberg,m+1,m);
		secondConvergenceTest(beta, e1, Hessenberg, y, m, size, s);
		
		x = vectSum(x, matVectProd(Krylov, size, m, y, m), size);
		//vectPrint("New x",x,size);
								
		if((fabs(hess[m][m-1])*beta)<=EPSILON)
		{
			convergence = true;
			restart++;
			printf("residue: %e\n",(fabs(hess[m][m-1])*beta));	
			printf("itération: %d\n",restart);
			//fputs(convergence ? "true\n" : "false\n", stdout);
			printf("\n\nAlgorithme de GMRES terminé\n");
			printf("Solution approchée du système\n");
			vectPrint("X",x,size);
		}else{	
			restart++;
			printf("residue: %e\n",(fabs(hess[m][m-1])*beta));	
			printf("itération: %d\n",restart);
			//fputs(convergence ? "true\n" : "false\n", stdout);
		}
	}
}
