#include "matrix.h"

void mat_init(double * A, int size)
{
	int i,j;
	double random;
	srand(RAND_INIT);
	
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			random = 10.*((double)rand())/RAND_MAX;
			A[(i*size)+j]=random;
		}
	}
}

void matPrint(const char * string, double * A, int n, int m)
{
	int i,j;
	printf("\nMatrix name:\t %s\n",string);
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
			printf("%e ",A[(i*m)+j]);
	
		printf("\n");
	}
	printf("\n");
}

void matPrint_2D(const char * string, double ** A, int n, int m)
{
	int i,j;
	printf("\nMatrix name:\t %s\n",string);
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
			printf("%lf ",A[i][j]);
	
		printf("\n");
	}
	printf("\n");
}

void nonsymSparseMatInit(double * A, int size)
{
	int i,j;

	//for(i=1;i<size;i++)
		//A[(i*size)+i-1]=-1.000000;
	
	for(i=0;i<size;i++)
	{
		A[(i*size)+i]=2.000000;
		A[(i*size)+i+1]=-1.000000;
	}
}		

void mat_assign(double * A, double * B, int size)
{
	int i,j;
	
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			B[(i*size)+j]=A[(i*size)+j];
		}
	}
}

void grabMatFromFile(const char * path, double * A, int size)
{
	FILE * f;
    f = fopen(path, "r");

    //read file into array
    int i,j;

    if (f == NULL)
    {
        printf("Error Reading File\n");
        exit (0);
    }
    
    for (i = 0; i < size; i++)
    {
		for (j = 0; j < size; j++)
		{
			fscanf(f, "%lf,", &A[(i*size)+j] );
		}
    }

	/*
    for (i = 0; i < size; i++)
    {
		for (j = 0; j < size; j++)
		{
			printf("%lf ", A[(i*size)+j]);
		}
		printf("\n");
    }
	*/
    
    fclose(f);	
}	
