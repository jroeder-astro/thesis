#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

void one_d_alloc(int N, double **x);
void one_d_free(int N, double **x);
void two_d_alloc(int N, int M, double ***x);
void two_d_free(int N, int M, double ***x);


int main(){


 // 1D alloc test

 int N = 3;
 int i1;
 double* t;
 one_d_alloc(N, &t);
 
 for (i1 = 0; i1 < N; i1++)
 {
   t[i1] = i1;
   cout << t[i1] << endl;
 }

 one_d_free(N, &t);


 // 2D alloc test

 int i2;
 int M = 4;
 double** s;
 two_d_alloc(N, M, &s);

 for (i1 = 0; i1 < N; i1++)
 {
   for (i2 = 0; i2 < M; i2++)
   {
     s[i1][i2] = i1;
     printf("%f \n", s[i1][i2]);
   }
 }

 two_d_free(N,M,&s);


 return 0;
}


void one_d_alloc(int N, double **x)
{
  if ((*x = (double*)malloc(N*sizeof(double))) == NULL)
  {
    cout << "Fehler!" << endl;
    exit(0);
  }
}

void two_d_alloc(int N, int M, double ***x)
{
  int i1;
  if( (*x = (double **)malloc(N*sizeof(double*))) == NULL)
  {
     cout << "Fehler 1. Instanz" << endl;
     exit(0);
  }
  for(i1 = 0; i1 < M; i1++)
  {
    if (((*x)[i1] = (double*)malloc(M*sizeof(double))) == NULL)
    {
      cout << "Fehler 2. Instanz" << endl;
      exit(0);
    }
  }
}

void one_d_free(int N, double **x)
{
  free(*x);
  *x = NULL;
  cout << "success one_d" << endl;
}

void two_d_free(int N, int M, double ***x)
{
  int i1;
  for(i1 = 0; i1 < N; i1++)
  {
    free((*x)[i1]);
  }
  cout << "first instance success" << endl;
  free(**x);
  cout << "almost\n";
  **x = NULL;
}

