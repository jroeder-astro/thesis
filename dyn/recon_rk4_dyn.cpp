#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;


// ***********
// * Physics *
// ***********


// Number of ODE system components
const int N = 2;

// Initial condition, y(t=0)
// y_0[N] = {p0, m0}

// double p0 = 0.000001 + 1000 * 0.0000001;
// double y_0[N] = { p0 , 0.0 };

// Equation of state
double eos(double p){
  return pow(p/10., 3./5.);
}


// *****************
// * RK parameters *
// *****************

const int order = 4; // RK order
const int num_steps_max = 100000000; // Max. RK steps
const double t_max = 10.0; // radius limit (needed?)
const double delta_abs_max = 0.001; // Max. error/step
double tau = 0.01; // Initial coarse stepsize


// double t[num_steps_max+1];     // Discrete radius "axis"
// double y[num_steps_max+1][N];  // Discrete ODE solution

// Calculates f(y(t),t) * tau.
void f_times_tau(double *y_t, double t, 
                 double *f_times_tau_, double tau);

// Do Runge-Kutta
void RK_step(double *y_t, double t, 
             double *y_t_plus_tau, double tau);


// **********
// * Memory *
// **********

void y_alloc(int num_steps_max, int N, double ***y);

void t_alloc(int num_steps_max, double **t);

void y_free(int num_steps_max, int N, double ***y);

void t_free(int num_steps_max, double **t);


// *****************
// * Main function *
// *****************

int main(/*int argc, char **argv*/)
{

  double d1, rho = 0.001;
  int i1, i2;

// *************************************

// Memory allocation for t "radius" axis (discrete)
// Also for te discrete ODE solution array

  vector<double> Eresult;
  vector<double> Presult;

  double *t;
  t_alloc(num_steps_max, &t);

  double **y;
  y_alloc(num_steps_max, N, &y);

// ****************************************

// Read MR input file 

  double M, R;

  FILE *TOV = fopen("MR.dat", "r");
  if (TOV == NULL) exit(0);
  int mcount = 0;

  while (M <= 1.28)
  {
     if(fscanf(TOV, "%lf,%lf", &M, &R) == EOF) break;
     mcount++;
  }
  mcount--;


// ***************************************

// Build initial EoS output vectors

  double p0 = 0.00001 + mcount * 0.00001;
  double y_0[N] = { p0 , 0.0 };

  // Initialize solution.

  t[0] = 0.000000001;

  for(i1 = 0; i1 < N; i1++)
    y[0][i1] = y_0[i1];

  // Do Runge-Kutta.

  //cout << "reached RK" << endl;

  for(i1 = 0;/* i1 < num_steps_max &&*/ 
		y[i1][0] > 0.0; i1++)
    {
      // RK steps

      double y_tau[N];

      // y(t)   --> \tau   y_{\tau}(t+\tau)
      RK_step(y[i1], t[i1], y_tau, tau);

      //cout << "about to vary stepsize" << endl;

      if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -9.))
	 // delta_abs <= delta_abs_max
	{
	  // Accepting the step

	  for(i2 = 0; i2 < N; i2++)
	    y[i1+1][i2] = y_tau[i2];

	  t[i1+1] = t[i1] + tau;
          
          Presult.push_back(y[i1][0]);
          Eresult.push_back(eos(y[i1][0]));
          // Building vectors
         
//	  tau = rho;
	}

      else
	// Adapt step size so that pressures 
	// are roughly evenly spaced.
	{
	  if( y_tau[0] > y[i1][0])
	  {
	     tau *= 1.1;
	     i1--;
	  }

	  if( y_tau[0] < y[i1][0])
	  {
	     tau /= 1.1;
	     i1--;
	  } 
	 // tau = tau_new;

	 // i1--;
	}
     }
/*
  for (i2 = 0; i2 <= num_steps_max+1; i2++)
  {
    free(y[i2]);
  }
  free(y); 
  y = NULL;

  free(t);
  t = NULL;

*/
  t_free(num_steps_max, &t);
  y_free(num_steps_max, N , &y);

  int num_steps = i1;

  // ***********************************************

  // Reconstruction algorithm

  double err = 0.0001;
  int num = 3;
  double DP = 0.0;
  double reos;
  double A;
  
  for (int it = 1; i1 <= num; it++)
  {
    DP += (Presult[i1]-Presult[it+1]);
  }

  double DP_av = DP/num;
  double NP = (Presult[0]-Presult[1])/DP_av;
  double **REOS = NULL;

  while (M <= 2.3)
  { 
    if (fscanf(TOV, "%lf,%lf", &M, &R) == EOF) break;
    reos = Eresult[0];

    while (fabs(y[][] - M) > err)
    { 
    
      // *********************************************    

      // Interpolation of EoS, 2D array
      
      y_alloc(NP, 2, &REOS);
     
      A = (reos-Eresult[0])/(Presult[0]-Presult[1]);
  
      for (int n = 0; n <= NP; n++)
      {
        REOS[n][0] = Presult[0] - n*DP_av;
        REOS[n][1] = A * (Presult[0] - n*DP_av) + 
		     (Eresult[0] - Presult[1]*A);
      }

      // *********************************************
  
      // Runge Kutta is back!

      double *t;
      t_alloc(num_steps_max, &t);

      double **y;
      y_alloc(num_steps_max, N, &y);


    }
  }




  fclose(TOV); 
  TOV = NULL;

  // Output
/*
  for (i1 = 0; i1 <= Mresult.size()-1; i1++)
  {
    printf("%2.15lf,%2.15lf\n", Rresult[i1], 
                                Mresult[i1]/1.4766);
  }
*/

 return EXIT_SUCCESS;
}

// *****************************************************



// *************
// * Functions *
// *************


// Allocating memory

void t_alloc(int num_steps_max,  double **t)
{
   *t = (double*)malloc((num_steps_max+1)*sizeof(double));
   if (*t == NULL)
   {
     cout << "Fehler! t_alloc" << endl; 
     exit(0);
   }  
}

void y_alloc(int num_steps_max, int N, double ***y)
{
  int A = num_steps_max + 1;
  if((*y = (double **)malloc(A*sizeof(double *))) == NULL)
 // if (*y == NULL);
  {
    cout << "Fehler! y_alloc 1. Instanz" << endl; 
    exit(0);
  }

  for (int i1 = 0; i1 <= num_steps_max+1; i1++)
  {
    (*y)[i1] = (double*)malloc(N*sizeof(double));
    if ((*y)[i1] == NULL)
    {
      cout << "Fehler! y_alloc 2. Instanz" << endl; 
      exit(0);
    }
  }
}


// Freeing the allocated memory

void t_free(int num_steps_max, double **t)
{
  free(*t);
  *t = NULL; 
}

void y_free(int num_steps_max, int N, double ***y)
{
  for (int i2 = 0; i2 <= num_steps_max+1; i2++)
  {
    free((*y)[i2]);
  }

  free(**y); 
  **y = NULL;
}


// RK4 function

void RK_step(double *y_t, double t, 
             double *y_t_plus_tau, double tau)
{
  int i1;

  // Calculate k1 = f(y(t),t) * tau.

  double k1[N];
  f_times_tau(y_t, t, k1, tau);

  // k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau.

  double y_2[N];

  for(i1 = 0; i1 < N; i1++)
    y_2[i1] = y_t[i1] + 0.5*k1[i1];

  double k2[N];
  f_times_tau(y_2, t + 0.5*tau, k2, tau);

  // k3

  double y_3[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_3[i1] = y_t[i1] + 0.5*k2[i1];

  double k3[N];
  f_times_tau(y_3, t + 0.5*tau, k3, tau);

  // k4

  double y_4[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_4[i1] = y_t[i1] + k3[i1];

  double k4[N];
  f_times_tau(y_4, t + tau, k3, tau);

  // final

  for(i1 = 0; i1 < N; i1++){
       y_t_plus_tau[i1] = y_t[i1] + 1./6. * 
       (k1[i1] + 2.*k2[i1] + 2.*k3[i1] + k4[i1]);
    }
}


// Calculates f(y(t),t) * tau

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau_, double tau)
{
  if(N != 2)
    {
      fprintf(stderr, "Error: N != 2!\n");
      exit(EXIT_FAILURE);
    }

  // TOV equation implementation
  f_times_tau_[0] = tau * (  -(eos(y_t[0]) + y_t[0]) 
                        * (y_t[1] + 4*M_PI*pow(t, 3.)*y_t[0]) 
                        / (-2*y_t[1]*t + pow(t, 2.0)));
  f_times_tau_[1] = tau * 4*M_PI * pow(t, 2.) *eos(y_t[0]);
}


// *****************************************************
 
/*    // Let's keep this here just to be sure.

	if ((REOS = (double**)malloc(NP*
	    sizeof(double*))) == NULL)
	{
	 cout << "Error" << endl;
	 exit(0); 
	} 

	for (i1 = 0; i1 < NP; i1++)
	{
	if ((REOS[i1] = (double*)malloc(2 *
		  sizeof(double))) == NULL)
	{
	  cout << "Error" << endl;
	  exit(0); 
	} 
	}
*/
/*    // Also, just to be sure.

	double* t;
	t = (double*)malloc((num_steps_max+1)*sizeof(double));
	if(t == NULL)
	{
	cout << "Fehler!" << endl; 
	exit(0);
	}

	double** y;
	y = (double**)malloc((num_steps_max+1)*sizeof(double*));
	if (y == NULL)
	{
	cout << "Fehler!" << endl; 
	exit(0);
	}
	for(i1 = 0; i1 <= num_steps_max+1; i1++)
	{
	y[i1] = (double*)malloc(N*sizeof(double));
	if (y[i1] == NULL)
	{
	cout << "Fehler!" << endl; 
	exit(0);
	}
	}
*/

