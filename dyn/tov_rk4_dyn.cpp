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


// Calculates f(y(t),t) * tau.
void f_times_tau(double *y_t, double t, double *f_times_tau_, double tau)
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


// *****************
// * RK parameters *
// *****************


// #define __RK_2ND__
#define __RK_4TH__

#ifdef __RK_2ND__
const int order = 2;
#endif

#ifdef __RK_4TH__
const int order = 4;
#endif

// Max. number of RK steps
const int num_steps_max = 100000000;

// Limit for radius? Nah......
const double t_max = 10.0;

// Max. allowed error per step (work on this) 
const double delta_abs_max = 0.001;

// Initial coarse stepsize
double tau = 0.01;


// **********


// double t[num_steps_max+1];     // Discrete radius "axis"
// double y[num_steps_max+1][N];  // Discrete ODE solution


// **********

#ifdef __RK_4TH__

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

#endif


// **********


// *****************
// * Main function *
// *****************


int main(/*int argc, char **argv*/)
{

  double d1, rho = 0.001;
  int i1, i2;

//  cout << "main" << endl;

  vector<double> Mresult;
  vector<double> Rresult;

//  double t[num_steps_max+1];     // Discrete radius "axis"
//  double y[num_steps_max+1][N];  // Discrete ODE solution

  double* t;
  t = (double*)malloc((num_steps_max+1)*sizeof(double));
  if(t == NULL)
  {
    cout << "Fehler!" << endl; 
    exit(0);
  }
/*
  else
  {
    cout << "t allocated" << endl;
  }
*/

  double** y;
  y = (double**)malloc((num_steps_max+1)*sizeof(double *));
  if (y == NULL)
  {
    cout << "Fehler!" << endl; 
    exit(0);
  }
/*
  else
  {
    cout << "y allocated" << endl;
  }
*/

  for(i1 = 0; i1 <= num_steps_max+1; i1++)
  {
    y[i1] = (double*)malloc(N*sizeof(double));
    if (y[i1] == NULL)
    {
      cout << "Fehler!" << endl; 
      exit(0);
    }
  }


  for (int P = 0; P <= 100; P++){
/*
	  double* t;
	  t = (double*)malloc((num_steps_max+1)*sizeof(double));
	  if(t == NULL)
          {
	    cout << "Fehler!" << endl; 
	    exit(0);
	  }
          else
          {
            cout << "t allocated" << endl;
          }

	  double** y;
	  y = (double**)malloc((num_steps_max+1)*sizeof(double));
	  if (y == NULL)
          {
	    cout << "Fehler!" << endl; 
	    exit(0);
	  }
          else
          {
            cout << "y allocated" << endl;
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

//	  double t[num_steps_max+1];     // Discrete radius "axis"  
//	  double y[num_steps_max+1][N];  // Discrete ODE solution

//        cout << "P loop" << endl;

	  double p0 = 0.00001 + P * 0.00001;
	  double y_0[N] = { p0 , 0.0 };

	  // Initialize solution.

	  t[0] = 0.000000001;

	  for(i1 = 0; i1 < N; i1++)
	    y[0][i1] = y_0[i1];

	  // *****

	  // Do Runge-Kutta.

//cout << "reached RK" << endl;

	  for(i1 = 0;/* i1 < num_steps_max &&*/ y[i1][0] > 0.0; i1++)
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
	 
	//  cout << "rho = " << rho << endl; 
	//  cout << "tau = " << tau << endl;

	     }

        Mresult.push_back(y[i1-1][1]);
        Rresult.push_back(t[i1-1]);

//cout << "M = " << Mresult[i1] << endl;
//cout << "R = " << Rresult[i1] << endl;

   }

        for (i2 = 0; i2 <= num_steps_max+1; i2++)
        {
          free(y[i2]);
        }
        free(y); 
        y = NULL;
        
        free(t);
        t = NULL;
    
/*        
        for (i2 = 1; i2 < i1; i2++)
        {
          y[i2][0] = 0;
          y[i2][1] = 0;
             t[i2] = 0;
        }           
*/

  int num_steps = i1;

  // *****

  // Output

  for (i1 = 0; i1 <= Mresult.size()-1; i1++)
  {
    printf("%2.15lf,%2.15lf\n", Rresult[i1], Mresult[i1]/1.4766);
  }


/*
  // single star output
  for(i1 = 0; i1 <= num_steps && y[i1][0] > 0.0; i1++)
    {
      printf("%2.6lf,%2.15lf\n", t[i1], y[i1][0]);
    }
*/  

  return EXIT_SUCCESS;
}

// **********


