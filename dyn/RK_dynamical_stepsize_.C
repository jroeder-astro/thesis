#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector.h>


// ********************
// Physics
// ********************


// TOV Eqs

double eos(double p)
{
  return pow(p/10., 3./5.);
}

double tov(double p, double m, double r)
{
  return - ((eos(p) + p) * (m + 4*M_PI*pow(r, 3.) * p)
         / (pow(r, 2.) - 2*m*r));
}


const int N = 2;  // Number of ODE components


// const int n = 2;
// const double alpha = 0.5;
// const int n = 20;
// const double alpha = 1.0;

double y_0[N] = { 1.0 , 0.0 };  // Starting condition y(t=0).


// Calculates f(y(t),t) * tau.

void f_times_tau(double *y_r, double r, 
                 double *f_times_tau_, double tau)
{
  if(N != 2)
    {
      fprintf(stderr, "Error: N != 2!\n");
      exit(EXIT_FAILURE);
    }

  f_times_tau_[0] = tov(y_r[0], y_r[1], r)  * tau;
  f_times_tau_[1] = 4*M_PI*pow(r, 2.0) * eos(y_r[0]) * tau;
}


// ********************
// RK parameters
// ********************


const int order = 4;


// Max. number of RK steps
const int num_steps_max = 10000;

// Above t_max, break
const double r_max = 20.0;

// Max. allowed error per step
const double d_abs_max = 0.001;

// Initial coarse stepsize
double tau = 1.0;


// **********


double r[num_steps_max+1];  // Diskretisierte Zeitachse.
double y[num_steps_max+1][N];  // Diskretisierte "Bahnkurven".


// **********

// RK step (2nd order) with size tau

void RK_step(double *y_r, double r, 
             double *y_r_plus_tau, double tau)
{
  int i1;

  // Calculate k1 = f(y(t),t) * tau.

  double k1[N];
  f_times_tau(y_r, r, k1, tau);

  // k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau.

  double y_2[N];

  for(i1 = 0; i1 < N; i1++)
    y_2[i1] = y_r[i1] + 0.5*k1[i1];

  double k2[N];
  f_times_tau(y_2, r + 0.5*tau, k2, tau);

  // k3

  double y_3[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_3[i1] = y_r[i1] + 0.5*k2[i1];

  double k3[N];
  f_times_tau(y_3, r + 0.5*tau, k3, tau);

  // k4

  double y_4[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_4[i1] = y_r[i1] + k3[i1];

  double k4[N];
  f_times_tau(y_4, r + tau, k3, tau);

  // final

  for(i1 = 0; i1 < N; i1++){
       y_r_plus_tau[i1] = y_r[i1] + 1./6. * 
       (k1[i2] + 2.*k2[i1] + 2.*k3[i2] + k4[i2]);
    }
}

// **********


int main(/*int argc, char **argv*/)
{
  double d1;
  int i1, i2;

  // Initialisiere "Bahnkurven" mit Anfangsbedingungen.

  t[0] = 0.0;

  for(i1 = 0; i1 < N; i1++)
    y[0][i1] = y_0[i1];


  // Führe Euler/RK-Schritte aus.


  for(i1 = 0; i1 < num_steps_max; i1++)
    {
      if(t[i1] >= t_max) break;

      // RK steps

      double y_tau[N], y_tmp[N], y_2_x_tau_over_2[N];

      // y(t)   --> \tau   y_{\tau}(t+\tau)
      RK_step(y[i1], t[i1], y_tau, tau);

      // y(t)   --> \tau/2   --> \tau_2   y_{2 * \tau / 2}(t+\tau)
      RK_step(y[i1], t[i1], y_tmp, 0.5*tau);
      RK_step(y_tmp, t[i1]+0.5*tau, y_2_x_tau_over_2, 0.5*tau);

      // *****

      // Fehlerabschätzung.

      double delta_abs = fabs(y_2_x_tau_over_2[0] - y_tau[0]);

      for(i2 = 1; i2 < N; i2++)
	{
	  d1 = fabs(y_2_x_tau_over_2[i2] - y_tau[i2]);

	  if(d1 > delta_abs)
	    delta_abs = d1;
	}

      delta_abs /= pow(2.0, (double)order) - 1.0;

      // *****

      // Schrittweitenanpassung (maximale Veränderung um Faktor 5.0).

      d1 = 0.9 * pow(delta_abs_max / delta_abs, 1.0 / (((double)order)+1.0));

      if(d1 < 0.2)
	d1 = 0.2;

      if(d1 > 5.0)
	d1 = 5.0;

      double tau_new = d1 * tau;

      // *****

      if(delta_abs <= delta_abs_max)
	{
	  // Akzeptieren des RK-Schrittes.

	  for(i2 = 0; i2 < N; i2++)
	    y[i1+1][i2] = y_2_x_tau_over_2[i2];

	  t[i1+1] = t[i1] + tau;

	  tau = tau_new;
	}
      else
	// Wiederholen des RK-Schrittes mit kleinerer Schrittweite.
	{
	  tau = tau_new;

	  i1--;
	}
    }

  int num_steps = i1;

  // Output

  for(i1 = 0; i1 <= num_steps; i1++)
    {
      printf("%9.6lf,%9.6lf\n", t[i1], y[i1][0]);
    }
 
  return EXIT_SUCCESS;
}

