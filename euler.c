#include<stdlib.h>
#include<math.h>
#include<stdio.h>



const int N = 2;            // number of components of ODE vector
const int num_steps = 10000;  // number of Runge-Kutta steps
const double tau = 0.01;     // stepsize
double y[N][num_steps+1];   // discretization



double p0;                  // initial pressure
double y_0[N] = {p0, 0.0};  // initial pressure & mass



double eos(double p){

  return pow(p/10., 3./5.);
}

double tov(double p, double m, double r){

  return - (eos(p)+ p) * (m + 4*M_PI*pow(r,3.) * p) / (-2*m*r+pow(r, 2.0));
}



int main()
{

  int i1, i2;

  for(int P = 0; P < 10001; P+=1  ){    // initial pressure loop

  double r = 0.000001;
  double p0 = 0.000001 + P * 0.00001;
  double y_0[N] = {p0, 0.0}; 

  for(i1 = 0; i1 < N; i1++)
    {
       y[i1][0] = y_0[i1];
    //   printf("y[i1][0] = %5.6lf\n", y[i1][0]);    // debug statement
    }
  for(i1 = 1; i1<=num_steps ; i1++)  // Euler algorithm (Runge Kutta 1st order)
    {
      double rho = r;
  
     // k1 = f(y(t) , t) * tau
      
      double k1[N];

      k1[0] = tov(y[0][i1-1], y[1][i1-1], r) * tau;
      k1[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1]) * tau;

      for(i2 = 0; i2 < N; i2++)
	y[i2][i1] = y[i2][i1-1] + k1[i2];
    
      r = rho + tau;
    } 

  for(i1 = 0; i1 <= num_steps && y[0][i1] > 0.0 ; i1++)
    {
      double radius = i1 * tau;
     // printf("%9.6lf,%9.6lf\n", radius, y[1][i1]);   // P(r) for one star
    }
 // printf("**********************************************************\n");   // debug
  
  printf("%5.8lf,%5.8lf\n",(i1-1)*tau, y[1][i1-1]); // mass-radius for multiple stars

}

  
  return 0;
}
