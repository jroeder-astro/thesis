#include<stdlib.h>
#include<math.h>
#include<stdio.h>



const int N = 2;            // number of components of ODE vector
const int num_steps = 10000;  // number of Runge-Kutta steps
const double tau = -0.01;     // stepsize
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

//  for(int P = 0; P < 10001; P+=1  ){    // initial pressure loop

  FILE *radii = fopen("radii_rk4.out", "r");
  if (radii == NULL) exit(0);
  FILE *mass = fopen("mass_rk4.out", "r");
  if (mass == NULL) exit(0);
 

  while(1){

    double r = 0.0;
    double m0 = 0.0;
  
    int d1 = fscanf(radii, "%lf", &r);
    if (d1 == EOF) break;
 // printf("r = %lf\n", r);
    int d2 = fscanf(mass, "%lf", &m0);
    if (d2 == EOF) break;
 // printf("m0 = %lf\n", m0);
    
    double R = r;
    double p0 = pow(10., -11.);
    double y_0[N] = {p0, m0}; 


 // while (m0 <= 1.0){

	  for(i1 = 0; i1 < N; i1++)
	    {
	       y[i1][0] = y_0[i1];
	    // printf("y[i1][0] = %5.12lf\n", y[i1][0]);    // debug statement
	    }

	  for(i1 = 1; i1<=num_steps ; i1++)  // actual Runge-Kutta 4th order algorithm
	    {
              // printf("entered rk4 loop\n"); 
              double rho = r;
	  
              // k1 = f(y(t) , t) * tau
	      
	      double k1[N];

	      k1[0] = tov(y[0][i1-1], y[1][i1-1], r) * tau;
	      k1[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1]) * tau;

	      // k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau

	      double k2[N];

	      k2[0] = tov(y[0][i1-1] + (1./2.)*k1[0], y[1][i1-1] + (1./2.)*k1[1], r + tau/2.)  * tau;
	      k2[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1] + (1./2.)*k1[0]) * tau;
	      r = rho;
	      // k3 = f(y(t)+(1/2)*k2 , t+(1/2)*tau) * tau

	      double k3[N];

	      k3[0] = tov(y[0][i1-1] + (1./2.)*k2[0], y[1][i1-1] + (1./2.)*k2[1], r + tau/2.)   * tau;
	      k3[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1] + (1./2.)*k2[0]) * tau;
	      r = rho;
	      // k4 = f(y(t)+k3 , t+tau) * tau

	      double k4[N];

	      k4[0] = tov(y[0][i1-1] + k3[0], y[1][i1-1] + k3[1], r + tau)  * tau;
	      k4[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1] + k3[0])*tau;

	      for(i2 = 0; i2 < N; i2++)
		{
                   y[i2][i1] = y[i2][i1-1] + 1./6. * (k1[i2] + 2.*k2[i2] + 2.*k3[i2] + k4[i2]);
	        // printf("y[i2][i1] = %lf\n", y[i2][i1]);
	        }
      r = rho + tau;
      if (y[0][i1] <= 0.0 || y[1][i1] <= 0.0) break; 	
    }
  
    // printf("i1 = %d\n", i1); 
    // printf("%5.8lf,%5.8lf\n", y[1][i1-1], eos(y[0][i1-1])); // given eos
//}

/*
while (m0 > 1.0){

       }

*/

  for(i1 = 0; i1 <= num_steps; i1++)
    {
      double radius = R + i1 * tau;
      if (radius < 0.0) break;
      if (y[1][i1] < 0.0) break;
     // printf("%9.6lf,%9.6lf\n", radius, y[1][i1]);   // m(r) for one star
    }

  if (R + i1 * tau < 0.0) break;
  if (y[0][i1] < 0.0) break;

//printf("%5.9lf,%5.9lf\n", R + (i1-1)*tau, y[0][i1-1]); // reconstructed eos

  printf("%5.9lf,%5.9lf\n", y[0][i1-1], eos(y[0][i1-1])); // reconstructed eos


 // printf("**********************************************************\n");   // debug
  
 // printf("%5.8lf,%5.8lf\n",(i1-1)*tau, y[1][i1-1]); // mass-radius for multiple stars
 // printf("%5.8lf\n",(i1-1)*tau); printf("%5.8lf\n",y[1][i1-1]);



//}

  }

  fclose(radii);
  fclose(mass);

  return 0;
}
