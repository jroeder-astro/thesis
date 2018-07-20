#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<iostream>
using namespace std;


const int N = 2;            // number of components of ODE vector
const int num_steps = 10000;  // number of Runge-Kutta steps
double tau = 0.5;     // stepsize
double y[N][num_steps+1];   // discretization



double eos(double p){

  return pow(p/10., 3./5.);
}

double tov(double p, double m, double r){

  return - (eos(p)+ p) * (m + 4*M_PI*pow(r,3.) * p) 
         / (-2*m*r+pow(r, 2.0));
}



int main()
{

  int i1, i2;
  int count = 1;
  int order = 4;
  double tau_new = 0.;

  // Three variables needed for pressure variation check
  double d1, d_abs;
  double d_abs_max = 0.001;

  // Variables needed for physics stuff
  double r, p0;
  double y_0[N];

 // for(int P = 0; P < 10001; P += 1){  // initial pressure loop

  r = 0.000001;
  p0 = 0.000001 + 2000 * 0.00001;
  y_0[0] = p0; y_0[1] = 0.0; 

	for(i1 = 0; i1 < N; i1++)
	{
	y[i1][0] = y_0[i1];
        }

	for(i1 = 1; i1<=num_steps ; i1++)  
        // actual Runge-Kutta 4th order algorithm
 	{
	double rho = r;

	// k1 = f(y(t) , t) * tau

	double k1[N];

	k1[0] = tov(y[0][i1-1], y[1][i1-1], r) * tau;
	k1[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1]) * tau;

	// k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau

	double k2[N];

	k2[0] = tov(y[0][i1-1] + (1./2.)*k1[0], y[1][i1-1] 
                + (1./2.)*k1[1], r + tau/2.)  * tau;
	k2[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1] 
                + (1./2.)*k1[0]) * tau;
	r = rho;
	
        // k3 = f(y(t)+(1/2)*k2 , t+(1/2)*tau) * tau

	double k3[N];

	k3[0] = tov(y[0][i1-1] + (1./2.)*k2[0], y[1][i1-1] 
                + (1./2.)*k2[1], r + tau/2.)   * tau;
	k3[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1] 
                + (1./2.)*k2[0]) * tau;
	r = rho;
	
        // k4 = f(y(t)+k3 , t+tau) * itau

	double k4[N];

	k4[0] = tov(y[0][i1-1] + k3[0], y[1][i1-1] 
                + k3[1], r + tau)*tau;
	k4[1] = 4*M_PI*pow(r, 2.0) * eos(y[0][i1-1] + k3[0])*tau;

	for(i2 = 0; i2 < N; i2++){
	y[i2][i1] = y[i2][i1-1] + 1./6. 
                    * (k1[i2] + 2.*k2[i2] + 2.*k3[i2] + k4[i2]);
        }
//	r = rho + tau;

//***************************************************************


cout << "y[0][i1] = " << y[0][i1] << endl;

cout << "y[0][i1-1] = " << y[0][i1-1] << endl;

cout << "y[0][i1+1] = " << y[0][i1+1] << endl;



    d_abs = fabs(y[0][i1-2]-y[0][i1-1])/(pow(2., (double)order)-1.0);
    d1 = 0.9 * pow(d_abs_max/d_abs, 1.0/(((double)order)+1.0));


//   d1 = fabs(y[0][i1]-y[0][i1+1]);

cout << "d1 = " << d1 << endl;
cout << "d_abs = " << d_abs <<endl;
cout << "tau = " << tau << endl;
    

    if (d1 < 0.2) d1 = 0.2;
    if (d1 > 2.0) d1 = 2.0;
    tau_new = d1*tau;
   
   
cout << "tau_new = " << tau_new << endl;

    if (d_abs <= d_abs_max)
      {
//	for(i2 = 0; i2 < N; i2++){
//	y[i2][i1] = y[i2][i1-1] + 1./6. 
//                    * (k1[i2] + 2.*k2[i2] + 2.*k3[i2] + k4[i2]);
//        }
	r = rho + tau;
        tau = tau_new;
      }

    else 
      {
        tau = tau_new;
        i1--;
      }

//***************************************************************
    
     } 

	for(i1 = 0; i1 <= num_steps && y[0][i1] > 0.0 ; i1++)
	    {
	      double radius = i1 * tau;
	      printf("%9.20lf\n", y[1][i1]);   
             // P(r) for one star
            }

// printf("****************************************************\n");   







  
//     printf("%5.8lf,%5.8lf,%d\n",(i1-1)*tau, y[1][i1-1],count); 
               // mass-radius for multiple stars
//     printf("%5.8lf,%5.8lf,%d\n",y[0][i1-1], eos(y[0][i1-1]),count); 
               // mass-radius for multiple stars
//     printf("%5.8lf\n",(i1-1)*tau); printf("%5.8lf\n",y[1][i1-1]);
     count += 1;

//}
  
  return 0;
}
