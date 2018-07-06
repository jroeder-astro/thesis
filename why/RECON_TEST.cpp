#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <array>

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
const int num_steps_max = 1000000; // Max. RK steps
const double t_max = 10.0; // radius limit (needed?)
const double delta_abs_max = 0.001; // Max. error/step
double tau = 0.01; // Initial coarse stepsize

// double t[num_steps_max+1];     // Discrete radius "axis"
// double y[num_steps_max+1][N];  // Discrete ODE solution

// Calculates f(y(t),t) * tau.
void f_times_tau(double *y_t, double t, 
                 double *f_times_tau_, double tau);

// Do I really need this?
void f_times_tau_rec(double *y_t, double t, 
                     double *f_times_tau_, double tau, 
                     double *recon);

// Do Runge-Kutta
void RK_step(double *y_t, double t, 
             double *y_t_plus_tau, double tau, 
             vector<double> EOS_tmp);

// Do Euler
void Euler_step(double *y_t, double t, 
       double *y_t_plus_tau, double tau, double *recon);


// **********
// * Memory *
// **********

void two_alloc(int num_steps_max, int N, double ***y);
void one_alloc(int num_steps_max, double **t);
void two_free(int num_steps_max, int N, double ***y);
void one_free(int num_steps_max, double **t);


// *****************
// * Main function *
// *****************

int main(/*int argc, char **argv*/)
{

  double d1, rho = 0.001;
  int i1, i2;

// *************************************

// Memory allocation for t "radius" axis (discrete)
// Also for the discrete ODE solution array

  vector<vector<double>> EoS;
  vector<double> Eresult;
  vector<double> Presult;  

  vector<array<double, N>> EOS_arr;
  array<double, N> eos_tmp;

  double *t;
  one_alloc(num_steps_max, &t);

  double **y;
  two_alloc(num_steps_max, N, &y);

// ****************************************

// Read MR input file 

  double M = 0.0, R = 0.0;

  FILE *TOV = fopen("MR.dat", "r");
  if (TOV == NULL) exit(0);
  int mcount = 0;

  while (M <= 1.28)
  {
     if(fscanf(TOV, "%lf,%lf", &R, &M) == EOF) break;
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
    

  double* y_tau;
  one_alloc(N, &y_tau);


  for(i1 = 0;/* i1 < num_steps_max &&*/ 
		y[i1][0] > 0.0; i1++)
    {
      // RK steps
      
      vector<double> EOS_tmp;
      vector<double> result_tmp;      

      // double* y_tau;
      // one_alloc(N, &y_tau);

      // y(t)   --> \tau   y_{\tau}(t+\tau)
   
      RK_step(y[i1], t[i1], y_tau, tau, EOS_tmp);

      //cout << "about to vary stepsize" << endl;

      // ******************************************
   
      // Check pressure difference between array points


/*
      if ( fabs(y_tau[0]-y[i1][0]) <= pow(10., -5.5) )
      {

          cout << "hi" << endl;

	  if( y_tau[0] > y[i1][0])
	  {
	     tau /= 1.1;
	     i1--;
	  }

	  if( y_tau[0] < y[i1][0])
	  {
	     tau *= 1.1;
	     i1--;
	  } 
      }

      else if ( fabs(y_tau[0]-y[i1][0]) >= pow(10., -4.5) )
      {

          cout << "yo2" <<endl;

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
      }

      else
      {
          cout << "accept" << endl;

     	  // Accepting the step
       
	  for(i2 = 0; i2 < N; i2++)
	    y[i1+1][i2] = y_tau[i2];

	  t[i1+1] = t[i1] + tau;

          // Building the vectors          

          EoS.push_back(EOS_tmp);
 
          eos_tmp = {y[i1][0], eos(y[i1][0])};
          EOS_arr.push_back(eos_tmp);
      }


*/ // FROM HERE, ONLY UPPER BOUND...


      if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -5.))	
      {
	  // Accepting the step

	  for(i2 = 0; i2 < N; i2++)
	    y[i1+1][i2] = y_tau[i2];

	  t[i1+1] = t[i1] + tau;

          // Building the vectors          

          EoS.push_back(EOS_tmp);

          // Presult.push_back(y[i1][0]);
          // Eresult.push_back(eos(y[i1][0]));
         
          eos_tmp = {y[i1][0], eos(y[i1][0])};
          EOS_arr.push_back(eos_tmp);

          // tau = rho;
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

      // one_free(N, &y_tau);
      // y_tau = NULL;    
 
    }

  one_free(N, &y_tau);
  y_tau = NULL;    
 
  one_free(num_steps_max, &t);

  cout << endl; // WHAT IS THIS FOR??

  two_free(num_steps_max, N , &y);



//  int num_steps = i1;
//  cout << "hello" << endl;    
//  cout << Presult[i1] << "  " << Eresult[i1] << endl; 

//  vector<std::array<double, N>> EOS_arr;


//  double **EoS_arr;
//  two_alloc(Eresult.size(), N, &EoS_arr);


/* // debug
  for(int j1 = 0; j1 < EOS_arr.size(); j1++)
  {  
     // out << Presult[j1] << "        " << Eresult[j1] << endl;     

     // array<double, N> eos { { Presult[j1], Eresult[j1] } };
     // EOS_arr.push_back(eos);
     
      cout << EOS_arr[j1][0] << "  "  << EOS_arr[j1][1] << endl;     

     // EoS_arr[j1][0] = Presult[j1];
     // EoS_arr[j1][1] = Eresult[j1];
  }
*/

//  cout << "end 1\n";

  // ***********************************************

  // Reconstruction algorithm

  double err = 0.0001;
  int num = 5;
  double DP = 0.0;
  double Preos = 0.0;
  double Ereos = 0.0;
  double A = 0.0;

  for (int it = 1; it <= num; it++)
  {
   // DP += (Presult[it]-Presult[it+1]);
   // cout << EOS_arr[it][0] << endl;
    DP += (EOS_arr[it][0] - EOS_arr[it+1][0]);
   // cout << DP << endl;    
   // First loop for interpolation setup 
   // Sum of Delta(epsilon)'s for mean D(eps)
  }


  double DP_av = DP/num;
//  double NPx = (EOS_arr[0][0]-EOS_arr[1][0])/DP_av;
  double NPx = 0.00001/DP_av; 
  int NP = (int)NPx; // typecast for later loop
//  double NP = (Presult[0]-Presult[1])/DP_av;
  double **REOS;


// control output
  cout << EOS_arr[0][0]-EOS_arr[1][0] << endl;
  cout << DP_av << endl;
  cout << NPx << endl; 
  cout << NP << endl;



// ******************************************************
// ******* Segmentation fault issues from here on *******
// ******************************************************



// allocation of new time and tov axes
  double *t_2;
  one_alloc(5, &t_2);

  double **y_2;
  two_alloc(8, N, &y_2);


// debug statement
  cout << "about to start interpolating\n";

// allocating the memory for interpolated line
  two_alloc(5, N, &REOS);

// First loop: only go up to mass 2.3 (whatever the unit)
  while (M <= 2.3)
  { 
    if (fscanf(TOV, "%lf,%lf", &M, &R) == EOF) break;
cout << "end 7\n";
 
    // interpolation "anchor", D(e) from last 
    Preos = EOS_arr[0][0] + 0.00001;
    Ereos = EOS_arr[0][1];
cout << "end 7\n";
 
    // Second loop: "error" checking, keeping D(p) even
    while (fabs(y[i1][1] - M) > err)
    { 
 cout << "end 7\n";
    
      // *********************************************    

      // Linear interpolation of EoS, 2D array
      // two_alloc(NP, N, &REOS);
     
      // A = (Ereos-EOS_arr[0][1])/(EOS_arr[0][0]-EOS_arr[1][0]);
      A = (Ereos - EOS_arr[0][1])/(Preos - EOS_arr[0][0]);
   
      for (int n = 0; n < NP; n++)
      {
        REOS[n][0] = EOS_arr[0][1] - n*DP_av;
        REOS[n][1] = A * (EOS_arr[0][1] - n*DP_av) + 
		     (EOS_arr[0][1] - EOS_arr[0][0]*A);

        cout << REOS[n][0] << " " << REOS[n][1] << endl;

      }
      // *********************************************
  
cout << "end 7\n";
 
      // Try reconstruction with Euler 
      // hopefully not too crappy

      double *y_tau;
      one_alloc(N, &y_tau);

      t_2[0] = 0.000000001;

      // ************************************************
      // PART ONE: Start calculating mass with the 
      // interpolated REOS array


      for(i1 = 0; i1 < num_steps_max &&   // !!
                 y_2[i1][0] > 0.0; i1++)
      {
      //  double *y_tau;
      //  one_alloc(N, &y_tau);

        p0 = EOS_arr[0][0];
        double y_0[N] = { p0 , 0.0 };

      // Initialize solution.

   cout << "end 8\n";

      // t_2[0] = 0.000000001;

        for(i2 = 0; i2 < N; i2++)
            y_2[0][i2] = y_0[i2];

   cout << y_2[i1] << " " << t_2[i1] << " " 
        << y_tau << " " <<  tau << " " << REOS[i1] << endl;


        Euler_step(y_2[i1], t_2[i1], y_tau, tau, REOS[i1]);
 
   cout << y_2[i1][0] << " " << y_2[i1][1] << endl;

   cout << "end 13\n";


      //  one_free(N, &y_tau);
      }

      one_free(N, &y_tau);
 
cout << "end 9\n";

      // **************************************************

      // PART TWO: calculate the rest of the mass with the
      // previously written EOS_arr
     
      for(i1 = 0; i1 < num_steps_max &&   // !!
                 y_2[i1][0] > 0.0; i1++)
      {
        double *y_tau;
        one_alloc(N, &y_tau);

        p0 = EOS_arr[0][0];
        double y_0[N] = { p0 , 0.0 };

        // Initialize solution.

cout << "end 10\n";


        t_2[0] = 0.000000001;

        for(i1 = 0; i1 < N; i1++)
            y_2[0][i1] = y_0[i1];
 
       // Euler_step(y_2[i1], t_2[i1], y_tau, tau, EOS_arr[i1]);

        one_free(N, &y_tau);
      }

cout << "end 11 \n";


    // ********************************************************

    // PART THREE: Check if calculated mass is close to the one
    // in the MR.dat file 
 
    if (fabs(y_2[i1][1] - M) < err)
    {
      // two_realloc(); //write this function
      for (i1 = 0; i1 < NP ;i1++)
      {
        //double A = REOS[j1][0];
        array<double, N> eos2 { { REOS[i1][0], REOS[i1][1] } };  
        EOS_arr.push_back(eos2);
      }
    }

    else 
    {
      if (y_2[i1][1] > M)
      {
         Ereos /= 1.2;
      }
      if (y_2[i1][1] < M)
      {
         Ereos *= 1.2;
      }
    }  
      


    //  two_free(NP, N, &REOS);

    }
  }

  one_free(num_steps_max, &t_2);
  two_free(num_steps_max, N, &y_2);
  two_free(NP, N, &REOS);

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


// Euler function

void Euler_step(double *y_t, double t, double *y_t_plus_tau, 
             double tau, double *recon)
{
  int i1;
  double k1[N];
  
  f_times_tau_rec(y_t, t, k1, tau, recon);

  for(i1 = 0; i1 < N; i1++){
       y_t_plus_tau[i1] = y_t[i1] + k1[i1];
  }

}


// RK4 function

void RK_step(double *y_t, double t, double *y_t_plus_tau, 
             double tau, vector<double> EOS_tmp)
{
  int i1;

  // Calculate k1 = f(y(t),t) * tau.

  double k1[N];
  f_times_tau(y_t, t, k1, tau);

  EOS_tmp.push_back(y_t[0]);
  EOS_tmp.push_back(eos(y_t[0]));  

  // k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau.

  double y_2[N];

  for(i1 = 0; i1 < N; i1++)
    y_2[i1] = y_t[i1] + 0.5*k1[i1];

  double k2[N];
  f_times_tau(y_2, t + 0.5*tau, k2, tau);
  
  EOS_tmp.push_back(y_2[0]);
  EOS_tmp.push_back(eos(y_2[0]));  

  // k3

  double y_3[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_3[i1] = y_t[i1] + 0.5*k2[i1];

  double k3[N];
  f_times_tau(y_3, t + 0.5*tau, k3, tau);
 
  EOS_tmp.push_back(y_3[0]);
  EOS_tmp.push_back(eos(y_3[0]));  

  // k4

  double y_4[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_4[i1] = y_t[i1] + k3[i1];

  double k4[N];
  f_times_tau(y_4, t + tau, k3, tau);

  EOS_tmp.push_back(y_4[0]);
  EOS_tmp.push_back(eos(y_4[0]));  

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

/*
// RK4 reconstruction function

void RK_step_rec(double *y_t, double t, 
        double *y_t_plus_tau, double tau, double *recon)
{
  int i1;

  // Calculate k1 = f(y(t),t) * tau.

  double k1[N];
  f_times_tau_rec(y_t, t, k1, 
                  tau, recon[0], recon[1]);

  // k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau.

  double y_2[N];

  for(i1 = 0; i1 < N; i1++)
    y_2[i1] = y_t[i1] + 0.5*k1[i1];

  double k2[N];
  f_times_tau_rec(y_2, t + 0.5*tau, k2, 
                  tau, recon[2], recon[3]);

  // k3

  double y_3[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_3[i1] = y_t[i1] + 0.5*k2[i1];

  double k3[N];
  f_times_tau_rec(y_3, t + 0.5*tau, k3, 
                  tau, recon[4], recon[5]);

  // k4

  double y_4[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_4[i1] = y_t[i1] + k3[i1];

  double k4[N];
  f_times_tau_rec(y_4, t + tau, k3, 
                  tau, recon[6], recon[7]);

  // final

  for(i1 = 0; i1 < N; i1++){
       y_t_plus_tau[i1] = y_t[i1] + 1./6. * 
       (k1[i1] + 2.*k2[i1] + 2.*k3[i1] + k4[i1]);
    }
}
*/

// Calculates f*tau with EoS array(s)

void f_times_tau_rec(double *y_t, double t, 
            double *f_times_tau_, double tau, 
            double *recon)
{
  if(N != 2)
    {
      fprintf(stderr, "Error: N != 2!\n");
      exit(EXIT_FAILURE);
    }

  // TOV equation implementation
  f_times_tau_[0] = tau * (  -( recon[1] + recon[0]) 
                        * (y_t[1] + 4*M_PI*pow(t, 3.)*recon[0]) 
                        / (-2*y_t[1]*t + pow(t, 2.0)));
  f_times_tau_[1] = tau * 4*M_PI * pow(t, 2.) * recon[1];
}


// Allocating memory

void one_alloc(int num_steps_max,  double **t)
{
   *t = (double*)malloc((num_steps_max+1)*sizeof(double));
   if (*t == NULL)
   {
     cout << "Fehler! t_alloc" << endl; 
     exit(0);
   }  
}

void two_alloc(int num_steps_max, int N, double ***y)
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

void one_free(int num_steps_max, double **t)
{
  free(*t);
  *t = NULL; 
}

void two_free(int num_steps_max, int N, double ***y)
{
  for (int i2 = 0; i2 <= num_steps_max+1; i2++)
  {
    free((*y)[i2]);
  }

  free(**y); 
  **y = NULL;
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

