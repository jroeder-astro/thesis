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

// Equation of state
double eos(double p, vector<double> *alpha){
  return pow(p/10., 3./5.);
}

// Reconstructed EoS
double line(double p, vector<double> *alpha){
// alpha[i] = {a,b,c,d}
  return ((*alpha)[3]-(*alpha)[1])/((*alpha)[2]-(*alpha)[0]) * p + 
          (*alpha)[1] - (*alpha)[0] * ((*alpha)[3]-
          (*alpha)[1])/((*alpha)[2]-(*alpha)[0]);
}


// *****************
// * RK parameters *
// *****************

const int num_steps_max = 1000000; // Max. RK step
double tau = 0.01;                 // Initial stepsize

// Calculates f(y(t),t) * tau
void f_times_tau(double *y_t, double t, double *f_times_tau_, 
                 double tau, vector<double> *alpha, 
                 double(*state)(double, vector<double> *));

// Do Runge-Kutta
void RK_step(double *y_t, double t, double *y_t_plus_tau, 
             double tau, vector<double> *alpha, 
             double(*state)(double, vector<double> *));

/* Do Euler (not really needed)
void Euler_step(double *y_t, double t, 
       double *y_t_plus_tau, double tau, double *recon);
*/


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

int main(){

  double d1, rho = 0.001;
  int i1, i2;
  double P = 0.0;
  double Preos = 0.0;
  double Ereos = 0.0;
  double err = 0.001;

// *************************************

// Memory allocation for t "radius" axis (discrete)
// Also for the discrete ODE solution array

//  vector<vector<double>> EoS;
//  vector<double> Eresult;
//  vector<double> Presult;  

  vector<array<double, N>> EOS_arr;
  array<double, N> eos_tmp;
  
  vector<vector<double>> recon_storage;
  vector<double> alpha;
     for(int x = 0; x < 4; x++){
          alpha.push_back(0);
     }
  vector<double> zero;

  double *t;
  one_alloc(num_steps_max, &t);

  double **y;
  two_alloc(num_steps_max, N, &y);

  vector<vector<double>> MR_rel;
  vector<double> one_MR;  

// ****************************************

// Read MR input file 

  double M = 0.0, R = 0.0;

  FILE *TOV = fopen("MR.dat", "r");
  if (TOV == NULL) exit(0);
  int mcount = 0;

  one_MR.push_back(0); one_MR.push_back(0);

  while (1)
  {
     if(fscanf(TOV, "%lf,%lf", &R, &M) == EOF) break;
     one_MR[0] = R; one_MR[1] = M;
     MR_rel.push_back(one_MR);
  }
  fclose(TOV);
  TOV = NULL;

  for (i2 = 0; i2 <= MR_rel.size(); i2++) 
  {
     if (MR_rel[i2][1] >= 1.28) break;
     mcount++;
  }
  mcount--;

// ***************************************

  double p_init = 0.00001 + mcount * 0.00001;
  double y_0[N] = { P , 0.0 };
  double p_dur = 0.0;
 
  mcount++;

  // Initialize solution.

  t[0] = 0.000000001;

  for(i1 = 0; i1 < N; i1++)
    {
      y[0][i1] = y_0[i1];
    }

  alpha[0] = p_init;  alpha[1] = eos(p_init, &zero);
  alpha[2] = Preos;   alpha[3] = Ereos;


// ***************************************


  // Do Runge-Kutta.

  double* y_tau;
  one_alloc(N, &y_tau);


  while(mcount <= MR_rel.size())
  {
      // Probably initiallization has to happen here!!

  for(i1 = 0; y[i1][0] > 0.0; i1++)
    {
      /* 
         I should probably write a function for doing an RK step
         with a certain eos or line, also checking the stepsize 
         (it is not convenient to have the same piece of code three
         times in a row)   
      */

      // *****************************************************
      // PART I.1: Reconstruction with an "interpolated" line

      while(y[i1][0] > p_dur)
       {  
	      // RK steps
	     
	      RK_step(y[i1], t[i1], y_tau, tau, &alpha, line);

	     // FROM HERE, ONLY UPPER BOUND...

	      if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -5.))	
	      {
		  // Accepting the step

		  for(i2 = 0; i2 < N; i2++)
		    y[i1+1][i2] = y_tau[i2];

		  t[i1+1] = t[i1] + tau;

		  // Building the vectors          

		  // EoS.push_back(EOS_tmp);
	    
		  // eos_tmp = {y[i1][0], line(y[i1][0], &alpha)};
		  // EOS_arr.push_back(eos_tmp);

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
	      }

         // ********************************************
       }

       // ****************************************************
       // PART I.2: Have fun with previous lines 


      if(recon_storage.size() > 1) 
      {
         for(int store = 0; store < recon_storage.size(); store++)
         {    
             /*
                How do I best determine the pressure intervals?
             */ 


            // RK steps
	     
	     RK_step(y[i1], t[i1], y_tau, tau, 
                     &recon_storage[store], line);


         }
      }

       
       // *************************************************
       // PART II: Calculate the rest with known eos
 
	      // RK steps
	     
	      RK_step(y[i1], t[i1], y_tau, tau, &alpha, eos);
  
	     // FROM HERE, ONLY UPPER BOUND...

	      if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -5.))	
	      {
		  // Accepting the step

		  for(i2 = 0; i2 < N; i2++)
		    y[i1+1][i2] = y_tau[i2];

		  t[i1+1] = t[i1] + tau;

		  // Building the vectors          

		  // EoS.push_back(EOS_tmp);
	    
		  // eos_tmp = {y[i1][0], eos(y[i1][0])};
		  // EOS_arr.push_back(eos_tmp);

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
	      }

         // **********************************************
         // PART III: Checking if calculated mass is ok

	   if (fabs(y[i1][1] - MR_rel[mcount][1]) < err)
	   {
 
              recon_storage.push_back(alpha);              
   
	      // two_realloc(); //write this function
	      // for (i2 = 0; i2 < NP; i2++)
	      // {

                 // ******************************************
                 // Write here the slope and intersect array!!
                 // ******************************************

		 // double A = REOS[j1][0];
		 // array<double, N> eos2 { { REOS[i1][0], REOS[i1][1] } };  
		 // EOS_arr.push_back(eos2);
	      // }
	    }

	    else 
	    {
	      if (y[i1][1] > MR_rel[mcount][1])
	      {
		 alpha[3] /= 1.2;
	      }
	      if (y[i1][1] < MR_rel[mcount][1])
	      {
		 alpha[3] *= 1.2;
	      }
              i1 = 0;
	    }  

         // ************************************************
     }

     alpha[0] = Preos; alpha [1] = Ereos;
     Preos += /* xxx */;
     alpha[2] = Preos; alpha[3] = Ereos;

     /*
         What to do with p_dur?
     */

   
     mcount++;

   }


  one_free(N, &y_tau);
  y_tau = NULL;    
 
  one_free(num_steps_max, &t);

  cout << endl; // WHAT IS THIS FOR??

  two_free(num_steps_max, N , &y);


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

/*


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
*/
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
/*
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
*/

// RK4 function

void RK_step(double *y_t, double t, double *y_t_plus_tau, 
             double tau, vector<double> *alpha, 
             double(*state)(double, vector<double> *))
{
  int i1;

  // Calculate k1 = f(y(t),t) * tau.

  double k1[N];
  f_times_tau(y_t, t, k1, tau, alpha, state);

  // k2 = f(y(t)+(1/2)*k1 , t+(1/2)*tau) * tau.

  double y_2[N];

  for(i1 = 0; i1 < N; i1++)
    y_2[i1] = y_t[i1] + 0.5*k1[i1];

  double k2[N];
  f_times_tau(y_2, t + 0.5*tau, k2, tau, alpha, state);
  
  // k3

  double y_3[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_3[i1] = y_t[i1] + 0.5*k2[i1];

  double k3[N];
  f_times_tau(y_3, t + 0.5*tau, k3, tau, alpha, state);

  // k4

  double y_4[N];
 
  for(i1 = 0; i1 < N; i1++)
    y_4[i1] = y_t[i1] + k3[i1];

  double k4[N];
  f_times_tau(y_4, t + tau, k3, tau, alpha, state);

  // final

  for(i1 = 0; i1 < N; i1++){
       y_t_plus_tau[i1] = y_t[i1] + 1./6. * 
       (k1[i1] + 2.*k2[i1] + 2.*k3[i1] + k4[i1]);
  }
}


// Calculates f(y(t),t) * tau

void f_times_tau(double *y_t, double t, double *f_times_tau_, 
                 double tau, vector<double> *alpha, 
                 double(*state)(double, vector<double> *))
{
  if(N != 2)
    {
      fprintf(stderr, "Error: N != 2!\n");
      exit(EXIT_FAILURE);
    }

  // TOV equation implementation
  f_times_tau_[0] = tau * (  -(state( y_t[0], alpha) + y_t[0]) 
                    * (y_t[1] + 4*M_PI*pow(t, 3.)*y_t[0]) 
                    / (-2*y_t[1]*t + pow(t, 2.0)));
  f_times_tau_[1] = tau * 4*M_PI * pow(t, 2.) *state(y_t[0], alpha);
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

 
// LOWER BOUND FOR RK STEP [unfinished] 

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
*/
