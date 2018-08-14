#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <array>
#include <fstream>

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

const int num_steps_max = 100000000; // Max. RK step
// double tau = 0.01;               // Initial stepsize

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

  int i1 = 0;              // iteration
  int i2 = 0;              // iteration
  double P = 0.0;
  double p_rec = 0.0;
  double e_rec = 0.0;
  double M_err = 0.34;     // 0.006 solar mass tolerance
  double R_err = 3.01;       // 100m radius tolerance
  double p_step = 0.00001;  // arbitrary step
  int n = 0;
  double tau = 0.01;
  int n_max = 0;

  // *************************************

  // Memory allocation for t "radius" axis (discrete)
  // Also for the discrete ODE solution array

  vector<array<double, N>> EOS_arr;
  array<double, N> eos_tmp;
  
  vector<vector<double>> recon_storage;
  vector<double> alpha(4);
  vector<double> zero;

  double *t;
  one_alloc(num_steps_max, &t);

  double **y;
  two_alloc(num_steps_max, N, &y);

  vector<vector<double>> MR_rel;
  vector<double> one_MR(2);  

  // either l_st. or pao_st. will be used in the
  // program's final version (to be tested)

  vector<int> l_storage;
  vector<double> pao(2);
  vector<vector<double>> pao_storage;

  vector<vector<double>> EOS_out;
  vector<double> ep(2);
 
  // ****************************************

  // Read MR input file 

  double M = 0.0, R = 0.0;

  FILE *TOV = fopen("MRn.dat", "r");
  if (TOV == NULL) exit(0);
  int mcount = 0;

  while (1)  // at some point: rewrite this in C++
  {
     if(fscanf(TOV, "%lf,%lf", &R, &M) == EOF) break;
     one_MR[0] = R; one_MR[1] = M;
     MR_rel.push_back(one_MR);
  }

  fclose(TOV);
  TOV = NULL;

  for (i2 = 0; i2 <= MR_rel.size(); i2++) 
  {
     if (MR_rel[i2][1] >= 1) break;
     mcount++;
  }


  // ***************************************

  double p_init = 0.00001 + mcount * 0.00001;
  double y_0[N];
  double p_dur = p_init;
 

  // Initialize solution.

  t[0] = 0.000000001;

  p_rec = p_init + p_step;
  e_rec = eos(p_init, &zero);

  alpha[0] = p_init;   alpha[1] = eos(p_init, &zero);
  alpha[2] = 5*p_init; alpha[3] = e_rec;

  // ***************************************************


  // Do Runge-Kutta.

  double* y_tau;
  one_alloc(N, &y_tau);

  i1 = 0;
  
  int l = 1;
  double p_it = 0.0;
  double p_end = 5*p_dur;



        int crap = 0;

  double A, B;

  while(mcount < 15)
  {
     /*log*/// cout << mcount << "  " << "next round\n";
     /*log*/// cout << "  p_dur  = " << p_dur << endl;

     // Initialization

     // y_0[0] = p_rec; 
     // y_0[1] = 0.0;
  
     // for(i1 = 0; i1 < N; i1++)
     // {
     //    y[0][i1] = y_0[i1];
         /*log*/// cout <<  "y[0][i1] = " 
         /*log*///      << y[0][i1] << endl;     
     // }

  while( p_end > p_dur )
  {

      p_end  = p_dur + l * pow(10, -6);
      y_0[0] = p_end; 
      y_0[1] = 0.0;

      for(i1 = 0; i1 < N; i1++)
      {
         y[0][i1] = y_0[i1];
         /*log*/// cout <<  "y[0][i1] = " 
         /*log*///      << y[0][i1] << endl;     
      }
 

  for(i1 = 0; y[i1][0] > 0.0; i1++)
    {
 
      // *****************************************************
      // PART I.1: Reconstruction with an "interpolated" line

      if(y[i1][0] > p_dur)
       {  
	      // RK steps

              /*log*/// cout << "started 1.1\n";
              /*log*/// cout << " alpha[3] = " << alpha[3] << endl;
              /*log*/// cout << "[1.1] p - p_dur  = " << y[i1][0] - p_dur << endl;
 
	      RK_step(y[i1], t[i1], y_tau, tau, &alpha, line);

              /*log*/// cout << "|p - p_cal| = " << fabs(y_tau[0]-y[i1][0]) << endl;
 
	      if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -4.))	
	      {
		  // Accepting the step

		  for(i2 = 0; i2 < N; i2++)
		  { 
                     y[i1+1][i2] = y_tau[i2];
                     /*log*/// cout <<  "y[i1+1][i2] = " 
                     /*log*///      << y[i1+1][i2] << endl;     
                  }

		  t[i1+1] = t[i1] + tau;
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

       }  // End of I.1








      // ****************************************************
      // PART I.2: Have fun with previous lines 


      /*log*/// cout << " y[i1-1][0] = " << y[i1-1][0] << endl;
      /*log*/// cout << " y[i1][0]   = " << y[i1][0] << endl;
      /*log*/// cout << "    p_dur   = " << p_dur << endl;
      /*log*/// cout << "   p_init   = " << p_init << endl;
      /*log*/// cout << "p_dur - y[i1][0]  = " << p_dur - y[i1][0] << endl; 
      /*log*/// cout << "y[i1][0] - p_init = " << y[i1][0] - p_init << endl;



      if(    crap > 0 &&     y[i1][0] <= p_dur && y[i1][0] > p_init) 
      {
        
        /*log*/// cout << "[1.2] p - p_dur  = " << y[i1][0] - p_dur << endl;
        /*log*/// cout << "started part 1.2" << endl;
        /*log*/// cout << "             n  = " << n << endl;
        /*log*/// cout << "y[i1][0] - (p_dur - p_step * n) = " 
        /*log*///      <<  y[i1][0] - (p_dur - p_step * n)  << endl;

        if(y[i1][0] > pao_storage[n][0])
         {    
             // RK steps
	     
             /*log*/// cout << " p - p_dur  = " 
             /*log*///      << y[i1][0] - p_dur << endl;

  	     RK_step(y[i1], t[i1], y_tau, tau, 
                     &recon_storage[n], line);
            
             /*log*/// cout << "|p - p_cal| = " 
             /*log*///      << fabs(y_tau[0]-y[i1][0]) << endl;
 
             // Again, upper bound stepsize check
 
              if( fabs(y_tau[0]-y[i1][0]) <= pow(10., +4.))	
	      {
		  // Accepting the step

		  for(i2 = 0; i2 < N; i2++)
		  {   
                     y[i1+1][i2] = y_tau[i2];
                     /*log*/// cout <<  "y[i1+1][i2] = " 
                     /*log*///      << y[i1+1][i2] << endl;     
                  }

		  t[i1+1] = t[i1] + tau;

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

         }
       
         else if (n < n_max - 1)
            {
               n++; i1--;
            }
         
       //  else continue; 

      }   // End of Part I.2






       
      // *************************************************
      // PART II: Calculate the rest with known eos
 
      if(y[i1][0] <= p_init)
         {
 	      // RK steps
              
              /*log*/// cout << "2 started\n";

              /*log*/// cout << "[ 2 ] p - p_init = " << y[i1][0] - p_init << endl;
 
              /*log*///   cout << "[ 2 ]         p  = " << y[i1][0] << endl;
 

              RK_step(y[i1], t[i1], y_tau, tau, &zero, eos);
  
	      if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -4.))	
	      {
		  // Accepting the step

		  for(i2 = 0; i2 < N; i2++)
		  {
                     y[i1+1][i2] = y_tau[i2];
                     /*log*/// cout <<  "y[i1+1][i2] = " 
                     /*log*///      << y[i1+1][i2] << endl;  
		  }

                  t[i1+1] = t[i1] + tau;

                  /*log*/// cout << "part 2 accepted" << endl;     
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
	      } //
   
         }

         // **************************************

     }   // Mass calculation loop ends here





         // **********************************************
         // PART III: Checking if calculated mass is ok

	   if (fabs(y[i1-1][1]/1.4766 - MR_rel[mcount][1]) < M_err  // mass check
               && fabs(t[i1-1]-MR_rel[mcount][0]) < R_err           // radius check
               && (alpha[3]-alpha[1])/(alpha[2]-alpha[0]) > 1.01)   // slope check
           {

              recon_storage.push_back(alpha);      // storing line parameters
              l_storage.push_back(l);              // storing l for p_end

              pao[0] = p_dur; pao[1] = p_end;      // storing p_alpha and p_omega 
              pao_storage.push_back(pao);          // of the line needed for the mass
              
              /*
	      for(i2 = 0; i2 < recon_storage.size(); i2++)
	      {
                // some debug output, may be deleted later
	        cout << recon_storage[i2][0] << " " << recon_storage[i2][1] << " "
		     << recon_storage[i2][2] << " " << recon_storage[i2][3] << endl;
	      }
              */

              /*log*/ cout << "part three  \n";              
              /*log*/// cout << "n = " << n << ", " << "n_m = " << n_max << endl;
              /*log*/ cout << "|M - M_dat| = " 
              /*log*/      << fabs(y[i1-1][1]/1.4766 - MR_rel[mcount][1]) << endl;
              /*log*/ cout << "       Mass = " << y[i1-1][1] << endl;
              /*log*/ cout << " Solar Mass = " << y[i1-1][1]/1.4766 << endl;
              /*log*/ cout << " Input Mass = " << MR_rel[mcount][1] << endl;
  
              /*log*/ cout << "|R - R_dat| = " 
              /*log*/      << fabs(t[i1-1] - MR_rel[mcount][0]) << endl;  
              /*log*/ cout << "   Slope =    "<<(alpha[3]-alpha[1])/(alpha[2]-alpha[0]) << endl;
             

	      A = p_end;                   
	      B = line(p_end, &alpha);
	      cout << " P = " << A << endl;
              cout << " E = " << B << endl;
              //EOS_out.push_back(ep);
	      //cout << EOS_out[i1][0] << "," << EOS_out[i1][1] << endl;

              p_dur = p_end;
              alpha[0] = p_end; alpha[1] = alpha[3];
              p_rec += p_step;
              alpha[2] = p_end + p_step;   // modified

              mcount++;
              p_dur += p_step;
              n = 0;
              n_max++;
              
                               crap++;
              // continue;

              p_end = 5 * p_dur;
              /*log*/// cout << "n = " << n << ", " << "n_m = " << n_max << endl;
      
              tau = 0.01;
           }

           else if (p_end <= 5*p_dur) // cutoff for p_end is 5*p_d0r
	   {
              /*log*/// cout << "part three (else if)\n";
              /*log*/// cout << " y[i1-1][1] = " << y[i1-1][1]/1.4766 << endl;
              /*log*/// cout << " y[i1-1][0] = " << y[i1-1][0] << endl;
              /*log*/// cout << "     i1 =     " << i1 << endl;
              /*log*/// cout << "      l =     " << l << endl;
              /*log*/// cout << " y[i1+1][1] = " << y[i1+1][1]/1.4766 << endl;
              /*log*/// cout << " y[i1][1]   = " << y[i1][1]/1.4766 << endl;
              /*log*/// cout << " MR[mcount] = " << MR_rel[mcount][1] << endl;
              /*log*/// cout << "|M - M_dat| = " 
              /*log*///      << fabs(y[i1-1][1]/1.4766 - MR_rel[mcount][1]) << endl;  
              /*log*/// cout << "|R - R_dat| = " 
              /*log*///      << fabs(t[i1-1] - MR_rel[mcount][0]) << endl;  
              /*log*/// cout << " alpha[3]_b = "<< alpha[3] << endl;              
              /*log*/// cout << " alpha[3]_l = "<< alpha[3] << endl;    
              /*log*/// cout << "n = " << n << ", " << "n_m = " << n_max << endl;
	      /*log*/// cout << " elif Slope = "<<(alpha[3]-alpha[1])/(alpha[2]-alpha[0]) << endl;

        
              l++;
              i1 = 0;
              n = 0;
           }
 
           else 
	   {
              /*log*/// cout << "part three (else)\n";      
              /*log*/// cout << "     i1 =     " << i1 << endl;
              /*log*/// cout << "      l =     " << l << endl;
              /*log*/// cout << "|M - M_dat| = " 
              /*log*///      << fabs(y[i1-1][1]/1.4766 - MR_rel[mcount][1]) << endl;  
              /*log*/// cout << "|R - R_dat| = " 
              /*log*///     << fabs(t[i1-1] - MR_rel[mcount][0]) << endl; 
              /*log*/// cout << " alpha[3]_b = " << alpha[3] << endl;              
              /*log*/   cout << " else Slope = "<<(alpha[3]-alpha[1])/(alpha[2]-alpha[0]) << endl;

              alpha[3] += pow(10, -6);
              l = 1;
              p_end = 5* p_dur; 

              /*log*/// cout << " alpha[3]_l = "<< alpha[3] << endl;    

              i1 = 0;
              n = 0;
             
              /*log*/// cout << "n = " << n << ", " << "n_m = " << n_max << endl;

	   }  
         // **********************************************
    }

  } // while loop for mass <= <value> ends here


  one_free(N, &y_tau);    
 
  one_free(num_steps_max, &t);

  cout << endl; // WHAT IS THIS FOR??

  two_free(num_steps_max, N , &y);



  // *************************************************

  // OUTPUT FILE 
  // a fun exercise in C++

  cout << "woop reached output section\n";

//  vector<vector<double>> EOS_out;
//  vector<double> ep;
  ep.push_back(0); ep.push_back(0);

  double DPx = p_init/100;

/*
  for (i1 = 0; i1 < 100 ; i1++)
  {
      cout << " if " << endl;
      ep[0] = DPx * i1;
      ep[1] = eos(ep[0], &zero);
      EOS_out.push_back(ep);
  }
*/

/*
  for (i1 = 0; i1 < recon_storage.size(); i1++)
  {  
      ep[0] = p_init + l_storage[i1] * p_step; 
      ep[1] = line(ep[0], &recon_storage[i1]);
      EOS_out.push_back(ep);
      cout << EOS_out[i1][0] << "," << EOS_out[i1][1] << endl;
  }
*/

  for (i1 = 0; i1 < recon_storage.size(); i1++)
  {
      ep[0] =  pao_storage[i1][1];
      ep[1] = line(ep[0], &recon_storage[i1]);
      EOS_out.push_back(ep);
      cout << EOS_out[i1][0] << "," << EOS_out[i1][1] << endl;
  }





/*
  ofstream output ("output.dat");
   
  for(i2 = 0; i2 < EOS_out.size(); i2++)
  {
     output << EOS_out[i2][0] << "," << EOS_out[i2][1] << endl;
  } 

  output.close();
*/

  // *************************************************

  return EXIT_SUCCESS;

}

// *****************************************************


// *************
// * Functions *
// *************


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
/*
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
  f_times_tau(y_4, t + tau, k4, tau, alpha, state);
*/
  // final

  for(i1 = 0; i1 < N; i1++){
    y_t_plus_tau[i1] = y_t[i1] + k1[i1]; 
    // y_t_plus_tau[i1] = y_t[i1] + 1./6. * 
    // (k1[i1] + 2.*k2[i1] + 2.*k3[i1] + k4[i1]);
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



/*


      if(y[i1][0] <= p_dur && y[i1][0] > p_init) 
      {
        
        cout << "started part 1.2" << endl;
        cout << "             n  = " << n << endl;
        cout << "y[i1][0] - (p_dur - p_step * n) = " 
             <<  y[i1][0] - (p_dur - p_step * n)  << endl;

        if(y[i1][0] > p_dur - p_step * n)
         {    
             // RK steps
	     
             // cout << " p - p_dur  = " 
             //      << y[i1][0] - p_dur << endl;

  	     RK_step(y[i1], t[i1], y_tau, tau, 
                     &recon_storage[n-1], line);
            
             // cout << "|p - p_cal| = " 
             //      << fabs(y_tau[0]-y[i1][0]) << endl;
 
             // Again, upper bound stepsize check
 
              if( fabs(y_tau[0]-y[i1][0]) <= pow(10., -4.))	
	      {
		  // Accepting the step

		  for(i2 = 0; i2 < N; i2++)
		  {   
                     y[i1+1][i2] = y_tau[i2];
                      cout <<  "y[i1+1][i2] = " 
                           << y[i1+1][i2] << endl;     
                  }

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

         }
       
         else
         {
            if (n < n_max)
            {
               n++; i1--;
            }
            else continue;
         }

      }   // End of Part I.2





*/



/*
  cout << "init eos loop " << endl; 


  for(i2 = 0; ep[0] < p_init ; i2++)
  {
     ep[0] = DPx * i2;
     ep[1] = eos(DPx*i2, &zero);
     cout << ep[0] << " " << ep[1] << endl;
     EOS_out.push_back(ep);
     cout << EOS_out[i2][0] << " " <<  EOS_out[i2][1] << endl;

  }
 
  int a = i2;

  cout << " storage loop " << endl;
 
  for(i1 = 0; i1 < recon_storage.size(); i1++)
  {
     for(i2 = a+100; ep[0] < p_init + i1 * p_step; i2++)
     {
        ep[0] = DPx*i2;
        ep[1] = line(DPx*i2, &recon_storage[i2]);
        EOS_out.push_back(ep);
        cout << EOS_out[i2][0] << " " <<  EOS_out[i2][1] << endl;
     }

  } 

*/


              /*
	      if (y[i1-1][1]/1.4766 > MR_rel[mcount][1])
	      {
	         alpha[3] *= 1.2;
	      }
	      if (y[i1-1][1]/1.4766 < MR_rel[mcount][1])
	      {
		 alpha[3] /= 1.2;
	      }
              */
              
