#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<fstream>
using namespace std;


// Global parameters

const int N = 2;
const int num_steps = 100000;


// Function heads

double eos(double p, vector<double> *alpha);

double line(double p, vector<double> *alpha);

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, vector<double> *alpha,
               double(*state)(double, vector<double> *));

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, vector<double> *alpha,
                 double(*state)(double, vector<double> *));


// Main function

main(){

  // Variables

  int    i1, i2;  
  double tau     = 0.01;
  double y[num_steps+1][N];
  double y_0[N];
  double y_tau[N];
  double t[num_steps+1];
  double p0      = 0.0;
  int    P       = 0;
  double M, R;
  int    mcount  = 0;

  vector<double> alpha(4);
  vector<double> one_MR(2);
  vector<vector<double>> MR_rel;

  double p_init  = 0.0;
  double p_end   = 0.0;
  double p_dur   = 0.0;
  double slope   = 0.07;
  double slope_step = 0.01;

  double alpha3_old;
  double e_rec;
  double pstep;

  double M_err   = 0.01;
  double R_err   = 0.0001;
  int    l       = -3;  

  vector<double> pao_store; 
  vector<vector<double>> reconstruction;
  bool   one     = true;
  int    n       = 0;
  int    m_init  = 0;
  double dial    = 2;


  // File I/O
 
  FILE *MRR = fopen("mr.out", "r");
  if (MRR == NULL)
    exit(0);
  while (1) {
    if (fscanf(MRR, "%lf,%lf", &R, &M) == EOF)
      break;
    one_MR[0] = R;
    one_MR[1] = M;
    MR_rel.push_back(one_MR);
  }
  fclose(MRR);
  MRR = NULL;

  for (i1 = 0; i1 < MR_rel.size(); i1++) {
    if (MR_rel[i1][1] >= 1)
      break;
    mcount++;
  }
  m_init = mcount;

  // Initialization

  p_init = 0.00001 + mcount * 0.00001;
  e_rec = eos(p_init, &alpha);
  t[0] = 0.0000000001;
  p_dur = p_init;
  p_end = 5 * p_dur;

  alpha[0] = p_init;
  alpha[1] = e_rec;
  alpha[2] = 5 * p_init;
  alpha[3] = e_rec; // + (alpha[2]-alpha[0]) / slope;

  while (mcount < 20) {
    // let's see

      while (p_end > p_dur) {
        p_end = p_dur + l*pow(10, -7);
	y_0[0] = p_end; 
	y_0[1] = 0.0;

        // Initializing
	for (i1 = 0; i1 < N; i1++) {
	  y[0][i1] = y_0[i1];
	}

        // Integration
        // cout << "doing Euler\n";

	for (i1 = 0; y[i1][0] > 0; i1++) {

          // PART I.1
          if (!one && y[i1][0] > pao_store[pao_store.size()-1]) {
            tov_euler(y[i1], t[i1], y_tau, tau, &alpha, line);
            // cout << " p I.1: " << y[i1][0] << " p_store: " 
            //     << pao_store[pao_store.size()-1] << endl;
	    for (i2 = 0; i2 < N; i2++) {
	      y[i1+1][i2] = y_tau[i2];
 	    } 

	    t[i1+1] = t[i1] + tau;
          }

          // PART I.2
          else if (one && y[i1][0] > p_init) {
	    tov_euler(y[i1], t[i1], y_tau, tau, &alpha, eos);
            // cout << "I.2 (one)\n";
	    for (i2 = 0; i2 < N; i2++) {
	      y[i1+1][i2] = y_tau[i2];
 	    } 

	    t[i1+1] = t[i1] + tau;
          }

          // PART II
          if (!one && y[i1][0] <= pao_store[pao_store.size()-1] 
                   && y[i1][0] > p_init) {

            if (y[i1][0] > pao_store[pao_store.size()-(n+2)]) {
	      tov_euler(y[i1], t[i1], y_tau, tau, &reconstruction[n], line);
              //cout << " p II : " << y[i1][0] << " p_store: " 
              //     << reconstruction[n][1] << endl;
	      for (i2 = 0; i2 < N; i2++) {
	        y[i1+1][i2] = y_tau[i2];
 	      } 

	      t[i1+1] = t[i1] + tau;
            }
   
            else if (n < pao_store.size()) {
              n++;
              i1--;
            }
          }
 
          // PART III
          if (y[i1][0] <= p_init) {
	    tov_euler(y[i1], t[i1], y_tau, tau, &alpha, eos);
            // cout << "III  " << y[i1][0]  << endl;
	    for (i2 = 0; i2 < N; i2++) {
	      y[i1+1][i2] = y_tau[i2];
 	    } 

	    t[i1+1] = t[i1] + tau;
          }

	} // end of mass calculation

        
     /* // Absolute error
        if (fabs(y[i1-1][1]/1.4766 - MR_rel[mcount][1]) < M_err
            && fabs(t[i1-1] - MR_rel[mcount][0]) < R_err) {
     */  
        // Relative error
        if (fabs(1 - MR_rel[mcount][1] / (y[i1-1][1]/1.4766)) < M_err
            && fabs(1 - MR_rel[mcount][0] / t[i1-1]) < R_err) {
 

          if (!one && (alpha[3]-alpha[1])/(alpha[2]-alpha[0]) > 1.) {

            cout << "part three (!one)\n";
            cout << "    Radius:  " << MR_rel[mcount][0] 
                 << " , " << t[i1-1] << endl;
            cout << "     Mass:   " << MR_rel[mcount][1] 
                 << " , " << y[i1-1][1]/1.4766 << endl;
            cout << " |R_c-R_f| = " 
                 << fabs(MR_rel[mcount][0] - t[i1-1]) << endl;
            cout << " |M_c-M_f| = " 
                 << fabs(MR_rel[mcount][1] - y[i1-1][1]/1.4766) << endl;
            cout << "   Slope   = " << (alpha[3]-alpha[1])/(alpha[2]-alpha[0]) 
                 << endl;
            
            reconstruction.push_back(alpha);
            pao_store.push_back(p_end);

            e_rec = line(p_end, &alpha);

            alpha[0] = p_end;
            alpha[1] = e_rec;
            alpha[2] = 5 * p_end;
            alpha[3] = e_rec + pow(10., -4.);
 
            mcount++;
            one = false;
            l = 1;
            p_end = 5 * p_dur;
            M_err += 0.02 * dial;
            R_err += 0.02 * dial;  
            cout << M_err << " , " << R_err << endl;
            dial++;
          }

          else if (one) {

            cout << "part three (one)\n";
            cout << "Radius: " << MR_rel[mcount][0] 
                 << " , " << t[i1-1] << endl;
            cout << "Mass:   " << MR_rel[mcount][1] 
                 << " , " << y[i1-1][1]/1.4766 << endl;
            cout << "|R_c-R_f| =  " 
                 << fabs(MR_rel[mcount][0] - t[i1-1]) << endl;
            cout << "|M_c-M_f| =  " 
                 << fabs(MR_rel[mcount][1] - y[i1-1][1]/1.4766) << endl;

            reconstruction.push_back(alpha);
            pao_store.push_back(p_end);

            e_rec = eos(p_end, &alpha);

            alpha[0] = p_end;
            alpha[1] = e_rec;
            alpha[2] = 5 * p_end;
            alpha[3] = e_rec + pow(10., -3.);
 
            one = false;
            l = 1;
            p_end = 5 * p_dur;
            //M_err += 0.01;
            //R_err += 0.01;
            //M_err = .5;
            
            mcount++;
            //cout << M_err << " , " << R_err << endl;    
          }
        }

        if (p_end <= 1.5*p_dur) {
         /* 
          if( mcount > m_init+1) {           
            cout << "else if\n";
            cout << "Radius: " << MR_rel[mcount][0] 
                 << " , " << t[i1-1] << endl;
            cout << "Mass:   " << MR_rel[mcount][1] 
                 << " , " << y[i1-1][1]/1.4766 << endl;
            cout << "   Slope   = " << (alpha[3]-alpha[1])/(alpha[2]-alpha[0]) 
                 << endl; 
          }
        */
          l++; 
          n = 0;
          p_end = 5 * p_dur;
        }
            
        else if (!one) {

          cout << "else\n";
          cout << "Radius: " << MR_rel[mcount][0] 
               << " , " << t[i1-1] << endl;
          cout << "Mass:   " << MR_rel[mcount][1] 
               << " , " << y[i1-1][1]/1.4766 << endl;
          cout << "   Slope   = " << (alpha[3]-alpha[1])/(alpha[2]-alpha[0])
               << endl;

          alpha[3] += pow(10, -6);

          l = -1; 
          p_end = 5 * p_dur;
          // i1 = 0;
        }

      } // end p_end < p_dur loop
/*
   cout << "Radius: " << MR_rel[mcount][0] 
                      << " , " << t[i1-1] << endl;
   cout << "Mass:   " << MR_rel[mcount][1] 
                      << " , " << y[i1-1][1]/1.4766 << endl;
*/
  } // end mcount loop
 
  // output...
 
  return 0;
}


// Functions

double eos(double p, vector<double> *alpha) {
  return pow(p/10., 3./5.);
}

double line(double p, vector<double> *alpha) {
   return p * ((*alpha)[3]-(*alpha)[1])/((*alpha)[2]-(*alpha)[0]) 
          + ((*alpha)[1] - (*alpha)[0] * 
          ((*alpha)[3]-(*alpha)[1])/((*alpha)[2]-(*alpha)[0]));
}

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, vector<double> *alpha,
               double(*state)(double, vector<double> *)) {
  int i1;
  
  double k1[N];
  f_times_tau(y_t, t, k1, tau, alpha, state);

  for (i1 = 0; i1 < N; i1++) {
    y_t_plus_tau[i1] = y_t[i1] + k1[i1];
  }
}

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, vector<double> *alpha,
                 double(*state)(double, vector<double> *)) {
  f_times_tau[0] = (-(state(y_t[0], alpha) + y_t[0])) *
                   (y_t[1] + 4 * M_PI * pow(t, 3.) * y_t[0]) /
                   (-2*y_t[1] * t + pow(t, 2.)) * tau;
  f_times_tau[1] = tau * 4 * M_PI * pow(t, 2.) * state(y_t[0], alpha);
}

