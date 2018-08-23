#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<time.h>
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
  int    m_max   = 0;

  vector<double> alpha(4);
  vector<double> one_MR(2);
  vector<vector<double>> MR_rel;

  double p_init  = 0.0;
  double p_end   = 0.0;
  double p_dur   = 0.0;
  double slope   = 0.03;
  double slope_step = 0.01;

  bool   flag_s  = false;
  double diff0_s = 0.0;
  double diff_s  = 0.0;
  bool   flag    = false;
  double diff0   = 0.0;
  double diff    = 0.0;

  double alpha3_old;
  double e_rec;
  double pstep;

  vector<double> store;
  vector<double> pao_store;
  vector<vector<double>> reconstruction;
  bool   one     = true;
  int    n       = 0;

  // File I/O
 
  FILE *MRR = fopen("mr.out", "r");
 
  if (MRR == NULL) {
    exit(0);
  }
 
  while (1) {
    if (fscanf(MRR, "%lf,%lf", &R, &M) == EOF) {
      break;
    }
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


  // Initialization

  p_init = 0.00001 + mcount * 0.00001;
  e_rec = eos(p_init, &alpha);
  t[0] = 0.0000000001;
  p_dur = p_init;

  alpha[0] = p_init;
  alpha[1] = e_rec;
  alpha[2] = 5 * p_init;
  alpha[3] = e_rec + (alpha[2]-alpha[0]) / slope;

  while (mcount < 20) {
    flag_s = false;
    
    while (slope < 1.0) {
      alpha3_old = alpha[3];
      slope = (alpha[2]-alpha[0]) / (alpha[3]-e_rec); 
      alpha[3] = e_rec + (alpha[2]-alpha[0]) 
                       / (slope + slope_step); 
      //   Does this make sense?  ^^^^^^^^^^

      //slope = (alpha[2]-alpha[0]) / (alpha[3]-e_rec); 
     
      flag = false;
      pstep = 1e-6;
      p_end = p_dur - 5 * pstep;
                        // this destroys the three if statements!!!
      while (p_end >= 0.8 * p_dur) {
        p_end += pstep;
	y_0[0] = p_end; 
	y_0[1] = 0.0;
	// t[0] = 0.00000000001; // ?

        // Initializing
	for (i1 = 0; i1 < N; i1++) {
	  y[0][i1] = y_0[i1];
	}

        // Integration
        // cout << "doing Euler\n";

	for (i1 = 0; y[i1][0] > 0; i1++) {
 
          // Part I.1
          if (!one && y[i1][0] > pao_store[pao_store.size()-1]) {
            tov_euler(y[i1], t[i1], y_tau, tau, &alpha, line);
            
	    for (i2 = 0; i2 < N; i2++) {
	      y[i1+1][i2] = y_tau[i2];
 	    } 

	    t[i1+1] = t[i1] + tau;
          }
            
          // Part I.2
          else if (one && y[i1][0] > p_init) {
            tov_euler(y[i1], t[i1], y_tau, tau, &alpha, eos);
	   
            for (i2 = 0; i2 < N; i2++) {
	      y[i1+1][i2] = y_tau[i2];
 	    } 

	    t[i1+1] = t[i1] + tau;
          }

          // Part II
          if (!one && y[i1][0] <= pao_store[pao_store.size()-1] 
                   && y[i1][0] > p_init
                /* && reconstruction.size() > 0 */)  {
           
            if (y[i1][0] > pao_store[pao_store.size()-(n+2)]) {
              tov_euler(y[i1], t[i1], y_tau, tau, &reconstruction[n], line);
              //cout << "II  " << n << endl;
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

          // Part III
          if (y[i1][0] <= p_init) {
	    tov_euler(y[i1], t[i1], y_tau, tau, &alpha, eos);

	    for (i2 = 0; i2 < N; i2++) {
	      y[i1+1][i2] = y_tau[i2];
 	    } 

	    t[i1+1] = t[i1] + tau;
          }

	} // end of mass calculation

        n = 0;
        cout << "Slope after mass: " <<slope<< endl;   
        cout << " P_end for mass:  " <<p_end<<endl;
 
        if (!flag) {
          diff0 = y[i1-1][1] / 1.4766 - MR_rel[mcount][1];
          printf("diff0 mass %g %g\n", pstep, diff0);
          n = 0;
          flag = true;
          continue;
        }

        else {
          diff = y[i1-1][1] / 1.4766 - MR_rel[mcount][1];
          printf("diff mass %g %g %g\n", pstep, diff, MR_rel[mcount][1]);
          n = 0;
          if (diff * diff0 > 0) {
            // cout << diff * diff0 << endl;
            continue;
          }

          pstep /= 10.0;

          if (pstep < 1e-9) {
            cout << "pstep < x breaking condition" << endl;
            break;
          }

          p_end -= 10.0 * pstep;
          continue;
        }

      } // end p_end < p_dur loop

      if (!flag_s) {
        diff0_s = t[i1 - 1] - MR_rel[mcount][0];
        printf("diff0 radius %f %f %f %f\n", slope_step, diff0_s, t[i1 - 1],
               MR_rel[mcount][0]);
        n = 0;
        flag_s = true;
        continue;
      }

      diff_s = t[i1 - 1] - MR_rel[mcount][0];
      printf("diff radius %f %f %f %f\n", slope_step, diff_s, t[i1 - 1],
             MR_rel[mcount][0]);

      if (fabs(diff_s) < 0.005) {
        cout << "fabs(diff_s) break condition" << endl;
        n = 0;
        break;
      }

      if (diff0_s * diff_s > 0) {  
        // cout << "diff0_s * diff_s > 0" << endl;
        n = 0;
        continue;
      }

      slope_step /= 20;

      if (slope_step < 1e-9) {
        cout << "slope_step break condition" << endl;
        n = 0;
        break;
      }

      n = 0;

      alpha[3] = alpha3_old;

    } // end slope loop

  //cout << "Radius: " << MR_rel[mcount][0] << " , " << t[i1-1] << endl;
  //cout << "Mass:   " << MR_rel[mcount][1] << " , " 
  //     << y[i1-1][1]/1.4766 << endl;
  //cout << t[i1-1] << "," << y[i1-1][1]/1.4766 << endl;
  
  n = 0;

  store.push_back(y[i1-1][1]/1.4766);

  if (store.size() > 1 && store[store.size()-1] == store[store.size()-2]
      /* && Error check?? Same problem as with my code then */    ) {
    mcount++;
    reconstruction.push_back(alpha);
    pao_store.push_back(p_end);
    n = 0;
 
    cout << "output: " << t[i1-1] << "," << y[i1-1][1]/1.4766 << endl;
 
    if (!one) {
      e_rec = line(p_end, &alpha);
    }
    else { 
      e_rec = eos(p_end, &alpha);
    }
 
    p_dur = p_end;

    alpha[0] = p_end;
    alpha[1] = e_rec;
    alpha[2] = 5 * p_end;
    alpha[3] = e_rec + 4*p_end / 0.03;

    cout << alpha[0] << " " << alpha[1] << " " 
         << alpha[2] << " " << alpha[3] << endl;

    // slope = 0.03; 
    slope_step = 0.01;
    one = false;
   }

  } // end mcount loop


  for (i1 = 0; i1 < store.size(); i1++) {
    cout << store[i1] << endl;
  }
  
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


/*
  for (i1 = 0; i1 < MR_rel.size(); i1++) {
    if (MR_rel[i1][1] <= 1) {
      mcount++;
    }
    if (MR_rel[i1][1] - MR_rel[i1+1][1]) {
      m_max = i1;
      break;
    }
  }
*/


