#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
using namespace std;

const int N = 2;
const int num_steps = 100000;

double eos(double p);

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau,
               double(*state)(double));

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau,
                 double(*state)(double));
 
main(){

  int i1, i2;  
  double tau = 0.01;
  double y[num_steps+1][N];
  double y_0[N];
  double y_tau[N];
  double t[num_steps+1];
  double p0;
  int P;
  vector<double> R;
  vector<double> M;
  vector<double> e;
  vector<double> p;
  
  for (P = 0; P < 100; P++) {
    p0 = 0.00001 + P * 0.00001;
    y_0[0] = p0; 
    y_0[1] = 0.0;
    t[0] = 0.00000000001;
    
    for (i1 = 0; i1 < N; i1++) {
      y[0][i1] = y_0[i1];
    }
    
    for (i1 = 0; y[i1][0] > 0; i1++) {
      tov_euler(y[i1], t[i1], y_tau, tau, eos);
      
      for (i2 = 0; i2 < N; i2++) {
        y[i1+1][i2] = y_tau[i2];
      }
    
      t[i1+1] = t[i1] + tau;
    }
    
    M.push_back(y[i1-1][1]/1.4766);
    R.push_back(t[i1-1]);
    e.push_back(eos(p0)/(1.51483228e-4));
    p.push_back(p0);
  }
  
  for (i1 = 0; i1 < M.size(); i1++) {
    cout << R[i1] << "," << M[i1] << "," 
         << e[i1] << "," << p[i1] << endl;
  }

  return 0;
}

double eos(double p) {
  return pow(p/10., 3./5.);
}

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau,
               double(*state)(double)) {
  int i1;
  
  double k1[N];
  f_times_tau(y_t, t, k1, tau, state);

  for (i1 = 0; i1 < N; i1++) {
    y_t_plus_tau[i1] = y_t[i1] + k1[i1];
  }
}

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau,
                 double(*state)(double)) {
  f_times_tau[0] = (-(state(y_t[0]) + y_t[0])) *
                   (y_t[1] + 4 * M_PI * pow(t, 3.) * y_t[0]) /
                   (-2*y_t[1] * t + pow(t, 2.)) * tau;
  f_times_tau[1] = tau * 4 * M_PI * pow(t, 2.) * state(y_t[0]);
}

