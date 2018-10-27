#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
using namespace std;

const int N = 2;
const int num_steps = 100000;
double m = 30;

double eos(double p, vector<double> *alpha);

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, vector<double> *alpha);

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, vector<double> *alpha);
 
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
  vector<double> enden;
  vector<double> press;
  vector<double> alpha;
  double e, p;
  double conv = 7.55616208*(1e+5);

  FILE *input = fopen("output.dat", "r");
  if (input == NULL) exit(0);

  while (1) {
    if (fscanf(input, "%lf,%lf", &e, &p) == EOF)
      break;
    alpha.push_back(p);
    alpha.push_back(e);
  }

  fclose(input);
  input = NULL;

  for (P = 0; P < 1000; P++) {
    p0 = alpha[0] + 0.00001 + P * 0.00001;
    y_0[0] = p0; 
    y_0[1] = 0.0;
    t[0] = 0.00000000001;
    
    for (i1 = 0; i1 < N; i1++) {
      y[0][i1] = y_0[i1];
    }
    
    for (i1 = 0; y[i1][0] > 0; i1++) {
      tov_euler(y[i1], t[i1], y_tau, tau, &alpha);
      
      for (i2 = 0; i2 < N; i2++) {
        y[i1+1][i2] = y_tau[i2];
      }
    
      t[i1+1] = t[i1] + tau;
    }
    
    M.push_back(y[i1-1][1]/1.4766);
    R.push_back(t[i1-1]);
    enden.push_back(eos(p0, &alpha)*conv);
    press.push_back(p0*conv);
  }
  
  for (i1 = 0; i1 < M.size(); i1++) {
    printf("%5.10lf,%5.10lf,%5.10lf,%5.10lf\n", 
           R[i1], M[i1], enden[i1], press[i1]);
    //cout << R[i1] << "," << M[i1] << endl;
  }

  return 0;
}

double eos(double p, vector<double> *alpha) {
  int i;

  for (i = 2; i < alpha->size(); i += 2) {
    if (p < (*alpha)[i])
      break;
  }  

  double e = (p - (*alpha)[i-2]) * 
             ((*alpha)[i+1] - (*alpha)[i-1]) / ((*alpha)[i] - (*alpha)[i-2])
             + (*alpha)[i-1];
  return e;
}

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, vector<double> *alpha) {
  int i1;
  
  double k1[N];
  f_times_tau(y_t, t, k1, tau, alpha);

  for (i1 = 0; i1 < N; i1++) {
    y_t_plus_tau[i1] = y_t[i1] + k1[i1];
  }
}

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, vector<double> *alpha) {
  f_times_tau[0] = (-(eos(y_t[0], alpha) + y_t[0])) *
                   (y_t[1] + 4 * M_PI * pow(t, 3.) * y_t[0]) /
                   (-2*y_t[1] * t + pow(t, 2.)) * tau;
  f_times_tau[1] = tau * 4 * M_PI * pow(t, 2.) * eos(y_t[0], alpha);
}

