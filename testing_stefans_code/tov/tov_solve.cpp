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
  vector<double> alpha/*(2)*/;
//  alpha[0]=0.0001; alpha[1]=0.0001;
  double e, p;
  double conv = 7.55616208*(1e+5);

  double d0, d1, d2, d3, d6;

  FILE *input = fopen("output.dat", "r");
  if (input == NULL) exit(0);
//  else cout << "input file open\n";

  while (1) {
//    if (fscanf(input, "%lf %lf %lf %lf %lf %lf %lf", 
//        &d0, &d1, &d2, &d3, &p, &e, &d6) == EOF)

    if (fscanf(input, "%lf,%lf", &e,&p) == EOF)
      break;
    alpha.push_back(p);
    alpha.push_back(e);
  }

//  for (i1 = 1; i1 < alpha.size(); i1+=2) {
//    cout << alpha[i1] << endl;
//  }

  fclose(input);
  input = NULL;

//  cout << "input read\n";

//  for (P = 0; P <= 1000000; P++) {
//    p0 = alpha[2] + P * 0.00001;

  for (P = 0; P <= alpha.size()-4; P += 4) {
    p0 = alpha[2+P];

    if (p0 > alpha[alpha.size()-2]) {
      break;
    }

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
/*
      if (y[i1][0] < alpha[0] && P == 500) {
        printf("%5.30lf,%5.20lf\n", y[i1][0], eos(y[i1][0], &alpha));
      }
*/
    }

//    cout << "pressure loop done\n";    
    
    M.push_back(y[i1-1][1]/1.4766);
    R.push_back(t[i1-1]);
    enden.push_back(eos(p0, &alpha)/* *conv */);
    press.push_back(p0/* *conv */);
  }

  for (i1 = 0; i1 < M.size(); i1++) {
    printf("%5.10lf,%5.10lf,%5.10lf,%5.10lf\n", 
           R[i1], M[i1], enden[i1], press[i1]);
  }

  return 0;
}

double eos(double p, vector<double> *alpha) {
  int i;
/*
  if (p < (*alpha)[0]) {
    return pow(p/10., 3./5.);
  }
*/
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

