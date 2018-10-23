#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<time.h>
using namespace std;


// Global parameters

const int N = 2;
const int num_steps = 1000000;
double *y;
double p_end   = 0.0;
const double tau  = 0.0001;
double t[num_steps+1];
double p_dur   = 0.0;
double Rcomp,Mcomp;
int indx = 2;
double eps0,epsp,epsm,slope_old=0,lambda;
vector<double> alpha(4);
vector<vector<double>> MR_rel;
int    mcount  = 0;
double p_init;


vector<double> mc(400);


// Function heads

// line contains the eos (both parts connected at a cutoff pressure pcut
double line(double p, double pcut, vector<double> *alpha);

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, double p_cut, 
               vector<double> *alpha);

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, double p_cut, 
                 vector<double> *alpha);

// function calculates one tov star for given eos (alpha) and central pressure (y_0)
int tov(double* y_0,double p_cut,vector<double> *alpha);

// calculates the optimum star for a given slope
void getR();

void monte_carlo(int n);



// Main function


main(){

  // Variables

  int    i1, i2;  
  double y_0[N];
  double y_tau[N];
  double p0      = 0.0;
  int    P       = 0;
  double M, R;
  double pbest,ebest;

  vector<double> one_MR(3);

  y=(double*)malloc(N*num_steps*sizeof(double));

  p_init  = 0.0;
  double slope   = 1.00;
  double slope_step = -0.05;

  bool   flag_s  = false;
  double diff0_s = 0.0;
  double diff_s  = 0.0;
  bool   flag    = false;
  double diff0   = 0.0;
  double diff    = 0.0;
  bool   first   = true;
  double alpha3_old;
  double e_rec;
  double pstep;
  int    indx    = 2;
  int n = 0;

  // File I/O

  FILE *fres = fopen("results","w"); 
  FILE *MRR = fopen("mr.out", "r");
  double pcenter;
  if (MRR == NULL)
    exit(0);
  while (1) {
    if (fscanf(MRR, "%lf %lf %lf %*f", &R, &M, &pcenter) == EOF)
      break;
    one_MR[0] = R;
    one_MR[1] = M;
    one_MR[2] = pcenter;
    MR_rel.push_back(one_MR);
  }
  fclose(MRR);
  MRR = NULL;

/*
  // check for maxima in MRR

  vector<int> maxmin;
  double maxdiff0 = MR_rel[1][1] - MR_rel[0][1];
  double maxdiff  = 0;
  int    mmindx   = 0;
  for (i1 = 2; i1 < MR_rel.size(); i1++) {
    maxdiff = MR_rel[i1][1] - MR_rel[i1-1];
    if (maxdiff * maxdiff0 < 0) {
      max.push_back(i1);
      maxdiff0 = maxdiff;
    }
  }
*/

  // set initial mcount

  for (i1 = 0; i1 < MR_rel.size(); i1++) {
    if (MR_rel[i1][1] >= 1)
      break;
  }
  mcount=i1+1; 


  // Initialization

  p_init   = MR_rel[mcount-1][2];
  p_dur    = MR_rel[mcount-1][2];
  e_rec    = line(p_dur, p_dur, &alpha);
  t[0]     = 0.0000000001;
  first    = true;
  alpha[0] = p_init;     // this has to be rewritten using indx
  alpha[1] = e_rec;
  indx     = alpha.size()-2;
  alpha[2] = 5 * p_init;
  alpha[3] = e_rec + (alpha[2]-alpha[0]) / slope;
  indx     = 2;
  
  lambda = 0;


  double psave, esave;


  while (mcount < 40) {

    pstep = 1e-7;
    flag_s = false;
    slope_step = -0.05;
    double d0,dplus,dminus;

    int num = 0;

    if (!first) {  
      alpha.push_back(alpha[alpha.size()-2] + 1000*pstep);
      alpha.push_back(alpha[alpha.size()-1] + 1000*pstep / slope);
    }
    indx = alpha.size()-2;

    cout << "about to enter n loop\n";
    n = 200;

      // generate initial diff0_s

      monte_carlo(n);
      cout << mc.size() << endl;
      alpha[indx+1] = mc[1];

      cout << "alpha[indx-2] = " << alpha[indx-2] << endl;
      getR();  
      diff0_s = Rcomp - MR_rel[mcount][0];

      printf("Rcomp_bef: %lf; MR_rel: %f\n", Rcomp, MR_rel[mcount][0]);
      printf("Mcomp_bef: %lf; MR_rel: %f\n", Mcomp, MR_rel[mcount][1]);

      printf("diff0_s = %lf\n", diff0_s);


    while (n <= mc.size()) {
      // generate points
      monte_carlo(n);
      cout << mc.size() << endl;
      alpha[indx+1] = mc[1];
      cout << "indx = " << indx << endl;
      //cout << "alpha[indx-2] = " << alpha[indx-2] << endl;
      //printf("Rcomp : %lf; MR_rel: %f\n", Rcomp, MR_rel[mcount][0]);
      //printf("Mcomp : %lf; MR_rel: %f\n", Mcomp, MR_rel[mcount][1]);

      //printf("diff0_s = %lf\n", diff0_s);

      for (num = 2; num < n; num += 2) {
        // run through all monte-carlo generated e/p pairs
        alpha[indx+1] = mc[num + 1];
        alpha[indx] = mc[num];
        printf("a[i+1] : %lf\n", alpha[indx+1]);
        printf("a[i] : %lf\n", alpha[indx]);

        getR();
        diff_s = Rcomp - MR_rel[mcount][0];

        //printf("Rcomp : %lf; MR_rel: %f\n", Rcomp, MR_rel[mcount][0]);
        //printf("Mcomp : %lf; MR_rel: %f\n", Mcomp, MR_rel[mcount][1]);
        printf("diff0_s = %lf\n", diff0_s);
        printf("diff_s  = %lf\n", diff_s);

        if (fabs(diff_s) < fabs(diff0_s)) {
          // if diff is better than before, save it
          diff0_s = fabs(diff_s);
          psave = mc[num];
          esave = mc[num+1];
          cout << "esave if: " << esave << "   psave if: " << psave << endl;
       }
      } 

      if (fabs(diff0_s) < 0.00001) {
        break;
      }
    
      n *= 2;
    }

    Rcomp = diff0_s + MR_rel[mcount][0];
    printf("Rcomp_aft: %lf; MR_rel: %f\n", Rcomp, MR_rel[mcount][0]);
    printf("Mcomp_aft: %lf; MR_rel: %f\n", Mcomp, MR_rel[mcount][1]);

    // use best e & p for future calculation
    alpha[indx] = psave;
    alpha[indx+1] = esave;

    // cout << "esave: " << esave << "   psave: " << psave << endl;

    first = false;
/* 
  while (slope > 0.0) {
      slope += slope_step;

      if(slope<=0) {
        slope_step -= slope_step;
        slope_step /= 10;
        slope += slope_step;
      }

      alpha3_old = alpha[indx+1];
      alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/slope; 

      getR();

      if(p_end>alpha[indx-2]*5.0) {
        cout << p_end-alpha[indx-2]*5.0 << endl;
        continue;
      }
      if(fabs(diff)>1e-4) 
        continue;

      diff_s = Rcomp - MR_rel[mcount][0];
      epsp = diff_s*diff_s +lambda*(slope-slope_old)*(slope-slope_old);
      printf("diff radius %g %g %g %g %g %g %g %g \n",epsp,slope,slope_step, 
      diff_s,Rcomp,MR_rel[mcount][0],Mcomp,MR_rel[mcount][1]);

      if (fabs(diff_s) < 0.000001) {
        break;
      }
      if (epsp < eps0) {  
        eps0=epsp;
        continue;
      }
      slope_step /= -10;
      eps0=epsp;


    } // end slope loop
*/
    // lambda=1e-3;
    // slope_old=slope;
    // alpha[indx]=y[0];
    // alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2]) 
    //                  / slope ; 

    // slope = (alpha[indx+1] - alpha[indx-1] ) / (alpha[indx]-alpha[indx-2]); 

    indx=alpha.size()-2;

    fprintf(fres,"%e %e %e %e %e %e %e\n",MR_rel[mcount][0],Rcomp,
            MR_rel[mcount][1],Mcomp,p_end,line(p_end,p_dur,&alpha),slope);
 
    fflush(fres);

    cout << "Radius: " << MR_rel[mcount][0] << " , " << Rcomp << endl;
    cout << "Mass:   " << MR_rel[mcount][1] << " , " << Mcomp << endl;
    cout << "central pressure: " << p_end << endl;
    cout << "slope: " << slope << endl;
    cout << "mcount: " << mcount << endl;

    mcount++;
  } // end mcount loop
 

  free(y); 
  fclose(fres);
  return 0;
}


// Functions

int tov(double* y_0,double p_cut,vector<double> *alpha) {
  int i1,i2;
  double y_tau[N];

  // Initializing
  for (i1 = 0; i1 < N; i1++) {
    y[i1] = y_0[i1];
  }

  for (i1 = 0; y[i1*N] > 0; i1++) {
    tov_euler(&y[i1*N], t[i1], y_tau, tau, p_cut, alpha);

    for (i2 = 0; i2 < N; i2++) {
      y[(i1+1)*N+i2] = y_tau[i2];
    }

    t[i1+1] = t[i1] + tau;
  }

  return i1;
}

double line(double p,double p_cut, vector<double> *alpha) {
  int i;

  // comment this out for external eos input
  if(p <= p_cut) {
    // cout << "normal eos\n";
    return pow(p/10.,0.6);
  }

  for(i=2;i<alpha->size();i+=2) {
    if(p < (*alpha)[i]) {
      // cout << "p = " << p << "  alpha[i] = " << (*alpha)[i] << endl;
      break;
    }
  }

  double e =  (p - (*alpha)[i-2]) * 
         ((*alpha)[i+1]-(*alpha)[i-1])/((*alpha)[i]-(*alpha)[i-2]) 
          + (*alpha)[i-1] ;
  return e;
}

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, double p_cut,  
               vector<double> *alpha) {
  int i1;
  
  double k1[N];
  f_times_tau(y_t, t, k1, tau, p_cut,alpha);

  for (i1 = 0; i1 < N; i1++) {
    y_t_plus_tau[i1] = y_t[i1] + k1[i1];
  }
}

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, double p_cut, 
                 vector<double> *alpha) {
  f_times_tau[0] = (-(line(y_t[0],p_cut, alpha) + y_t[0])) *
                   (y_t[1] + 4 * M_PI * pow(t, 3.) * y_t[0]) /
                   (-2*y_t[1] * t + pow(t, 2.)) * tau;
  f_times_tau[1] = tau * 4 * M_PI * pow(t, 2.) * line(y_t[0],p_cut, alpha);
}

void getR() {
  double y_0[N];
  int i1;
//  bool flag = false;
//  double diff0,diff=0.0;
//  double pstep = 1e-7;
  p_end = alpha[indx];  // hmmm not sure about this
  //  while (p_end >= 0.8 * p_dur) {



  for (i1 = 0; i1 < alpha.size(); i1++) {
    cout << "alpha[i1] = " << alpha[i1] << endl;
  }

    cout << "p_end = " << p_end << endl;

    y_0[0] = p_end; 
    y_0[1] = 0.0;

    i1 = tov(y_0, p_dur, &alpha);   

    // end of mass calculation
    Rcomp = t[i1-1]+y[(i1-1)*N]*tau/(y[(i1-2)*N]-y[(i1-1)*N]);
    Mcomp = y[(i1-1)*N+1] +  
            (y[(i1-1)*N+1]-y[(i1-2)*N+1])/tau*(Rcomp-t[i1-1]);
    Mcomp /= 1.4766;

/*
    if (!flag) {
      diff0 = Mcomp - MR_rel[mcount][1];

      if(fabs(diff0)<1e-6) 
        break;

      if(diff0>0) {
         printf("masses seem to go down with increasing p_center.\n");
         exit(0);
      }
      flag = true;
      continue;
    }

    else {
      double diffx = Mcomp - MR_rel[mcount][1];
      if(diff>0.0) {
        if(fabs(diffx-diff)>0) {
          diff=diffx;
          break;
        }
      }
      diff = diffx;
      if(fabs(diff)<1e-6) 
        break;

      if (diff * diff0 > 0) {
          // cout << diff * diff0 << endl;
        continue;
      }

      pstep /= 10.0;
      diff = 0.0;
      if (pstep < 1e-9) {
          // cout << "pstep < x breaking condition" << endl;
        break;
      }
 
      p_end -= 10.0 * pstep;
      continue;
    }

//  } // end p_end < p_dur loop
*/
/*
  Rcomp = t[i1-1]+y[(i1-1)*N]*tau/(y[(i1-2)*N]-y[(i1-1)*N]);
  Mcomp = y[(i1-1)*N+1] + (y[(i1-1)*N+1]-y[(i1-2)*N+1])/tau*(Rcomp-t[i1-1]);
  Mcomp /= 1.4766;
*/
  return;
}

void monte_carlo(int n) {

  time_t t;
  struct tm tm;
  srand(time(NULL)); 

  cout << "mc.size() , n : " << mc.size() << " , "<< n << endl;

  for (int mc_c = 0; mc_c < n; mc_c += 2) {
    mc[mc_c] = (double)rand()/(double)RAND_MAX * 2e-5 + alpha[indx-2];
    mc[mc_c+1] = (double)rand()/(double)RAND_MAX * 0.06 + alpha[indx-1];
  }
  
  return;
}


