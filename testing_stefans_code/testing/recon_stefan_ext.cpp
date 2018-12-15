#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<fstream>
using namespace std;


// Global parameters

const int N = 2;
const int num_steps = 10000000;
double *y;
double p_end   = 0.0;
const double tau  = 0.01;
double t[num_steps+1];
double p_dur   = 0.0;
double Rcomp,Mcomp;
int indx = 2;
double eps0,epsp,epsm,slope_old=0,lambda;
vector<double> alpha;
vector<vector<double>> MR_rel;
int    mcount  = 0;
double p_init;
double x, e, p;

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

  vector<double> one_MR(2);

  y=(double*)malloc(N*num_steps*sizeof(double));

  p_init  = 0.0;
  double slope   = 1;
  double slope_step = -0.06;

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
//  int    indx    = 2;
  double convert = (1.3234)*(1e-6);

  // File I/O

  FILE *fres = fopen("results_el2_6_p","w"); 
  FILE *MRR  = fopen("mr_eos_alpha.out", "r");
  double pcenter;

  if (MRR == NULL)
    exit(0);

  while (1) {
    if (fscanf(MRR, "%lf,%lf,%lf,%lf", &R, &M, &e, &p) == EOF)
      break;
    one_MR[0] = R;
    one_MR[1] = M;
//    alpha.push_back(p*convert);
//    alpha.push_back(e*convert);
    MR_rel.push_back(one_MR);
  }

  fclose(MRR);
  MRR = NULL;
  cout << "read mrr\n";
  cout << "a.size(): " << alpha.size() << endl;

  cout << "MR_rel.size(): " << MR_rel.size() << endl;

  for (i1 = 0; i1 < MR_rel.size(); i1++) {
    if (MR_rel[i1][1] >= 1.7) {
//      p_init = alpha[2*i1];
//      p_dur = p_init;
      break;
    }
  }

  mcount = i1+1; 
  cout << "initial mcount: " << mcount << endl;

  // read external eos

  FILE *exteos = fopen("output.dat", "r");
  if (exteos == NULL)
    exit(0);
  while (1) {
    if (fscanf(exteos, "%lf,%lf", &e, &p) == EOF)
      break;
    alpha.push_back(p);
    alpha.push_back(e);
  }
  cout << "eos successfully read\n";
  fclose(exteos);
  exteos = NULL;

/* 
  // scan EoS

  for (i2 = 0; i2 < alpha.size()-2; i2+=2) {
    y_0[0] = alpha[indx-i2];
    y_0[1] = 0.0;
    tov (y_0, 0, &alpha);
    if (y[(i1-1)*N+1] < 1)
      break;
  }
*/

  for (int i3 = 2*i1+2+44; i3 < alpha.size(); i1++) {
    // for (int i3 = 2*i2+2; i3 < alpha.size(); i1++) {
    // removing all the parts from alpha that we do not need
    alpha.pop_back();
  }

  cout << "a.size(): " << alpha.size() << endl;

  //  p_init = alpha[2] /* + 0.00001 */ + mcount * 0.00001;
  p_init = alpha[alpha.size()-2];
  p_dur = p_init;

  indx = alpha.size() - 2;

  // adjusting alpha so that the last two fixed values are the central
  // values of the last star that we assume to know the EOS for
  alpha[indx] = p_init;
  alpha[indx+1] = line(p_dur, p_dur, &alpha);

  cout << "indx: " << indx << endl;
  alpha.push_back(5 * alpha[indx]);
  alpha.push_back(alpha[indx+1] + (alpha[indx+2]-alpha[indx]) / slope);
  indx += 2;
  cout << "pushed alpha\n";

  // Initialization

  //  p_init   = MR_rel[mcount-1][2];
  //  p_dur    = MR_rel[mcount-1][2];
  //  p_init = alpha[indx-2];
 
/*
  for (i2 = 0; i2 < alpha.size(); i2++) {
    cout << alpha[i2] << endl;
  }
*/
  e_rec    = line(p_dur, p_dur, &alpha);
  t[0]     = 0.0000000001;
  first    = true;
//  alpha[0] = p_init;     // this has to be rewritten using indx
//  alpha[1] = e_rec;
  indx     = alpha.size()-2;
//  alpha[2] = 5 * p_init;
//  alpha[3] = e_rec + (alpha[2]-alpha[0]) / slope;
//  indx     = 2;
  cout << "initialized\n";



// code optimizes a function eps = (R-R_wanted)**2 + lambda*(slope-slope_before)**2 do avoid big jumps in the slope
// as the slope_before is not known in the first round, use lambda=0 for optimization

  lambda = 0;

  while (mcount < 40) {

  /*
    if (mcount == maxmin[mmindx]) {
      // set lambda differently or destroy lambda completely
      // mmindx += 2; 
      // += 1 would mean the minimum... do we want that or the next max.?
      // only skip unstable branch?
    }
  */

    for (i2 = 0; i2 < alpha.size(); i2++) {
      cout << alpha[i2] << endl;
    }

    cout << "in mcount while loop\n";

    pstep  = 1e-7;
    flag_s = false;
    slope_step  =-0.06;
    double pold = alpha[alpha.size()-2];
    double eold = alpha[alpha.size()-1];
    double d0, dplus, dminus;


    if (!first) {  
      alpha.push_back(pold + 1000*pstep);
      alpha.push_back(eold + 1000*pstep / slope);
    }
    indx=alpha.size()-2;

// calculate eps for the initial slope (slope from previous mcount solution)  
// and a value below and above to determine whether slope_step should be positive or negative 

    cout << "before getR()\n";

    getR();
    diff0_s=Rcomp-MR_rel[mcount][0];
    eps0 = diff0_s*diff0_s +lambda*(slope-slope_old)*(slope-slope_old);
    printf("diff0 slope radius %f %f\n", slope, diff0_s);

    alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/(1.2*slope) ; 
    getR();
    dplus=Rcomp-MR_rel[mcount][0];
    epsp = dplus*dplus +lambda*(slope*1.2-slope_old)*(slope*1.2-slope_old);

    alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/(0.8*slope) ; 
    getR();
    dminus=Rcomp-MR_rel[mcount][0];
    epsm = dminus*dminus +lambda*(slope*0.8-slope_old)*(slope*0.8-slope_old);

    cout << "before ifs\n";

    if((epsm-eps0)*(epsp-eps0)<0) { 
      if(epsp<eps0) 
        slope_step = 0.06; 
      else 
        slope_step = -0.06; 
    }

    else {
      if(epsm>epsp) 
        slope_step = 0.06; 
      else 
        slope_step = -0.06; 
    }


   first = false;
// check wich direction to go
   
    while (slope > 0) {

      slope += slope_step;
      cout << "slope = " << slope << endl;
      cout << "slope_step = " << slope_step << endl;
     // cout << "in slope loop\n";

      if(slope < 0) {
        cout << slope << endl;
        cout << slope_step << endl;
        slope += fabs(slope_step);
        cout << slope << endl;
        slope_step /= 10;
        cout << slope_step << endl;
        slope += slope_step;
        cout << slope << endl;
      }

      alpha3_old = alpha[indx+1];
//      slope = (alpha[indx]-alpha[indx-2]) / (alpha[indx+1]-alpha[indx-1]); 
      alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/slope; 

      getR();
/*
      if(p_end>alpha[indx-2]*5.0) {
        cout << "wee : " <<  p_end-alpha[indx-2]*5.0 << endl;
        continue;
      }*/
      if(fabs(diff)>1e-4) 
        continue;

      diff_s = Rcomp - MR_rel[mcount][0];
      epsp = diff_s*diff_s +lambda*(slope-slope_old)*(slope-slope_old);
      printf("diff radius %g %g %g %g %g %g %g %g\n",epsp,slope,slope_step,
      diff_s, Rcomp,MR_rel[mcount][0],Mcomp,MR_rel[mcount][1]);
                      // 0.000001
      if (fabs(diff_s) < 0.0001) {
//      cout << "fabs(diff_s) break condition" << endl;
        break;
      }

      if (epsp < eps0) {  
      eps0=epsp;
      cout << "diff0_s * diff_s > 0" << endl;
        continue;
      }
// if eps is getting bigger, decrease slope step size and flip the sign
      slope_step /= -10;
      eps0=epsp;
/*
      if (fabs(slope_step) < 1e-7) {
        cout << "yo\n";
        break;
      }
*/

    } // end slope loop

// after first round set lambda to some reasonable value (one might play around with this number)
    lambda = (6)*(1e-2);
    slope_old = slope;
    alpha[indx] = y[0];
    alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2]) / slope ; 
//  alpha[indx+1]=line(y[0],p_dur,&alpha);


    indx=alpha.size()-2;

    fprintf(fres,"%e %e %e %e %e %e %e\n",MR_rel[mcount][0],Rcomp,MR_rel[mcount][1],Mcomp,p_end,line(p_end,p_dur,&alpha),slope);
 
    fflush(fres);

    cout << "Radius: " << MR_rel[mcount][0] << " , " << Rcomp << endl;
    cout << "Mass:   " << MR_rel[mcount][1] << " , " << Mcomp << endl;
    cout << "central pressure: " << p_end << endl;
    cout << "slope: " << slope << endl;
    cout << "mcount: " << mcount << endl;

    mcount++;
  } // end mcount loop
 
  // output...
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

  // Integration
  //  cout << "doing Euler\n";

  for (i1 = 0; y[i1*N] > 0; i1++) {
    
    //cout << "    p_cut: " << p_cut << endl;
    //cout << "a.size(): " << alpha->size() << endl;

    tov_euler(&y[i1*N], t[i1], y_tau, tau, p_cut, alpha);

    for (i2 = 0; i2 < N; i2++) {
      y[(i1+1)*N+i2] = y_tau[i2];
// //      cout << "y[(i1+1)*N+i2]: " << y[(i1+1)*N+i2] << endl; 
    }

/*
    Weird things are being calculated here: the radius increases
    too fast whereas the mass almost does not increase at all
*/

    t[i1+1] = t[i1] + tau;
    //cout << "t[i1+1]: " << t[i1+1] << endl;
  }
  
  return i1;
}

// end of mass calculation

double line(double p, double p_cut, vector<double> *alpha) {
  int i;

  // comment this out for external eos input
  // if(p <= p_cut) 
  //   return pow(p/10.,0.6);

  for(i=2;i < alpha->size(); i+=2) {
    if(p<(*alpha)[i]) break;
  }

  double e =  (*alpha)[i-1] + (p - (*alpha)[i-2]) * 
         ((*alpha)[i+1]-(*alpha)[i-1])/((*alpha)[i]-(*alpha)[i-2]);
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
  bool flag = false;
  double diff0,diff=0.0;
  double pstep = 1e-7;

/*
  for (int i2 = 0; i2 < alpha.size(); i2++) {
    cout << alpha[i2] << endl;
  }
*/

  indx = alpha.size()-2;

  p_end = alpha[indx-2]-pstep;
  
  //cout << "p_end = " << p_end << endl; 
  //cout << "p_dur = " << p_dur << endl; 

  while (p_end >= 0.8 * p_dur) {
 
   // cout << "In p while loop in getR()\n";
   // cout << "p_dur = " << p_dur << endl;
    if(p_end>alpha[indx-2]*5.0) 
      break;
    p_end += pstep;

    y_0[0] = p_end; 
    y_0[1] = 0.0;

    i1 = tov(y_0, p_dur, &alpha);    

    // cout << "mass calculated\n";

    Rcomp = t[i1-1]+y[(i1-1)*N]*tau/(y[(i1-2)*N]-y[(i1-1)*N]);
    Mcomp = y[(i1-1)*N+1] + 
            (y[(i1-1)*N+1]-y[(i1-2)*N+1])/tau*(Rcomp-t[i1-1]);
    Mcomp /= 1.4766;

    if (!flag) {
      diff0 = Mcomp - MR_rel[mcount][1];
      //cout << "diff0 set\n";
      if(fabs(diff0)<1e-6) 
        break;
/*
      if(diff0>0) {
         printf("masses seem to go down with increasing p_center.not yet implemented\n");
         exit(0);
      }*/
      flag = true;
      continue;
    }

    else {
      double diffx = Mcomp - MR_rel[mcount][1];
      //cout << "diffx set";
      if(diff>0.0) {
        if(fabs(diffx-diff)>0) {
          diff=diffx;
          break;
        }
      }
      diff = diffx;
      if(fabs(diff)<1e-6) 
        break;

  // printf("diff mass %g %g %g %g %g \n",pstep,
  // y[(i1-1)*N+1] / 1.4766, MR_rel[mcount][1],t[i1-1],MR_rel[mcount][0]);

      if (diff * diff0 > 0) {
          // cout << diff * diff0 << endl;
        continue;
      }

      pstep /= 10.0;
      diff = 0.0;
      if (pstep < 1e-9) {
       //pstep *= 10;
       // cout << "pstep < x breaking condition" << endl;
       break;
       //continue;
      }
 
      p_end -= 10.0 * pstep;
      continue;
    }

  } // end p_end < p_dur loop

  Rcomp = t[i1-1]+y[(i1-1)*N]*tau/(y[(i1-2)*N]-y[(i1-1)*N]);
  Mcomp = y[(i1-1)*N+1] + (y[(i1-1)*N+1]-y[(i1-2)*N+1])/tau*(Rcomp-t[i1-1]);
  Mcomp /= 1.4766;

//  cout << "getR() finished\n";

  return;
}

