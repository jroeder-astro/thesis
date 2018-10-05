#include<iostream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<fstream>
using namespace std;


// Global parameters

const int N         = 2;
const int num_steps = 10000000;
double    **y;
double    *t;
double    p_end     = 0.0;
const double tau    = 0.0001;
//double    t[num_steps+1];
double    p_dur     = 0.0;
double    Rcomp, Mcomp;  // what exactly are these
int       indx      = 2;
double    eps0, epsp, epsm;
double    slope_old = 0, lambda;
vector<double>         alpha(4);
vector<vector<double>> MR_rel;
int       mcount    = 0;
int what;

// Function heads

// line contains the eos (both parts connected at a cutoff pressure pcut)
double line(double p, double pcut, vector<double> *alpha);

void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, double p_cut, 
               vector<double> *alpha);

void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, double p_cut, 
                 vector<double> *alpha);

// function calculates one tov star for given eos (alpha) 
// and central pressure (y_0)
int tov(double* y_0, double p_cut,vector<double> *alpha);

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
  double pbest, ebest;

  vector<double> one_MR(3);
  vector<double> result_one(7);
  vector<vector<double>> result;

  // moved a few variables to the global section so that 
  // I can use them in functions without passing them on (just lazy)

  // y no longer 2D, "superindex"?  
  // other memory area...

  y = (double**)malloc((num_steps+1) * sizeof(double*));
  for(i1 = 0; i1 <= num_steps+1; i1++) {
    y[i1] = (double*)malloc(N*sizeof(double));
  }

  t = (double*)malloc((num_steps+1)*sizeof(double));

  double p_init  = 0.0;
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
  double pold, eold;
  double d0, dplus, dminus;


  // File I/O

//  FILE *fres = fopen("results_sep18.out","w"); 


  FILE *MRR = fopen("mr_sep18.out", "r");
  double pcenter;
  if (MRR == NULL)
    exit(0);

  while (1) {
    // mr.out also contains the corresponding pressures now
    if (fscanf(MRR, "%lf %lf %lf %*f", &R, &M, &pcenter) == EOF)
      break;

    one_MR[0] = R;
    one_MR[1] = M;
    one_MR[2] = pcenter;
    MR_rel.push_back(one_MR);
  }

  fclose(MRR);
  MRR = NULL;

  for (i1 = 0; i1 < MR_rel.size(); i1++) {
    if (MR_rel[i1][1] >= 1)
      break;
  }
  mcount   = i1 + 1; 


  // Initialization

  p_init   = MR_rel[mcount-1][2];
  p_dur    = MR_rel[mcount-1][2];
  e_rec    = line(p_dur,p_dur, &alpha);
  t[0]     = 0.0000000001;
  first    = true;
  alpha[0] = p_init;
  alpha[1] = e_rec;
  indx     = alpha.size()-2;
  alpha[2] = 5 * p_init;
  alpha[3] = e_rec + (alpha[2]-alpha[0]) / slope;
  indx     = 2; // only for curent alpha

  // code optimizes a function 
  // eps = (R-R_wanted)**2 + lambda*(slope-slope_before)**2 
  // to avoid big jumps in the slope; as the slope_before is not known 
  // in the first round, use lambda=0 for optimization
 
  lambda   = 0;

  cout << "about to start mcount loop\n";

  while (mcount < 50) {
    pstep  = 1e-7;
    flag_s = false;
    slope_step  = -0.05;
  
    cout << "alpha.size() = " << alpha.size() << endl;
    cout << "alpha[alpha.size()-2] = " << alpha[alpha.size()-2] << endl;
    cout << "alpha[alpha.size()-1] = " << alpha[alpha.size()-1] << endl;

    pold = alpha[alpha.size()-2];
    eold = alpha[alpha.size()-1];
  
    cout << "pold, eold set\n";

    // alpha contains more parameters now
    if (!first) {  
      alpha.push_back(pold + 1000 * pstep);
      alpha.push_back(eold + 1000 * pstep / slope);
    }
  
    indx = alpha.size()-2;

    // *********************************************************************
    // This happens only once per mass chosen from the input file

    // calculate eps for the initial slope (slope from previous mcount 
    // solution) and a value below and above to determine whether slope_step 
    // should be positive or negative 
    // eps is a variable not a function! 

    getR();
    diff0_s = Rcomp - MR_rel[mcount][0];            // vvvvv
    eps0 = diff0_s*diff0_s + lambda*(slope-slope_old)*(slope-slope_old);
    printf("diff0 slope radius %f %f\n", slope, diff0_s);

    cout << "1st getR() done\n";

    alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/(1.2*slope); 
    getR();                                         //           ^^^
    dplus = Rcomp - MR_rel[mcount][0];              // vvvvvvvvv
    epsp = dplus*dplus + lambda*(slope*1.2-slope_old)*(slope*1.2-slope_old);

    cout << "2nd getR() done\n";

    alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/(0.8*slope); 
    getR();                                         //           ^^^
    dminus = Rcomp - MR_rel[mcount][0];             //   vvvvvvvvv
    epsm = dminus*dminus + lambda*(slope*0.8-slope_old)*(slope*0.8-slope_old);

    cout << "3rd getR() done\n";

    if ((epsm - eps0) * (epsp - eps0) < 0) { 
      if (epsp < eps0) 
        slope_step = 0.05; 
      else 
        slope_step = -0.05; 
    }

    else {
      if (epsm > epsp) 
        slope_step = 0.05; 
      else 
        slope_step = -0.05; 
    }

    cout << "slope_step set\n";

    first = false;
    
    // *********************************************************************
   
    while (slope > 0.0) {
      slope += slope_step;

      if (slope <= 0) {
        // break;
        cout << "slope_step sign flip bc. slope <= 0\n";
        slope_step /= -10; 
        slope = slope_old;
      }

      alpha3_old = alpha[indx+1];
      // slope = (alpha[indx]-alpha[indx-2]) / (alpha[indx+1]-alpha[indx-1]); 
      alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2])/slope; 

      cout << "alphas set\n";

      // actual calculation
      getR();
      cout << "getR() done\n";

      if (p_end > alpha[indx-2]*5.0) {
        cout << "p_end > alpha[indx-2]*5.0\n";
        continue;
      }

      if (fabs(diff) > 1e-4) {
        cout << "fabs(diff) > 1e-4\n";
        continue;
      }

      diff_s = Rcomp - MR_rel[mcount][0];
      epsp   = diff_s*diff_s + lambda*(slope-slope_old)*(slope-slope_old);
      printf("diff radius %g %g %g %g %g %g %g %g \n",
             epsp, slope, slope_step, diff_s, Rcomp, MR_rel[mcount][0],
             Mcomp, MR_rel[mcount][1]);

      cout << "diff_s set, epsp re-evaluated\n";

      if (fabs(diff_s) < 0.00001) {
        // cout << "fabs(diff_s) break condition" << endl;
        break;
      }

      if (epsp < eps0) {  
        eps0 = epsp;
        // cout << "diff0_s * diff_s > 0" << endl;
        continue;
      }

      // if eps is getting bigger, decrease slope step size and flip the sign
      slope_step /= -10;
      eps0 = epsp;

      if (fabs(slope_step) < 1e-6) {
        // cout << "slope_step break condition" << endl;
        break;
      }

    } // end slope loop

    // after first round set lambda to some reasonable value 
    // (one might play around with this number)

    lambda = 1e-2;
    slope_old = slope;
    alpha[indx] = p_end;
    alpha[indx+1] = alpha[indx-1] + (alpha[indx]-alpha[indx-2]) / slope ; 
    indx = alpha.size()-2;

    result_one[0] = MR_rel[mcount][0];
    result_one[1] = Rcomp;
    result_one[2] = MR_rel[mcount][1];
    result_one[3] = Mcomp;
    result_one[4] = p_end;
    result_one[5] = line(p_end, p_dur, &alpha);
    result_one[6] = slope;

    result.push_back(result_one);

//  fprintf(fres,"%e,%e,%e,%e,%e,%e,%e\n", MR_rel[mcount][0], Rcomp, 
//          MR_rel[mcount][1], Mcomp, p_end, line(p_end,p_dur,&alpha), slope);
//  fflush(fres);

    cout << "Radius: " << MR_rel[mcount][0] << " , " << Rcomp << endl;
    cout << "Mass:   " << MR_rel[mcount][1] << " , " << Mcomp << endl;
    cout << "central pressure: " <<  p_end  << endl;
    cout << "slope: "  << slope  << endl;
    cout << "i1 = " << what << endl;
    
    cout << "mcount = " << mcount << endl;
    mcount++;
    cout << "mcount increased. mcount = " << mcount << endl;

  } // end mcount loop


  FILE *fres = fopen("results_sep18.out","w"); 

  for (i2 = 0; i2 < result.size(); i2++) {
    fprintf(fres, "%e,%e,%e,%e,%e,%e,%e\n", result[i2][0], result[i2][1], 
            result[i2][2], result[i2][3], result[i2][4], 
            result[i2][5], result[i2][6]);
  }
  
  fclose(fres);

//  fprintf(fres,"%e,%e,%e,%e,%e,%e,%e\n", MR_rel[mcount][0], Rcomp, 
//          MR_rel[mcount][1], Mcomp, p_end, line(p_end,p_dur,&alpha), slope);

//  fflush(fres);







  for (i2 = 0; i2 <= num_steps+1; i2++) 
    free(y[i2]);
  free(y); y = NULL;
  free(t); t = NULL;
  return 0;
}


// Functions


int tov(double* y_0, double p_cut, vector<double> *alpha) {
  int i1, i2;
  double y_tau[N];

  // initializing
  for (i1 = 0; i1 < N; i1++) {
    y[0][i1] = y_0[i1];
  }

  // integration
  for (i1 = 0; y[i1][0] > 0; i1++) {
    tov_euler(y[i1], t[i1], y_tau, tau, p_cut, alpha);

    for (i2 = 0; i2 < N; i2++) {
      y[i1+1][i2] = y_tau[i2];
    }

    t[i1+1] = t[i1] + tau;
  }

  // star is complete at i1

  return i1;
}


double line(double p, double p_cut, vector<double> *alpha) {
  int i;

  // initial eos
  if (p <= p_cut) 
    return pow(p/10., 0.6);

  // work with previously determined lines
  // only interval to be found
  // for (i = 2; i < alpha->size(); i += 2) {
  for (i = alpha->size(); i >= 2; i -= 2) {
    if (p < (*alpha)[i]) 
      break;
  }

  // line but with the "new" alpha
  double e = (p - (*alpha)[i-2]) * 
             ((*alpha)[i+1]-(*alpha)[i-1])/((*alpha)[i]-(*alpha)[i-2]) + 
             (*alpha)[i-1];
  return e;
}


void tov_euler(double *y_t, double t, 
               double *y_t_plus_tau, double tau, double p_cut,  
               vector<double> *alpha) {
  int i1;
  
  double k1[N];
  f_times_tau(y_t, t, k1, tau, p_cut, alpha);

  for (i1 = 0; i1 < N; i1++) {
    y_t_plus_tau[i1] = y_t[i1] + k1[i1];
  }
}


void f_times_tau(double *y_t, double t, 
                 double *f_times_tau, double tau, double p_cut, 
                 vector<double> *alpha) {
  f_times_tau[0] = (-(line(y_t[0], p_cut, alpha) + y_t[0])) *
                   (y_t[1] + 4 * M_PI * pow(t, 3.) * y_t[0]) /
                   (-2*y_t[1] * t + pow(t, 2.)) * tau;
  f_times_tau[1] = tau * 4 * M_PI * pow(t, 2.) * line(y_t[0], p_cut, alpha);   
}


void getR() {
  double y_0[N];
  int i1;
  bool flag = false;
  double diffx, diff0, diff = 0.0;
  double pstep = 1e-7;

  p_end = alpha[indx-2]-pstep;
 
  while (p_end >= 0.8 * p_dur) {
    // cutoff
    if (p_end > alpha[indx-2]*5.0) 
      break;
   
    p_end += pstep;
    y_0[0] = p_end; 
    y_0[1] = 0.0;
            
    // i1 for which star is complete
    i1 = tov(y_0, p_dur, &alpha);
    
    // why the mash-up of t and y arrays / masses and radii?
    // interpolation to hit the actual radius
    Rcomp = t[i1-1] + y[i1-1][0] * tau / (y[i1-2][0]-y[i1-1][0]);
    Mcomp = y[i1-1][1] + (y[i1-1][1]-y[i1-2][1]) / tau * (Rcomp-t[i1-1]);
    
    // conversion to solar masses
    Mcomp /= 1.4766;

    if (!flag) {
      diff0 = Mcomp - MR_rel[mcount][1];
      
      if (fabs(diff0) < 1e-6) 
        break;
   /*
      if (diff0 > 0) {
        printf("masses drop with increasing p_center. not yet implemented\n");
        exit(0);
      }
   */
      flag = true;
      continue;
    }

    else {
      diffx = Mcomp - MR_rel[mcount][1];
      
      if (diff > 0.0) {
        if (fabs(diffx - diff) > 0) {
          diff = diffx;
          break;
        }
      }

      diff = diffx;
      if (fabs(diff) < 1e-6) 
        break;

      // printf("diff mass %g %g %g %g %g %g\n",pstep,slope,
      // y[(i1-1)*N+1] / 1.4766, MR_rel[mcount][1],t[i1-1],
      // MR_rel[mcount][0]);

      if (diff * diff0 > 0) {
        // cout << diff * diff0 << endl;
        continue;
      }

      pstep /= 10.0;
      diff = 0.0;

      if (pstep < 1e-10) {
        // cout << "pstep < x breaking condition" << endl;
        break;
      }

      p_end -= 10.0 * pstep;
      continue;
    }

  } // end p_end < p_dur loop


Rcomp = t[i1-1] + y[i1-1][0] * tau / (y[i1-2][0]-y[i1-1][0]);
Mcomp = y[i1-1][1] + (y[i1-1][1]-y[i1-2][1]) / tau * (Rcomp-t[i1-1]);
Mcomp /= 1.4766;
what = i1;

  return;
}

