#include<iostream>                 
#include<math.h>                   
#include<vector> 
#include<fstream>
#include<omp.h>
#include<iomanip>
#include<stdlib.h>

using namespace std;                

double eos(double p){
    return pow(p/10,3.0/5);
}  

double tov(double p, double m, double r){
    return -(eos(p)+p)*(m+4*M_PI*r*r*r*p)/(r*(r-2*m));
}

int main(void){
   
    int Pmax = 1001;

    vector<double> Eresult;    
   //   int Esize = pow(Pmax,2.)-Pmax;
   //   Eresult.resize(Esize);
    vector<double> Presult;   
   //   int Psize = pow(Pmax,2.)-Pmax;
   //   Presult.resize(Psize);
/*
    vector<double> Mresult;
      int Msize = pow(Pmax,2.)-Pmax;
      Mresult.resize(Msize);
    vector<double> Rresult;
      int Rsize = pow(Pmax,2.)-Pmax;
      Rresult.resize(Rsize); 
*/
//  double Eresult[Pmax], Presult[Pmax], Mresult[Pmax], Rresult[Pmax];   

    double m,p,e,r,dm,dp,de,dr;
    double eos(double);
    dr = 0.00001;    

    unsigned int i = 0;
    int P = 0;
    
    double M, R;

    FILE *TOV = fopen("tov.out", "r");
    if (TOV == NULL) exit(0);
    int mcount = 0;
    while (1) {

      if (fscanf(TOV, "%lf,%lf", &M, &R) == EOF) break; 
      if (M > 1.6) break;
      mcount++;
    } 

//    cout << "mcount = " << mcount << "\n";
//    cout << "M = " << M << "\n";

    fclose(TOV);
    TOV == NULL;

//    #pragma omp parallel for private(P,m,p,e,r,dm,dp,de)
    
    m = 0;
    r = pow(10, -14);
    p = 0.00000001 + mcount * 0.0000001;

    int i1,i2; 
    const int N = 2;
    const int nS = 100000;
    double tau = 0.01;
    double y[N][nS+1];    
    double y0[N] = {p, m};

    for (i1 = 0; i1 < N; i1++){
      y[i1][0] = y0[i1];
    }

    for (i1 = 1; i1 <= nS && y[0][i1-1] > 0.; i1++){
      double rho = r; double k1[N];

      k1[0] = tov(y[0][i1-1], y[1][i1-1], r) * tau; 
      k1[1] = 4*M_PI*pow(r, 2.) * eos(y[0][i1-1]) * tau;
     
      Eresult.push_back(eos(y[0][i1-1]));     
      Presult.push_back(y[0][i1-1]);

      for (i2 = 0; i2 < N; i2++){
        y[i2][i1] = y[i2][i1-1] + k1[i2];
      }

    r = rho + tau;
    }

//   for (unsigned int i = 0; i < Presult.size(); i++) {
//          cout << Presult[i] << "," << Eresult[i]  << "\n";
//   }


    vector<double>::iterator itP;
    itP = Presult.begin();
    itP = Presult.insert(itP, Presult[0]+5.); 

    for (unsigned int i = 0; i < 10; i++) {
          cout << setprecision(10) << Presult[i] 
               << "," << setprecision(10) << Eresult[i]  << "\n";
    }


/* 
    fstream f;
    f.open("eos.out", ios::out);
    
      for (unsigned int i = 0; i < 100; i++) {
          f << Presult[i] << "," << Eresult[i]  << "\n";
      }

    f.close();
*/
/*
    cout << "*********************" << endl;

    vector<double>::iterator itP;
    itP = Presult.begin();
    itP = Presult.insert(itP, Presult[0]+5.); 

    for (unsigned int i = 0; i < 10; i++) {
          cout << setprecision(10) << Presult[i] 
               << "," << setprecision(10) << Eresult[i]  << "\n";
    }

*/
/* 
    fstream MR;
    MR.open("tov.out", ios::out);
     
       for (unsigned int j = 0; j < Pmax; j++) {
          MR << Mresult[j] << "," << Rresult[j]  << "\n";
       }
    MR.close();
*/

/*
    cout<<"Neutronensternradius [km]          = "<<r<<"\n";
    cout<<"Neutronensternmasse [Sonnenmassen] = "<<m/1.4766<<"\n";
    
    fstream f;
    f.open("eos.out", ios::out);
    
      for (unsigned int i = 0; i < Eresult.size(); i++) {
          f << Presult[i] << "," << Eresult[i]  << "\n";
      }

    f.close();
*/
    return 0;   
}

