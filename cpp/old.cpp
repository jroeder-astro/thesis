#include<iostream>                 
#include<math.h>                   
#include<vector>
#include<fstream>
#include<omp.h>

using namespace std;                

double eos(double p){
    return pow(p/10,3.0/5);
}  

double tov(double p, double m, double r){
    return -(eos(p)+p)*(m+4*M_PI*r*r*r*p)/(r*(r-2*m));
}

int main(void){
   
    int Pmax = 15001;

    vector<double> Eresult;    
     // int Esize = pow(Pmax,2.)-Pmax;
     // Eresult.resize(Esize);
    vector<double> Presult;   
     // int Psize = pow(Pmax,2.)-Pmax;
     // Presult.resize(Psize);
    vector<double> Mresult;
     // int Msize = pow(Pmax,2.)-Pmax;
     // Mresult.resize(Msize);
    vector<double> Rresult;
     // int Rsize = pow(Pmax,2.)-Pmax;
     // Rresult.resize(Rsize); 

//  double Eresult[Pmax], Presult[Pmax], Mresult[Pmax], Rresult[Pmax];   

    double m,p,e,r,dm,dp,de,dr;
    double eos(double);
    dr = 0.00001;    

    unsigned int i = 0;
    int P = 0;
    
//    #pragma omp parallel for private(P,m,p,e,r,dm,dp,de)

    int i1, i2;
    const int N = 2;
    const int nS = 1000;
    double y[N][nS-1];     
    double tau = 0.01; 

    for (int P = 0; P <= Pmax; P++){

    m = 0;
    r = pow(10,-14);
    p = 0.00000001 + P * 0.0000001;

/*
    vector<double> y0;
      y0.push_back(p);
      y0.push_back(m);

    vector<vector<double>> y;
    y.push_back(y0);

*/

    double y0[N] = {p, m};
    double rho = 0.0;    
    double k1[N];

    for (i1 = 0; i1 < N; i1++){
      y[i1][0] = y0[i1];
    }

    for (i1 = 1; i <= nS; i1++){
      rho = r;
      k1[0] = tov(y[0][i1-1], y[1][i1-1], r) * tau; 
      k1[1] = 4*M_PI*pow(r, 2.) * eos(y[0][i1-1]) * tau;

      for (i2 = 0; i2 < N; i2++){
        y[i2][i1] = y[i1][i1-1] + k1[i2];
      }

    r = rho + tau;
    }

    Rresult.push_back((i1-1)*tau);    
    Mresult.push_back(y[1][i1-1]);
 
    i++;

    }
 
    fstream MR;
    MR.open("old.out", ios::out);
       for (unsigned int j = 0; j < Mresult.size(); j++) {
          MR << Mresult[j] << "," << Rresult[j]  << "\n";
       }
    MR.close();

/*  
    fstream f;
    f.open("eos.out", ios::out);
    
      for (unsigned int i = 0; i < Eresult.size(); i++) {
          f << Presult[i] << "," << Eresult[i]  << "\n";
      }

    f.close();
*/

    return 0;   
}

