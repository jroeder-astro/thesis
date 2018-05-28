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
      int Esize = pow(Pmax,2.)-Pmax;
      Eresult.resize(Esize);
    vector<double> Presult;   
      int Psize = pow(Pmax,2.)-Pmax;
      Presult.resize(Psize);
    vector<double> Mresult;
      int Msize = pow(Pmax,2.)-Pmax;
      Mresult.resize(Msize);
    vector<double> Rresult;
      int Rsize = pow(Pmax,2.)-Pmax;
      Rresult.resize(Rsize); 

//  double Eresult[Pmax], Presult[Pmax], Mresult[Pmax], Rresult[Pmax];   

    double m,p,e,r,dm,dp,de,dr;
    double eos(double);
    dr = 0.00001;    

    unsigned int i = 0;
    int P = 0;
    
//    #pragma omp parallel for private(P,m,p,e,r,dm,dp,de)

    int i1, i1;
    const int N = 2;
    
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
 
    for (i1 = 0; i1 < N; i1++){

    }






















    do {
        e = eos(p);
        //  Eresult.push_back(e);
        //  Presult.push_back(p);                                     
        dm = 4*M_PI*e*r*r*dr;                         
        dp = tov(p,m,r)*dr; 
        r = r+dr; m = m+dm; p = p+dp;                                       
    }
    while (p>0);
   
    Mresult[P] = m;
    Rresult[P] = r; 

//    Mresult.push_back(m);    
//    Rresult.push_back(r);
 
    i++;

    //cout << Mresult[i] << "," << Rresult[i]  << "\n";
    //cout<<"Neutronensternradius [km]          = "<<r<<"\n";
    //cout<<"Neutronensternmasse [Sonnenmassen] = "<<m/1.4766<<"\n";
    }
 
    fstream MR;
    MR.open("tov.out", ios::out);
     
       for (unsigned int j = 0; j < Pmax; j++) {
          MR << Mresult[j] << "," << Rresult[j]  << "\n";
       }
    MR.close();


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

