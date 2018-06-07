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

    while (M <= 1.6) {

      if (fscanf(TOV, "%lf,%lf", &M, &R) == EOF) break; 
     // if (M > 1.6) break;
      mcount++;
    } 
    
//    cout << "mcount = " << mcount << "\n";
//    cout << "M = " << M << "\n";

//    fclose(TOV);
//    TOV == NULL;

//    #pragma omp parallel for private(P,m,p,e,r,dm,dp,de)
    
    m = 0;
    r = pow(10, -14);
//    p = 0.00000001 + mcount * 0.0000001;

//    for (int a = 0; a <= mcount; a++){
//        p = 0.00000001 + a * 0.0000001;     
//        Eresult.push_back(eos(p));     
//        Presult.push_back(p);
//    }
//cout <<mcount<<endl;
    p = 0.00000001 + mcount * 0.0000001;

    int i1,i2; 
    const int N = 2;
    const int nS = 100000;
    double tau = 0.001;
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


/* 
    fstream f;
    f.open("eos.out", ios::out);
    
      for (unsigned int i = 0; i < Eresult.size(); i++) {
          f << setprecision(10) << Presult[i] 
            << "," << setprecision(10) << Eresult[i]  << "\n";
      }

    f.close();

*/


/*   // Massive debug 

	   for (unsigned int i = 0; i < 5; i++) {
		  cout << Presult[i] << "," << Eresult[i]  << "\n";
	   }
	   cout << setprecision(12) << Presult.back() << "," 
		<< setprecision(12) << Eresult.back() << "\n";
	   cout << Presult.size() << " , " << Eresult.size() << "\n";

	   vector<double>::iterator itP;
	    itP = Presult.begin();
	    itP = Presult.insert(itP, Presult[0]+5); 
	 
	   vector<double>::iterator itE;
	    itE = Eresult.begin();
	    itE = Eresult.insert(itE, 6); 

	   cout << "*******************\n";

	   for (unsigned int i = 0; i < 5; i++) {
		  cout << Presult[i] << "," << Eresult[i]  << "\n";
	   }
	   cout << setprecision(12) << Presult.back() << "," 
		<< setprecision(12) << Eresult.back() << "\n";
	   cout << Presult.size() << " , " << Eresult.size() << "\n";

*/




   vector<double>::iterator itP;
    itP = Presult.begin();
    itP = Presult.insert(itP, Presult[0]+0.0000001); 


/*    // debug
    for (unsigned int i = 0; i < 10; i++) {
          cout << setprecision(10) << Presult[i] 
               << "," << setprecision(10) << Eresult[i]  << "\n";
    }
*/

  mcount -= 1;
  double err = 0.0001;
  double reos;
   
//  double y0[N];
 
  while (M <= 1.7) {

      if (fscanf(TOV, "%lf,%lf", &M, &R) == EOF) break; 
     // if (M > 1.7) break;
//      mcount++;
//    }  

         cout << "M = " << M << "\n";

         reos = Eresult[0];
   
         cout << "reos = " << reos << "\n";
//         y[1][i1] = 1.5;
 
         vector<double> mass;
         mass.push_back(M);

	 while (fabs(y[1][i1] - M) > err){           
           

	    p = Presult[0];
            r = pow(10., -14.);
            cout << "p = " << p << "\n";

	    y0[0] = p; y0[1] = 0.0;

	    for (i1 = 0; i1 < N; i1++){
	      y[i1][0] = y0[i1];
	    }
//cout<<"hi"<<endl;
	      double rho = r; double k1[N];
                       // y[0]
	      k1[0] = tov(Presult[0], y[1][0], r) * tau; 
	      k1[1] = 4*M_PI*pow(r, 2.) * reos * tau;
	     
	      for (i2 = 0; i2 < N; i2++){
		y[i2][1] = y[i2][0] + k1[i2];
	      }

	      r = rho + tau;
//cout<<"hi"<<endl;
 
	    for (i1 = 2; i1 <= Eresult.size() 
                 && y[0][i1] > 0. && y[1][i1] > 0.; i1++){

	      rho = r;
                       // y[0]
	      k1[0] = tov(Presult[i1-1], y[1][i1-1], r) * tau; 
	      k1[1] = 4*M_PI*pow(r, 2.) * Eresult[i1-1] * tau;
	     
	      for (i2 = 0; i2 < N; i2++){
		y[i2][i1-1] = y[i2][i1-2] + k1[i2];
	      }
         //          cout << y[0][i1] << "  " << y[1][i1] << endl;
	      r = rho + tau;
	    }
//cout<<"hi"<<endl;

	    if (y[1][i1] < M) reos += 0.00001; //variation
	    if (y[1][i1] >= M) reos -= 0.00001; //variation

            vector<double>::iterator itM;
	    itM = mass.begin();
	    itM = mass.insert(itM, y[1][i1]); 
/*
            cout << mass[0] << "  "  << mass[1] << endl;

            if (mass[0]-mass[1] > err) reos += 0.00001;
            if (mass[0]-mass[1] <= err) reos -= 0.00001;
*/
            cout << "var. reos = " << reos << endl;
           
	    cout << "M2 = " << y[1][i1] << "\n";	 
         
      //     if (fabs(y[1][i1] - M) > err) break;

           // mass.push_back(y[1][i1]);
	   
           // cout << reos << "\n";
             }
   

//   cout << mcount << "\n"; 
	   
   vector<double>::iterator itP;
    itP = Presult.begin();
    itP = Presult.insert(itP, Presult[0]+0.0000001); 
 
   vector<double>::iterator itE;
    itE = Eresult.begin();
    itE = Eresult.insert(itE, reos); 

   mcount++;
}

/* 
    fstream f;
    f.open("eos.out", ios::out);
    
      for (unsigned int i = 0; i < Eresult.size(); i++) {
          f << setprecision(10) << Presult[i] 
            << "," << setprecision(10) << Eresult[i]  << "\n";
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

   fclose(TOV); 
   TOV = NULL;
   return 0;   
}

