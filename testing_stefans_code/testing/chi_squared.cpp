#include<iostream>
#include<math.h>
#include<vector>
using namespace std;

double eos(double p);

main() {

// *** VARIABLES ***

  double R=0., M=0., p=0., e=0., Rd=0., Md=0., s=0.;
  vector<double> data_one(4);
  vector<vector<double>> data;
  vector<double> rdata_one(7);
  vector<vector<double>> rdata;
  int i1=0, i2=0;
  double cs=0.;

// *** FILE I/O: reconstruction output***

  FILE *recon = fopen("results_l2_7", "r");
  if (recon == NULL) 
    exit(0);
  else
    cout << "rdata open\n";

  while (1) {
  
    if (fscanf(recon, "%lf %lf %lf %lf %lf %lf %lf", &Rd, &R, &Md, &M, &p, &e ,&s) == EOF)
      break;

    printf("%lf %lf %lf %lf %lf %lf %lf\n", Rd, R, Md, M, p, e ,s);

    rdata_one[0] = Rd; rdata_one[1] = R; 
    rdata_one[2] = Md; rdata_one[3] = M; 
    rdata_one[4] = p;  rdata_one[5] = e; 
    rdata_one[6] = s; 
    rdata.push_back(rdata_one);

  }

  fclose(recon); 
  recon = NULL;

// *** FILE I/O: default values ***

  FILE *input = fopen("mr.out", "r");
  if (input == NULL) 
    exit(0);
  else
    cout << "data open\n";

  while (1) {
  
    if (fscanf(input, "%lf %lf %lf %lf", &R, &M, &p, &e) == EOF)
      break;

    printf("%lf %lf %lf %lf\n", R, M, p, e);

    data_one[0] = R; data_one[1] = M; 
    data_one[2] = p; data_one[3] = e;
    data.push_back(data_one);

  }

  fclose(input); 
  input = NULL;

// *** CHI SQUARED ***

  for (i1 = 0; i1 < rdata.size(); i1++){
    cs += pow( (data[i1+7][2]-rdata[i1][4]), 2. ) / data[i1+7][2];
    cout << cs << endl;
  }

// *** OUTPUT ***

  FILE *output = fopen("cs.out", "a");
  if (output == NULL) 
    exit(0);

  fprintf(output, "%5.10lf\n", cs);

  fclose(output);
  output = NULL;

  return 0;
}

double eos(double p) {
  return pow(p/10., 3./5.);
}


