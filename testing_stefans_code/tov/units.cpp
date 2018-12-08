#include<iostream>
#include<stdlib.h>
#include<stdio.h>
using namespace std;

main() {

  double e, p, u1, u2;
  double convert = (1.3234)*(1e-6);
 
  // MeV/fm^3 = convert * 1/km^2
  
  FILE *eos = fopen("table.csv", "r");
  FILE *out = fopen("output.dat", "w");

  if (eos == NULL) exit(0);
 
  while (1) {
    if (fscanf(eos, "%lf %lf %lf %lf", &e, &p, &u1, &u2) == EOF )
      break;

    e *= convert;
    p *= convert;
    fprintf(out, "%5.20lf,%5.20lf\n", e, p);
  }

  fclose(eos); eos = NULL;
  fclose(out); out = NULL;

  return 0;
}
