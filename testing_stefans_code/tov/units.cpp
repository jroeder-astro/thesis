#include<iostream>
#include<stdlib.h>
#include<stdio.h>
using namespace std;

main() {

  double e, p, u1, u2;
  double mfm = 0; // MeV/fm^3
  double cgs = 0;

  FILE *eos = fopen("table.csv", "r");
  FILE *out = fopen("output.dat", "w");

  if (eos == NULL) exit(0);
 
  while (1) {
    if (fscanf(eos, "%lf %lf %lf %lf", &e, &p, &u1, &u2) == EOF )
      break;
   // e *= ;
   // p *= ;
    printf("%3.6lf,%3.6lf\n", e, p);
  }

  fclose(eos); eos = NULL;
  fclose(out); out = NULL;

  return 0;
}
