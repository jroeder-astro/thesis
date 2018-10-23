#include<iostream>
#include<stdlib.h>
#include<stdio.h>
using namespace std;

main() {

  double e, p, u1, u2;

  FILE *eos = fopen("table.csv", "r");
  FILE *out = fopen("output.dat", "w");

  if (eos == NULL) exit(0);
 
  while (1) {
    if (fscanf(eos, "   %lf          %lf         %lf         %lf", &e, &p, &u1, &u2))
      break;
    e *= 5; p += 6;
    printf("%3.6lf,%3.6lf\n", e, p);
  }

  fclose(eos); eos = NULL;
  fclose(out); out = NULL;

  return 0;
}
