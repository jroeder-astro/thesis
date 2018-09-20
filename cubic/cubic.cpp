#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>

using namespace std;



main() {

// **** File i/o ********************

  vector<vector<double>> data;
  vector<double> point(2);

  FILE *infile = fopen("stuff.dat", "r");
  if (infile == NULL) exit(0);
  while (1) {
    if (fscanf(infile, "%lf,%lf", &point[0], &point[1]) == EOF) break;
    data.push_back(point);
  }
  fclose(infile);

  // for (int i = 0; i < data.size(); i++)
  //   cout << data[i][0] << "," << data[i][1] << endl;


// **** Interpolation ***************





  return 0;
}

