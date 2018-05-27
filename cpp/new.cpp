#include<iostream>
#include<fstream>
using namespace std;

main(void){
  
  fstream f;
  f.open("tov.out", ios::in);
    string item;
    int count = 0;
    while (count < 20){
      f >> item;
      cout << item << "\n";
      count++;
    }
  
  return 0;
}
