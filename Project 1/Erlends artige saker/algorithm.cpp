#include <iostream>
#include <ctime>
#include <stdlib.h>

using namespace std;

int main() {
  int n;
  cout << "What should n be?\n";
  cin >> n;

  double* a  = new double[n];
  double* d  = new double[n];
  double* c  = new double[n];
  double* v  = new double[n];
  double* d_ = new double[n];
  double* b  = new double[n];
  double* b_ = new double[n];

  /*
  Av = b

    d0 c0  0
A = a0 d1 c1  0
     0 a1 d2 ..
    ..  0 .. ..

    v0
v = ..
    vn-1

    v0
b = ..
    bn-1
  */

  v[0] = 0;
  v[n-1] = 0;
  b[0] = 100;
  b[n-1] = 100;

  // Setting up arrays
  for(int i = 0; i < n-1; i++) {
    a[i] = -1;
    d[i] = 2;
    c[i] = -1;
  }

clock_t start, stop;
start = clock();
for(int i = 1; i < n-1; i++) {
  d_[i] = d[i] -(a[i-1]*c[i-1])/d_[i-1];
  b_[i] = b[i] - (b_[i-1]*a[i-1])/d_[i-1];
}
for(int i = n-1; i > 1; i--) {
  v[i] = (b_[i] - c[i]*v[i+1])/d[i];
}
stop = clock();
double timeused;
timeused = (double) (start-stop)/(CLOCKS_PER_SEC);
cout << "Time used for loop = " << timeused << "s" << endl;

cout << "v is:\n";
for(int i = 0; i < n-1; i++) {
  cout << v[i] << endl;
}


}
