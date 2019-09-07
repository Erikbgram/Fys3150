#include <iostream> // cout, cin, etc.
#include <ctime> // clock-stuff

using namespace std;


int main() {
  int n = 10;
  double* a = new double[n];

  a[0] = 0;
  a[n-1] = 0;

  clock_t start, stop;
  start = clock();

  for(int i = 1; i < n-1; i++) {
    a[i] = 85 + i;
    cout << i << "  " << a[i] << endl;
  }

  stop = clock();
  double timeused = (double) (stop-start)/(CLOCKS_PER_SEC);
  cout << "Time used for loop = " << timeused << "s" << endl;

  cout << "index  position  value\n";
  for(int i = 0; i < n; i++) {
    cout << "  " << i << "      " << i+1 << "       " << a[i] << endl;
  }

}
