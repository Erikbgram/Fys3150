#include <iostream>
#include <vector>

#include <cstdlib>
#include <iomanip>
#include "time.h"

using namespace std;

std::vector<int> even = {2, 4, 3, 6, 1, 9};

int main() {

    int n = 10;
    double* a[n]  = {};
    //double* b[n]  = {};
    double* c[n]  = {};
    double* v[n]  = {};
    double* b_[n] = {};
    double* g[n]  = {};
    double* g_[n] = {};

    vector<double> b[n];

    //b[0] = 1;
    //b[n] = 1;

    int sum = 0

    clock_t start, finish;
    start = clock();
    /*
    for(int i = 1; i<n-1; i++) {
        b_[i] = b[i] - (a[i-1]*c[i-1])/b_[i-1];
        g_[i] = g[i] - (g[i-1]*a[i-1])/b_[i-1];
    }
    for(int i = n-1; i>1; i--) {
        v[i]  = (g_[i] - c[i]*v[i+1])/b_[i];
    }
    */
    for(int i = 1; i<n-1; i++) {
      sum += i;
    }
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC);
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(10) << setw(20) << "Time used for algorithm = " << timeused  << endl;
}
