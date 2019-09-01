#include <iostream>

using namespace std;

int main() {

    int n = 10;
    int *a  = new int[n];
    int *b  = new int[n];
    int *c  = new int[n];
    int *v  = new int[n];
    int *b_ = new int[n];
    int *g  = new int[n];
    int *g_ = new int[n];



    for (int i = 1; i<n-1;i++) {
        b_[i] = b[i] - (a[i-1]*c[i-1])/b_[i-1];
        g_[i] = g[i] - (g[i-1]*a[i-1])/b_[i-1];
        v[i]  = (g_[i] - c[i]*v[i+1])/b_[i];
    }
}
