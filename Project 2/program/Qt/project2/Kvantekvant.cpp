#include <iostream>
#include <armadillo>

using namespace std;


int main(){


int n = 10;
int pn = 10;
int po = 0; 
double h = pn - po/n;




double* p = new double [10];




for (int i = 0; i < n; i++){

p[i] = p[0] + h*i;

cout << p[i] << endl; 

}




}




































