#include <iostream>
#include <armadillo>

using namespace std;


int main(){


int n = 5;
int pn = 10;
int po = 0; 
double h = pn - po/n;




double* p = new double [10];




for (int i = 0; i < n; i++){

p[i] = p[0] + h*i;



}

arma::mat A(n,n,arma::fill::zeros);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(i == j) {
          A(i,j) = 2/(h*h) + p[i]*p[i];
        }
        else if(i == j+1) {
            A(i,j) =(-1)/(h*h) ;
        }
        else if(i == j-1) {
            A(i,j) = (-1)/(h*h);
        }
    }
    }



A.print("A: ");



}





































