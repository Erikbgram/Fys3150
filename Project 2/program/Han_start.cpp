#include <chrono> // timer
#include <iostream> // 
#include <armadillo>




/*

c = 1/


if (A[i][j] != 0)

tau = A[i][i] - A[j][j] / (2*A[i][j])

    if (tau > 0){

    t = 1/ (tau - sqrt(1 + tau*tau) )
    }
    else
    {

    t = -1/(-tau - sqrt(1+tau*tau))

 
    }

else
{
    c = 1 
    s = 0 
}


A[i][i] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll 
A[j][j] =
A[i][j] = 0
A[j][i] = 0

*/


int main(){


int n;

std::cin >> n; 



     
    
    arma::mat A(n,n,arma::fill::zeros);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(i == j) {
          A(i,j) = 2;
        }
        else if(i == j+1) {
            A(i,j) = -1;
        }
        else if(i == j-1) {
            A(i,j) = -1;
        }
    }
    }



std::cout << A;


}

  arma::mat B = arma::mat(n, n, arma::fill::zeros);
  for(int i = 0; i < n; i++) {
    for 
    if(i != k && i != l) {
      B(i,i) = A(i,i);
      B(i,k) = A(i,k)*std::cos(theta) - A(i,l)*std::sin(theta);
      B(i,l) = A(i,l)*std::cos(theta) + A(i,k)*std::sin(theta);
    }
    
    else
    {
     B =
    }
*/













/*
ch::steady_clock::time_point start = ch::steady_clock::now();






ch::steady_clock::time_point stop = ch::steady_clock::now();


ch::duration<double> time_span = ch::duration_cast<ch::nanoseconds>(stop - start);
  std::cout << std::setprecision(10) << std::setw(20) << "Time used for algorithm = " << time_span.count()  << "s" << std::endl;

*/


