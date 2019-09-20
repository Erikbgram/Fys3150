/*
Last edited by Erlend 20.09.2019 16:05
*/

#include <armadillo>
#include <cmath>

void rotate(arma::mat A, int n, int k, int l, double theta) {
  /*
  For en gitt rad k, og kolonne l
  Gjør en loop for å fylle matrisen B
  */
  arma::mat B = mat(n,n,arma::fill::zeros);
  for(int i = 0; i < n; i++) {
    if(i != k && i != l) {
      B(i,i) = A(i,i);
      B(i,k) = A(i,k)*std::cos(theta) - A(i,l)*std::sin(theta);
      B(i,l) = A(i,l)*std::cos(theta) + A(i,k)*std::sin(theta);
    }
    double temp = 2*A(k,l)*std::cos(theta)*std::sin(theta);
    B(k,k) = A(k,k)*std::cos(theta)*std::cos(theta)
           - temp
           + A(l,l)*std::sin(theta)*std::sin(theta);
    B(l,l) = A(l,l)*std::cos(theta)*std::cos(theta)
           + temp
           + A(k,k)*std::sin(theta)*std::sin(theta);
    B(k,l) = (A(k,k) - A(l,l))*std::cos(theta)*std::sin(theta)
           + A(k,l)*(std::cos(theta)*std::cos(theta) - std::sin(theta)*std::sin(theta));
  }
}

int main() {
  double eps = 1.0E-10;
  std::cout << eps << std::endl;
}
