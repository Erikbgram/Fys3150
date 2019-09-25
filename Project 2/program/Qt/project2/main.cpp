/*
Last edited by Erlend 20.09.2019 16:05
*/

#include <iostream>
#include <armadillo>
#include <cmath>

int* maxdiag(arma::mat A, int n) {
    //finds the largest element in A and returns its indices
    double max = 0;
    int p = 0;
    int q = 0;
    int* indices = new int[2];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a_ij = fabs(A(i,j));
            if (a_ij > max) {
                max = a_ij;
                p = i;
                q = j;
            }
        }
    }
    indices[0] = p;
    indices[1] = q;
    return indices;
}

void rotate(arma::mat A, int n, int k, int l, double theta) {
    /*
    For en gitt k og l
    Gjør en loop for å fylle matrisen B
    */
    arma::mat B = arma::mat(n, n, arma::fill::zeros);
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
    /*B(k,l) = (A(k,k) - A(l,l))*std::cos(theta)*std::sin(theta)
           + A(k,l)*(std::cos(theta)*std::cos(theta) - std::sin(theta)*std::sin(theta));*/
    }
}

int main() {
    double eps = 1.0E-10;
    int n;
    std::cout << "The tolerance is: " << eps << std::endl;
    std::cout << "Choose n: ";
    std::cin >> n;

    arma::mat A = arma::mat(n, n, arma::fill::randu);
    A.print("A: ");
    /*
     * Fill matrix A here
    */
    arma::mat S = diagmat(A);

    int* index = maxdiag(A, n);
    int k = index[0];
    int l = index[1];
    std::cout << k << " " << l << std::endl;
    /*while(A(k, l)*A(k, l) > eps) {
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0) {

        }
    }*/
}
