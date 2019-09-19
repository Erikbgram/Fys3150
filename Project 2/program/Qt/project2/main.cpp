#include <iostream>
#include <cmath>
#include <armadillo>

double s(double theta) {
    // s = sin(theta)
    return std::sin(theta);
}

double c(double theta) {
    // c = cos(theta)
    return std::cos(theta);
}

double t(double theta) {
    // t = tan(theta)
    return std::tan(theta);
}

int* maxdiag(arma::mat A, int n) {
    double max;
    int p;
    int q;
    int* r[2];
    for (int i; i < n; i++) {
        for (int j; j < n; i++) {
            double a_ij = fabs(A(i,j));
            if (a_ij > max) {
                max = a_ij;
                p = i;
                q = j;
            }
        }
    }
    r[0] = p;
    r[1] = q;
    return r;
}

int main() {
    int n = 4;
    double theta = 0;
    arma::mat A(n, n, arma::fill::randu); // Create a random matrix with size nxn
    arma::mat S = diagmat(A);

    A.print("A: ");
    S.print("S: ");


}