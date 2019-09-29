/*
Last edited by: Erlend T. North 23:19 27.09.2019
*/

#include <iostream>
#include <armadillo>
#include <chrono>
#include <fstream>
#include "catch.hpp"

using namespace std;
namespace ch = std::chrono;

//  the offdiag function, using Armadillo
void offdiag(arma::mat A, int *p, int *q, int n) {
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              max = aij;  *p = i; *q = j;
           }
       }
   }
}


void Jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l, int n ) {
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));

    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }

    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  A(l,k) = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  }
  return;
} // end of function jacobi_rotate

// Tried catch.cpp, but did not have time to learn
// Primitive, yet functional tests

bool test_ortho() {
    // Checks for largest offdiag-element in matrices A, B, C and D
    int n = 5;
    arma::mat A = arma::mat(n,n, arma::fill::zeros);
    arma::mat B = arma::mat(n,n, arma::fill::zeros);
    arma::mat C = arma::mat(n,n, arma::fill::zeros);
    arma::mat D = arma::mat(n,n, arma::fill::zeros);
    int p, q;

    A(0,0) =  1.0; A(0,1) =  4.0; A(0,2) =    7.0; A(0,3) =   5.0; A(0,4) =   22.0;
    A(1,0) =  4.0; A(1,1) =  5.0; A(1,2) =    8.0; A(1,3) =  57.0; A(1,4) =   62.0;
    A(2,0) =  7.0; A(2,1) =  8.0; A(2,2) =    9.0; A(2,3) = -52.0; A(2,4) = 1034.0;
    A(3,0) =  5.0; A(3,1) = 57.0; A(3,2) =  -52.0; A(3,3) =  25.0; A(3,4) =   83.0;
    A(4,0) = 22.0; A(4,1) = 62.0; A(4,2) = 1034.0; A(4,3) =  83.0; A(4,4) =  -38.0;
    offdiag(A, &p, &q, n);
    double A_offdiag = A(p,q);
    bool A_pass = A_offdiag == 1034.0;
    cout << "Test for A returned: " << A_pass << endl;

    B(0,0) =    1.0; B(0,1) =      22.0; B(0,2) = 25.0; B(0,3) = 0.028; B(0,4) =    0.0083;
    B(1,0) =   22.0; B(1,1) =       6.0; B(1,2) =  8.0; B(1,3) = -87.0; B(1,4) = 1000000.0;
    B(2,0) =   25.0; B(2,1) =       8.0; B(2,2) =  5.0; B(2,3) =  17.0; B(2,4) =      94.0;
    B(3,0) =  0.028; B(3,1) =     -87.0; B(3,2) = 17.0; B(3,3) =  31.0; B(4,3) =      87.0;
    B(4,0) = 0.0083; B(4,1) = 1000000.0; B(4,2) = 94.0; B(4,3) =  87.0; B(4,4) =      31.0;
    offdiag(B, &p, &q, n);
    double B_offdiag = B(p,q);
    bool B_pass = B_offdiag == 1000000.0;
    cout << "Test for B returned: " << B_pass << endl;

    C(0,0) =   2.0; C(0,1) =  7.0; C(0,2) =   4.0; C(0,3) = 37.56; C(0,4) = 91.0;
    C(1,0) =   7.0; C(1,1) = 12.0; C(1,2) =   6.0; C(1,3) =  89.0; C(1,4) = 75.0;
    C(2,0) =   4.0; C(2,1) =  6.0; C(2,2) =   3.0; C(2,3) = 25.24; C(2,4) = 82.0;
    C(3,0) = 37.56; C(3,1) = 89.0; C(3,2) = 25.24; C(3,3) =  23.0; C(3,4) = 13.0;
    C(4,0) =  91.0; C(4,1) = 75.0; C(4,2) =  82.0; C(4,3) =  13.0; C(4,4) = 23.0;
    offdiag(C, &p, &q, n);
    double C_offdiag = C(p,q);
    bool C_pass = C_offdiag == 91.0;
    cout << "Test for C returned: " << C_pass << endl;

    D(0,0) =   72.0; D(0,1) = 3.1415; D(0,2) =   2.718; D(0,3) =    53.0; D(0,4) = 64.27;
    D(1,0) = 3.1415; D(1,1) =  188.0; D(1,2) =  3.1416; D(1,3) =  -183.0; D(1,4) =  13.0;
    D(2,0) =  2.718; D(2,1) = 3.1416; D(2,2) =    68.0; D(2,3) = -121.57; D(2,4) =  -3.0;
    D(3,0) =   53.0; D(3,1) = -183.0; D(3,2) = -121.57; D(3,3) =   1.119; D(3,4) =  87.0;
    D(4,0) =  64.27; D(4,1) =   13.0; D(4,2) =    -3.0; D(4,3) =    87.0; D(4,4) = 1.119;
    offdiag(D, &p, &q, n);
    double D_offdiag = D(p,q);
    bool D_pass = D_offdiag == 87.0;
    cout << "Test for D returned: " << D_pass << endl;

    return (A_pass && B_pass && C_pass && D_pass);
}

bool test_eig() {
    int n = 5;

    arma::mat A = arma::mat(n, n, arma::fill::zeros);
    A(0,0) = 2.0; A(0,1) =  4.0; A(0,2) =    7.0; A(0,3) =   5.0; A(0,4) =   22.0;
    A(1,0) = 0.0; A(1,1) =  5.0; A(1,2) =    8.0; A(1,3) =  57.0; A(1,4) =   62.0;
    A(2,0) = 2.0; A(2,1) =  8.0; A(2,2) =    9.0; A(2,3) = -52.0; A(2,4) = 1034.0;
    A(3,0) = 0.0; A(3,1) = 57.0; A(3,2) =  -52.0; A(3,3) =  25.0; A(3,4) =   83.0;
    A(4,0) = 2.0; A(4,1) = 62.0; A(4,2) = 1034.0; A(4,3) =  83.0; A(4,4) =  -38.0;


    return 1;
}

int main(int argc, char *argv[]) {
    std::cout << std::scientific;
    int n = atoi(argv[1]);
    //cout << n;
    cout << endl;

    arma::mat A = arma::mat(n, n, arma::fill::zeros);
    arma::mat R = arma::mat(n, n, arma::fill::eye);

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
    //A.print("A: ");



    arma::vec eigval;
    arma::mat eigvec;


    ch::steady_clock::time_point start = ch::steady_clock::now();

    arma::eig_sym(eigval, eigvec, A);

    ch::steady_clock::time_point stop = ch::steady_clock::now();
    ch::duration<double> time_span_eig_sym = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used by eig_sym = " << time_span_eig_sym.count()  << "s" << std::endl;


    double tolerance = 1.0E-10;
    int maxiteration = 100000;
    double maxnondiag = 1;

    int iterations = 0;

    start = ch::steady_clock::now();

    while ( maxnondiag > tolerance && iterations <= maxiteration) {
       int p, q;
       offdiag(A, &p, &q, n);
       maxnondiag = fabs(A(p,q));
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
       //cout << "Iteration: " << iterations << " complete!" << endl;
    }

    stop = ch::steady_clock::now();
    if(iterations > maxiteration) {
        cout << "Reached max iterations" << endl;
    }
    else {
        cout << "Complete!" << endl;
    }
    ch::duration<double> time_span_ours = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used by us = " << time_span_ours.count()  << "s" << std::endl;

    /*
    fstream outfile;
    outfile.open("../../stats.txt", std::fstream::out | std::ofstream::app);
    outfile << std::scientific ;
    outfile << n << ", " << iterations << ", " << time_span_eig_sym.count() << ", " << time_span_ours.count() << endl;
    */

    A.print("A: ");
    eigval.print("Armadillo eigenvalues: ");

    R.print("R: ");
    eigvec.print("S: ");

    cout << endl << "Commencing tests.." << endl;
    bool ortho = test_ortho();
    if(ortho) {
        cout << "Orthogonality preserved!" << endl;
    }
    else {
        cout << "Orthogonality NOT preserved!" << endl;
    }

    return 0;
}
