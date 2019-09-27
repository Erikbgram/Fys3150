/*
Last edited by: Erlend T. North 22:13 27.09.2019
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

using namespace std;

void offdiag(arma::mat A, int *p, int *q, int n);

void Jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l, int n );

TEST_CASE( "offdiags are tested", "[offdiag]" ) {

    int n = 3;
    arma::mat A, B, C, D = arma::mat(n,n, arma::fill::eye);
    int p, q;


    A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
    A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;
    A(2,0) = 7.0; A(2,1) = 8.0; A(2,2) = 9.0;
    offdiag(A, &p, &q, n);
    double A_offdiag = A(p,q);

    B(0,0) =  1.0; B(0,1) = 12.0; B(0,2) = 6.0;
    B(1,0) = 22.0; B(1,1) =  6.0; B(1,2) = 4.0;
    B(2,0) = 25.0; B(2,1) =  8.0; B(2,2) = 5.0;
    offdiag(B, &p, &q, n);
    double B_offdiag = B(p,q);

    C(0,0) = 2.0; C(0,1) =  1.0; C(0,2) = 0.0;
    C(1,0) = 7.0; C(1,1) = 12.0; C(1,2) = 6.0;
    C(2,0) = 4.0; C(2,1) =  5.0; C(2,2) = 3.0;
    offdiag(C, &p, &q, n);
    double C_offdiag = C(p,q);

    D(0,0) =   72.0; D(0,1) =1043.0; D(0,2) = sqrt(2);
    D(1,0) = 3.1415; D(1,1) = 22/89; D(1,2) =   516.0;
    D(2,0) =  2.718; D(2,1) =   8.0; D(2,2) =   23*23;
    offdiag(D, &p, &q, n);
    double D_offdiag = D(p,q);

    REQUIRE( A_offdiag == 8.0 );
    REQUIRE( B_offdiag == 25.0 );
    REQUIRE( C_offdiag == 7.0 );
    REQUIRE( D_offdiag == 1043.0 );
}
