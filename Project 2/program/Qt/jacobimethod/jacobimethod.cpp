#include <iostream>
#include <armadillo>
#include <chrono>
#include <fstream>

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


void Jacobi_rotate (arma::mat &A, arma::mat &R, int k, int l, int n ) {
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

int main(int argc, char *argv[])
{
    int n = atoi(argv[1]);
    cout << n;
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
    A.print("A: ");


    //  The final matrix R has the eigenvectors in its row elements, it is set to one
    //  for the diagonal elements in the beginning, zero else.

    arma::vec eigval;
    arma::mat eigvec;


    ch::steady_clock::time_point start = ch::steady_clock::now();

    arma::eig_sym(eigval, eigvec, A);

    ch::steady_clock::time_point stop = ch::steady_clock::now();
    ch::duration<double> time_span = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used by eig_sym = " << time_span.count()  << "s" << std::endl;


    double tolerance = 1.0E-10;
    int maxiteration = 100000;
    double maxnondiag = 1;

    int iterations = 0;

    start = ch::steady_clock::now();

    while ( maxnondiag > tolerance && iterations <= maxiteration)
    {
       cout << "Iteration: " << iterations << endl;
       int p, q;
       offdiag(A, &p, &q, n);
       maxnondiag = fabs(A(p,q));
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
    }

    stop = ch::steady_clock::now();
    time_span = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used by us = " << time_span.count()  << "s" << std::endl;


    fstream outfile;
    outfile.open("../../stats.txt", std::fstream::out | std::ofstream::app);
    outfile << n << ", " << iterations << ", " << time_span.count() << endl;

    A.print("A: ");
    eigval.print("Armadillo eigenvalues: ");

    R.print("R: ");
    eigvec.print("S: ");
    return 0;
}
