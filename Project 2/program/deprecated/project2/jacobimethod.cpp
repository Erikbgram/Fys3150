/*
Last edited by Erlend 27.09.2019 13:40
*/

#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <fstream>

namespace ch = std::chrono;

int* maxnondiag(arma::mat A, int n) {
    //finds the largest element in A and returns its indices
    double max = 0;
    int p = 0;
    int q = 0;
    int* indices = new int[2];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a_ij = fabs(A(i,j));
            if (i!=j && a_ij > max) { //if not diagonal & element is larger than current max
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

arma::mat Jacobi_rotate ( arma::mat A, arma::mat R, int k, int l, int n )
{
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
  return A;
} // end of function jacobi_rotate

int main() {
    double eps = 1.0E-8;
    int n;
    std::cout << "The tolerance is: " << eps << std::endl;
    std::cout << "Choose n: ";
    std::cin >> n;
    std::cout << std::endl;
    arma::mat A = arma::mat(n, n, arma::fill::zeros);
    arma::mat R = arma::mat(n, n, arma::fill::eye);

    //Creates the matrix A
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


    arma::mat S = diagmat(A); //"exact" solution

    int* index = maxnondiag(A, n);
    int k = index[0];
    int l = index[1];
    delete[] index;

    int iteration = 0;

    //A.print("A: ");
    std::cout << "row: " << k << " col: " << l << std::endl;
    std::cout << "iteration: " << iteration << std::endl;
    std::cout << "A(k,l): " << A(k,l) << "\nA(k,l)^2: " << fabs(A(k,l)*A(k,l)) << std::endl;
    std::cout << std::endl;
    int iterations = 0;
    while (fabs(A(k,l)*A(k,l)) > eps) {
        A = Jacobi_rotate(A, R, k, l, n);
        int* index = maxnondiag(A, n);
        int k = index[0];
        int l = index[1];
        iterations++;
        A.print("A: ");
        std::cout << "row: " << k << " col: " << l << std::endl;
        std::cout << "Rotation " << iterations << " complete!" << std::endl;
    }



    /*
    ch::steady_clock::time_point start = ch::steady_clock::now();
    while( fabs(A(k,l)*A(k,l)) > eps) {
        //defines tau, tan, cos and sin
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));
        double t;
        if (tau > 0) {
            t = 1.0/(tau+std::sqrt(1.0+tau*tau)); //tan
        }
        else {
            t = -1.0/(-tau+std::sqrt(1.0+tau*tau)); //tan
        }
        double c = 1.0/std::sqrt(1.0+tau*tau); //cos
        double s = t*c; //sin


        //create B
        arma::mat B = arma::mat(n, n, arma::fill::zeros);
        for(int i = 0; i < n; i++) {
            if(i != k && i != l) {
                B(i,i) = A(i,i);
                B(i,k) = A(i,k)*c - A(i,l)*s;
                B(i,l) = A(i,l)*c + A(i,k)*s;
                B(k,i) = B(i,k);
                B(l,i) = B(i,l);
            }
            double temp = 2*A(k,l)*c*s;
            B(k,k) = A(k,k)*c*c - temp + A(l,l)*s*s;
            B(l,l) = A(l,l)*c*c + temp + A(k,k)*s*s;
            B(k,l) = 0;
            //B(k,l) = (A(k,k) - A(l,l))*c*s + A(k,l)*(c*c - s*s);
            B(l,k) = B(k,l);
        }
        A = B;
        iteration++;
        index = maxnondiag(A, n);
        k = index[0];
        l = index[1];
        delete[] index;
        //A.print("A: ");
        std::cout << "row: " << k << " col: " << l << std::endl;
        std::cout << "iteration: " << iteration << std::endl;
        std::cout << "A(k,l): " << A(k,l) << "\nA(k,l)^2: " << fabs(A(k,l)*A(k,l)) << std::endl;
        std::cout << std::endl;
    }
    ch::steady_clock::time_point stop = ch::steady_clock::now();

    ch::duration<double> time_span = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used = " << time_span.count()  << "s" << std::endl;
    std::cout << "Iterations = " << iteration << std::endl;
    */

    /* FILEWRITING
    std::string filename = "../../stats.txt";
    std::fstream outfile;

    outfile.open(filename, std::fstream::out | std::ofstream::app);
    outfile << n << " " << iteration << std::endl;
    */

    return 0;
}

//12: 1.31E-05
