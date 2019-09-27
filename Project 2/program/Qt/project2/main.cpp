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
    double eps = 1.0E-15;
    int n;
    std::cout << "The tolerance is: " << eps << std::endl;
    std::cout << "Choose n: ";
    std::cin >> n;

    arma::mat A = arma::mat(n, n, arma::fill::randn);

    //Creates the matrix A
    /*for(int i = 0; i < n; i++) {
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
    }*/


    arma::mat S = diagmat(A); //"exact" solution

    int* index = maxdiag(A, n);
    int k = index[0];
    int l = index[1];
    delete[] index;

    int iteration = 0;
    while( fabs(A(k,l)*A(k,l)) > eps) {
        A.print("A: ");
        std::cout << "row: " << k+1 << " col: " << l+1 << std::endl;

        //defines tau, tan, cos and tan
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));
        double t;
        if (tau > 0) {
            t = 1/(tau+std::sqrt(1+tau*tau)); //tan
        }
        else {
            t = -1/(-tau+std::sqrt(1+tau*tau)); //tan
        }
        double c = 1/std::sqrt(1+tau*tau); //cos
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
        std::cout << "iteration: " << iteration << std::endl;
        std::cout << "A(k,l): " << std::scientific << A(k,l) << "\nA(k,l)^2: " << fabs(A(k,l)*A(k,l)) << std::endl;
        iteration++;
        index = maxdiag(A, n);
        k = index[0];
        l = index[1];
        delete[] index;
    }
    return 0;
}
