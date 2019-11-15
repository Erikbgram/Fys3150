/*
  Last edited: 15.11.2019 13:04 by Erlend T. North
*/

#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>
#include <chrono>
#include <fstream>
//#include <mpi.h>
#define kB 1
#define J 1
using namespace std;
namespace ch = std::chrono;

// Functions
void initLattice(arma::Mat<int> &M, int L) { // Initializes the lattice
  random_device rd;
  mt19937_64 generator(rd());
  uniform_int_distribution<int> distribution(0,1);

  for(int i = 0; i < L; i++) {
    for(int j = 0; j < L; j++) {
      M(i,j) = distribution(generator);
      if(M(i,j) == 0) { // We want spins -1 and 1, not 0 and 1
        M(i,j) = -1;
      }
    }
  }
  M.print();
}

void flip(int &spin) { // Flips spin-value
  if(spin==1) {
    spin = -1;
  }
  else {
    spin = 1;
  }
}

int localE(arma::Mat<int> M, int L, int k, int l) { // Calculates energy around <kl>
  int sum = 0;
  for(int i = -1; i < 2; i++) {
    // k
    if(k==0) {
      sum += M(k+abs(i),l)*M(k,l); // Periodic boundary condition.  abs(i) makes k-1 into k+1
    }
    else if(k==L-1) {
      sum += M(k-abs(i),l)*M(k,l); // Periodic boundary condition. -abs(i) makes k+1 into k-1
    }
    else {
      sum += M(k+i,l)*M(k,l);
    }
    // l
    if(l==0) {
      sum += M(k,l+abs(l))*M(k,l); // Periodic boundary condition.  abs(i) makes l-1 into l+1
    }
    else if(l==L-1) {
      sum += M(k,l-abs(l))*M(k,l); // Periodic boundary condition. -abs(i) makes l+1 into l-1
    }
    else {
      sum += M(k,l+i)*M(k,l);
    }
  }
  return sum;
}

int globalE(arma::Mat<int> M, int L) {
  int sum = 0;
  for(int k = 0; k < L; k++) {
    for(int l = 0; l < L; l++) {
      for(int i = -1; i < 2; i++) {
        // k
        if(k==0) {
          sum += M(k+abs(i),l)*M(k,l); // Periodic boundary condition.  abs(i) makes k-1 into k+1
        }
        else if(k==L-1) {
          sum += M(k-abs(i),l)*M(k,l); // Periodic boundary condition. -abs(i) makes k+1 into k-1
        }
        else {
          sum += M(k+i,l)*M(k,l);
        }
        // l
        if(l==0) {
          sum += M(k,l+abs(l))*M(k,l); // Periodic boundary condition.  abs(i) makes l-1 into l+1
        }
        else if(l==L-1) {
          sum += M(k,l-abs(l))*M(k,l); // Periodic boundary condition. -abs(i) makes l+1 into l-1
        }
        else {
          sum += M(k,l+i)*M(k,l);
        }
      }
    }
  }
  return sum;
}

int sumNeighbours(arma::Mat<int> lattice, int L) {

  return 0;
}

int moro(arma::Mat<int> lattice, int L) {
  //int *energies = new double[-8, 0, 0, -4, 0, 0, 0, 0, 0, 4, 0, 0, 8];
  return 0;
}

int magnetization(arma::Mat<int> lattice, int L) {
  int M = 0;
  for(int k = 0; k < L; k++) {
    for(int l = 0; l < L; l++) {
      M += lattice(k,l);
    }
  }
  return abs(M);
}

int excpectMagnetization(arma::Mat<int> lattice, int L) {

}

double partitionFunction(int *E, int L) {
  double Z = 0;
  for(int i = 0; i < L*L; i++) {
    Z += exp(-E[i]);
  }
  return Z;
}

double susceptibility() {
  return 1;
}

void MonteCarlo(int n, double a, double b, double  &integral, double  &std) {
  // Do stuff
}

bool Metropolis(double delE,double T=1.0) {
  random_device rd;
  mt19937_64 generator(rd());
  uniform_real_distribution<double> distribution(0.0,1.0);

  return 0;
}

int main(int argc, char *argv[]) { // Main function
  int L = atoi(argv[1]);
  int n = atoi(argv[2]);
  int* E = new int[n];
  int* M = new int[n];
  int* X = new int[n];



  //----------------------------------------------------------------------------
  // Initialize lattice
  arma::Mat<int> lattice(L,L,arma::fill::zeros);
  initLattice(lattice, L);

  /// PRECALCULATE E AND M
  E[0] = globalE(lattice, L);
  M[0] = magnetization(lattice, L);

  //----------------------------------------------------------------------------
  // Initiate Monte Carlo Markov Chain
  double T = 1.0;
  bool metropolis_bool = 0;
  int Ei = 8;
  int Ej = 0;
  int deltaE = 0;
  random_device rd;
  mt19937_64 generator(rd());
  uniform_int_distribution<int> distribution(0,L-1);
  uniform_real_distribution<double> Metro(0,1);

  fstream outfile;
  outfile.open("../../energy.txt", fstream::out);
  outfile << "Energy" << endl;

  //----------------------------------------------------------------------------
  // Monte Carlo Markov Chain
  for(int i = 1; i < n; i++) {

    // Choose random spin to flip
    int k = distribution(generator);
    int l = distribution(generator);

    Ej = localE(lattice, L, k, l); // Energy pre-flip

    flip(lattice(k,l)); // Flipping

    Ei = localE(lattice, L, k, l); // Energy post-flip

    deltaE = Ei-Ej; // Energy difference

    // Metropolis
    double r = Metro(generator);
    cout << r << ", " << exp(deltaE/(kB*T)) << endl;
    if(r <= exp(deltaE/(kB*T))) {
      metropolis_bool = 1;
    }
    else  {
      metropolis_bool = 0;
    }

    if(metropolis_bool) {
      cout << "Accepted" << endl;
    }
    else {
      flip(lattice(k,l)); // Flip back
      cout << "Rejected" << endl;
    }

    // Calculating averages
    E[i] = globalE(lattice, L);
    M[i] = magnetization(lattice, L);


    cout << "delE=" << deltaE << endl;
    outfile << deltaE << endl;
    cout << endl;
    lattice.print();
  }
  //cout << "Test" << endl;

  /// THINGS WE NEED
  ///
  /// mean energy E
  /// mean magnetization |M|
  /// specific heat Cv
  /// susceptibility chi (X)

  outfile.close();

}
