/*
  Last edited: 01.11.2019 14:50 by Erlend T. North
*/

#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>
#include <chrono>
#include <fstream>
//#include <mpi.h>
#define kB 1
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

int localE(arma::Mat<int> M, int L, int k, int l) { // Calculates energy
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

void MonteCarlo(int n, double a, double b, double  &integral, double  &std) {
  // Do stuff
}

bool Metropolis(double delE,double T=1.0) {
  random_device rd;
  mt19937_64 generator(rd());
  uniform_real_distribution<double> distribution(0.0,1.0);

  double r = distribution(generator);
  cout << r << ", " << exp(delE/(kB*T)) << endl;
  if(r <= exp(delE/(kB*T))) {
    return 1;
  }
  else  {
    return 0;
  }
}

int main(int argc, char *argv[]) { // Main function
  int L = atoi(argv[1]);
  int n = atoi(argv[2]);

  // Initialize lattice
  arma::Mat<int> lattice(L,L,arma::fill::zeros);
  initLattice(lattice, L);



  //----------------------------------------------------------------------------
  // Initiate Monte Carlo Markov Chain
  int Ei = 8;
  int Ej = 0;
  int delE = 0;
  random_device rd;
  mt19937_64 generator(rd());
  uniform_int_distribution<int> distribution(0,L-1);

  fstream outfile;
  outfile.open("../../energy.txt", fstream::out);
  outfile << "Energy" << endl;

  for(int i = 0; i < n; i++) {
    // Choose and flip spin
    int k = distribution(generator);
    int l = distribution(generator);
    Ej = localE(lattice, L, k, l); // Energy pre-flip
    flip(lattice(k,l)); // Flipping
    Ei = localE(lattice, L, k, l); // Energy post-flip
    delE = Ei-Ej;
    //Metropolis(delE(k,l));
    bool temp = Metropolis(delE);
    cout << "Metropolis: " << temp << endl;
    // Update averages
    outfile << exp(delE/(kB*1.0)) << endl;
    //lattice.print();
  }
  cout << "Test" << endl;

  /// THINGS WE NEED
  ///
  /// mean energy E
  /// mean magnetization |M|
  /// specific heat Cv
  /// susceptibility chi (X)

  outfile.close();

}
