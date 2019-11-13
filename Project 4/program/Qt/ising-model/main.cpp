/*
  Last edited: 01.11.2019 14:50 by Erlend T. North
*/

#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>
#include <chrono>
//#include <mpi.h>
#define kB 1
using namespace std;
namespace ch = std::chrono;

// Functions
void initLattice(arma::Mat<int> &M,int n,int L) {
  random_device rd;
  mt19937_64 generator(rd());
  uniform_int_distribution<int> distribution(1,L);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      M(i,j) = distribution(generator);
      cout << M(i,j) << " "; // For checking the distribution
    }
    cout << endl; // For checking the distribution
  }
}

void flip(int &spin,int L) { // Generates a random spin-value
  random_device rd;
  mt19937_64 generator(rd());
  uniform_int_distribution<int> distribution(1,L);
  spin = distribution(generator);
}

double delE(int k,int l) {
  // Do Stuff
  return 0;
}

void MonteCarlo(int n, double a, double b, double  &integral, double  &std) {
  // Do stuff
}

int Metropolis(double delE,double T=1.0) {
  random_device rd;
  mt19937_64 generator(rd());
  uniform_real_distribution<double> distribution(0.0,1.0);

  double r = distribution(generator);
  if(r <= exp(delE/(kB*T))) {
  return 1;
  }
  return 0;
}

int main(int argc, char *argv[]) { // Main function
  int n = atoi(argv[1]);
  int L = atoi(argv[2]);
  int MC_cycles = atoi(argv[3]);

  // Initialize lattice
  arma::Mat<int> lattice(n,n,arma::fill::zeros);
  initLattice(lattice, n, L);

  // Choose and flip spin




  //----------------------------------------------------------------------------
  // Initiate Monte Carlo Markov Chain
  int Ei = 8;
  int Ej = 0;
  random_device rd;
  mt19937_64 generator(rd());
  uniform_int_distribution<int> distribution(0,n-1);
  for(int i = 0; i < MC_cycles; i++) {
    int k = distribution(generator);
    int l = distribution(generator);
    flip(lattice(k,l), L);
    Ej = 0; // CALCULATE
    Metropolis(delE(k,l));
    // Update averages
  }
  cout << "Test" << endl;

}
