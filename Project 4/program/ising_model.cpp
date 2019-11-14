/*
  Last edited: 01.11.2019 14:50 by Erlend T. North
*/

#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>
#include <chrono>
#include <mpi.h>
#define kB=1
using namespace std;
namespace ch = std::chrono;

// Functions
void flip(&spin, L) { // Generates
  random_device rd;
  mt9937_64 generator(rd());
  uniform_int_distribution<double> distribution(1,L);
  spin = distribution(generator)
}

double delE(k,l) {

}

void MonteCarlo(int n, double a, double b, double  &integral, double  &std) {
  // Do stuff
}

int Metropolis(delE, T=1.0) {
  random_device rd;
  mt9937_64 generator(rd());
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

  random_device rd;
  mt9937_64 generator(rd());
  uniform_int_distribution<double> distribution(1,L);

  // Initialize lattice
  arma::Mat<int> lattice(n,n,arma::fill::zeros);
  for(int i = 0; i < n; i++) {
    for(int j = 0; i < n; i++) {
      lattice(i,j) = distribution(generator)
      cout << lattice(i,j) << " "; // For checking the distribution
    }
    cout << endl; // For checking the distribution
  }

  // Choose and flip spin





  //----------------------------------------------------------------------------
  // Initiate Monte Carlo Markov Chain
  int Ei = 8;
  int Ej = 0;
  uniform_int_distribution<double> distribution(1,n);
  for(int i = 0; i < MC_cycles; i++) {
    int k = distribution(generator);
    int l = distribution(generator);
    flip(lattice(k,l));
    Ej = CALCULATE
    Metropolis(delE(k,l));
    // Update averages
  }

  for(int i = 0; i < n; i++) {
    // Flip random spin
    // Compute delta E and perform Metropolis test
  }

}
