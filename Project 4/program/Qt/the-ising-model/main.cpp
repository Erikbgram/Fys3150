/*
    Last edited: 15.11.2019 18:54 by Erlend T. North
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <random>
#include <armadillo>
#include <iomanip>

// kB = 1, J = 1

using namespace std;
namespace ch = std::chrono;

ofstream outfile;

long long seed = ch::high_resolution_clock::now().time_since_epoch().count(); //random_device did not work (always made the same sequence)
mt19937_64 mt_gen(seed);
uniform_real_distribution<double> rand_frac(0,1);

int periodic(int i, int limit, int add) { // Periodic Boundary Conditions
    return (i+limit+add) % (limit);
}

void initialize(int L, double temp, arma::Mat<int> &lattice, double& E, double& M) {  // Initialize energy and magnetization
    // Setup lattice and initial magnetization
    double spin;
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {

            spin = rand_frac(mt_gen);
            if(spin<0.5) {
                lattice(x,y) = -1;
            }
            else {
                lattice(x,y) = 1;
            }

            /*
            if(temp < 1.5) {
                lattice(x,y) = 1; // Spin orientation for the ground state
            }
            */

            M += lattice(x,y);
        }
    }

    // Setup initial energy
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            E -= lattice(x,y) * (
                        lattice(x, periodic(y,L,-1)) +
                        lattice(periodic(x,L,-1), y)
                        );
        }
    }
}

void Metropolis(int L, arma::Mat<int> &lattice, double& E, double& M, double *w) { // The Metropolis algorithm
    // Loop over all spins
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            //Find random position
            int ix = rand_frac(mt_gen)*L; // If problem put "(int)" in front of rand_frac
            int iy = rand_frac(mt_gen)*L;
            int deltaE = 2*lattice(ix,iy) * (
                        lattice(ix,periodic(iy,L,-1)) +
                        lattice(periodic(ix,L,-1),iy) +
                        lattice(ix,periodic(iy,L, 1)) +
                        lattice(periodic(ix,L, 1),iy)
                        );

            // Here we perform the Metropolis test

            double r = rand_frac(mt_gen);

            double E_ = w[deltaE+8];

            if(r <= E_) {
                lattice(ix,iy) *= -1; // Flip one spin and accept new spin config

                // Update energy and magnetization
                M += 2*lattice(ix,iy);
                E += deltaE;
            }
        }
    }
}

void output(int L, int n, double temperature, double *average) { // Prints to file the results of the calculations
    double norm = 1.0/n; // Divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;

    // all expectation values are per spin, divide by 1/L/L
    double Evariance = (E2average- Eaverage*Eaverage)/L/L;
    double Mvariance = (M2average - Maverage*Maverage)/L/L;
    double Mabsvariance = (M2average - Mabsaverage*Mabsaverage)/L/L;
    outfile << temperature << ", ";
    outfile << Eaverage/L/L << ", ";
    outfile << Evariance/temperature/temperature << ", ";
    // outfile << setw(15) << setprecision(8) << Maverage/L/L;
    outfile << Mabsvariance/temperature << ", ";
    outfile << Mabsaverage/L/L << endl;
}

int main(int argc, char *argv[]) { // Main function
    // Read in initial values
    double w[17], average[5], E, M, E2, M2;
    int L = atoi(argv[1]);
    int n = atoi(argv[2]);
    double initial_temp = atof(argv[3]);
    double final_temp = atof(argv[4]);
    double temp_step = atof(argv[5]);

    // Initialize lattice
    arma::Mat<int> lattice(L,L,arma::fill::zeros);

    for(double temp = initial_temp; temp <= final_temp; temp += temp_step) {
        // Initialize energy and magnetization
        E = M = 0;

        // setup array for possible energy changes
        for(int de = -8; de <= 8; de++) {
            w[de+8] = 0;
        }
        for(int de = -8; de <= 8; de+=4) {
            w[de+8] = exp(-de/temp);
        }

        // Initialize array for expectation values
        for(int i = 0; i < 5; i++) {
            average[i] = 0;
        }
        initialize(L, temp, lattice, E, M);
        lattice.print();

        // Start Monte Carlo Computation
        outfile.open("../../output.txt");
        outfile << "final_temp , Eaverage/L/L , Evariance/final_temp/final_temp , Mabsvariance/final_temp , Mabsvariance/L/L1" << endl;
        for(int cycles = 1; cycles <= n; cycles++) {
            cout << "cycle = " << cycles << endl;
            Metropolis(L, lattice, E, M, w);

            // Print to console
            cout << "E = " << E/L/L << ", M = " << M/L/L << endl;
            cout << endl;

            // Update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
            output(L, n, temp, average);
        }

        // Print results

        outfile.close();
        lattice.print();
    }
    return 0;
}
