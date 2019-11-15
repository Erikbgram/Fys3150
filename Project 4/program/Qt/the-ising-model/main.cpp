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

#define kB 1
#define J 1

using namespace std;
ofstream outfile;

random_device rd;
mt19937_64 mt_gen(rd());
uniform_real_distribution<double> rand_frac(0,1);

int periodic(int i, int limit, int add) { // Periodic Boundary Conditions
    return (i+limit+add) % (limit);
}

void initialize(int L, double temp, arma::Mat<int> lattice, double& E, double& M) {  // Initialize energy and magnetization
    // Setup lattice and initial magnetization
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            if(temp < 1.5) {
                lattice(x,y) = 1;
            }
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

void Metropolis(int L, arma::Mat<int> lattice, double& E, double& M, double *w) { // The Metropolis algorithm
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
            if(rand_frac(mt_gen) <= w[deltaE+8]) {
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
    outfile << setiosflags(ios::showpoint | ios::uppercase);
    outfile << setw(15) << setprecision(8) << temperature;
    outfile << setw(15) << setprecision(8) << Eaverage/L/L;
    outfile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    // outfile << setw(15) << setprecision(8) << Maverage/L/L;
    outfile << setw(15) << setprecision(8) << Mabsvariance/temperature;
    outfile << setw(15) << setprecision(8) << Mabsaverage/L/L << endl;
}

int main(int argc, char *argv[]) { // Main function
    char *outfilename;
    int **mcs;
    double w[17], average[5], E, M;
    int L = atoi(argv[2]);
    int n = atoi(argv[3]);
    double initial_temp = atof(argv[4]);
    double final_temp = atof(argv[5]);
    double temp_step = atof(argv[6]);
    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ) {
        cout << "Bad Usage: " << argv[0] <<
        " read also output file on same line" << endl;
        exit(1);
    }
    else {
        char* outfilename = argv[1];
    }

    // Read in initial values such as size of lattice, temp and cycles
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

        // Start Monte Carlo Computation
        for(int cycles = 1; cycles <= n; cycles++) {
            Metropolis(L, lattice, E, M, w);

            // Update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }
        // Print results
        outfile.open(outfilename);
        output(L, n, temp, average);
        outfile.close();
    }
    outfile.close(); // Close output file
    return 0;
}
