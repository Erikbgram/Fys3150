/*
    Last edited: 12.12.2019 14:05 by Erlend T. North
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <armadillo>
#include <cmath>
#include <body.h>
#include <vector>

using namespace std;

// Functions
double mod(arma::rowvec vec) { // Calculates the modulus of a vector
    double val = 0;
    for(int i = 0; i < vec.n_elem; i++) {
        val += vec(i)*vec(i);
    }
    return sqrt(val);
}


void forwardEuler(Body &body, vector<Body> system, int i, double dt) { // Performs 1 iteration of Forward Euler
    body.acceleration(system, i);
    body.new_vel(body.get_vel() + body.get_acc()*dt);
    body.new_pos(body.get_pos().row(i) + body.get_vel()*dt, i+1);
    //body.write_data(i);
}


int main() {
    ifstream infile("../../input.txt");

    int n, yr;
    double dt;

    string trash;
    if(infile.is_open()) {
        infile >> trash >> trash >> trash;
        infile >> n >> trash >> yr;
    }
    else {
        cout << "Unable to open file" << endl;
        n = 100;
        yr = 100;
    }

    dt = (1.0*(double)yr)/(double)n;

    double dtv = dt/2.0;
    double dtr = (dt*dt)/2.0;
    arma::vec t(n, arma::fill::zeros);

    Body Sun("Sun", n);
    Body Earth("Earth", n);
    vector<Body> system = {Sun, Earth};

    ofstream data;
    data.open("../../bodyOutput/" + Earth.get_name() + ".txt", fstream::app);
    data << Earth.get_pos()(0,0) << " , " << Earth.get_pos()(0,1) << " , " << Earth.get_pos()(0,2) << endl;

    for(int i = 0; i < n-1; i++) { /// CHANGE THIS BACK TO n-1
        forwardEuler(Earth, system, i, dt);
        cout << "Iteration: " << i << " complete!" << endl;
        data << Earth.get_pos()(i+1,0) << " , " << Earth.get_pos()(i+1,1) << " , " << Earth.get_pos()(i+1,2) << endl;

    }
    data.close();



    //system("cd ..");
    //system("cd ..");
    //system("python ../../plot.py");

    return 0;
}
