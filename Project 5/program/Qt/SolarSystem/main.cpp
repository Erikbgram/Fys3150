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

const double GMs = 4*M_PI*M_PI;

using namespace std;


// Functions
double mod(arma::rowvec vec) { // Calculates the modulus of a vector
    double val = 0;
    for(int i = 0; i < vec.n_elem; i++) {
        val += vec(i)*vec(i);
    }
    return sqrt(val);
}


void forwardEuler(Body body, vector<Body> system, int i, double dt) { // Performs 1 iteration of Forward Euler
    body.acceleration(system, i);
    body.new_vel(body.get_vel() + body.get_acc()*dt);
    body.new_pos(body.get_pos()(i) + body.get_vel()*dt, i+1);
    body.write_data(i);
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

    cout << "hey" << endl;

    Body Sun("Sun", n);
    Body Earth("Earth", n);
    vector<Body> system = {Sun, Earth};
    cout << system[0].get_name() << endl;

    Earth.acceleration(system, 0);
    Earth.get_acc().print();



    //system("cd ..");
    //system("cd ..");
    //system("python ../../plot.py");

    return 0;
}
