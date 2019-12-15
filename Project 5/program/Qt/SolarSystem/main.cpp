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

void velocityVerlet(Body &body, vector<Body> system, int i, double dt, double dt_pos, double dt_vel) {
    body.acceleration(system, i); // Cannot re-use new_acc as other bodies have moved
    arma::rowvec acc_old = body.get_acc();
    body.new_pos(body.get_pos().row(i) + body.get_vel()*dt + acc_old*dt_pos, i+1);
    // Once we have more dynamic bodies this will be wrong(?). Will have to split vVerlet in two parts. acc_old and new_pos, then acc_new and new_vel
    body.acceleration(system, i+1);
    body.new_vel(body.get_vel() + (acc_old+body.get_acc())*dt_vel);
    //body.new_vel(body.get_vel() + )

}


int main() {
    ifstream infile("../../input.txt");

    int n, yr;
    double dt, dt_pos, dt_vel;

    string trash, planet;
    vector<string> planetList;
    if(infile.is_open()) {
        infile >> trash >> trash >> trash; // "n" >> " , " >> "yr"
        infile >> n >> trash >> yr; // 100 >> " , " >> 100
        infile >> trash; // "planets"
        while(getline(infile, planet))  {
            planetList.push_back(planet);
        }
        planetList.erase(planetList.begin());

    }
    else {
        cout << "Unable to open file" << endl;
        cout << "Quitting program..." << endl;
        return 1;
    }

    dt = (1.0*(double)yr)/(double)n;
    dt_pos = (dt*dt)/2.0;
    dt_vel = dt/2.0;

    arma::vec t(n, arma::fill::zeros);

    vector<Body> system;
    for(int i = 0; i < planetList.size(); i++) {
        system.push_back(Body(planetList[i], n));
    }


    //Body Sun("Sun", n);
    //Body Earth("Earth", n);
    //vector<Body> system = {Sun, Earth};


    ofstream data;

    // Main loop
    for(int bodyCount = 0; bodyCount < system.size(); bodyCount++) {
        Body body = system[bodyCount];
        data.open("../../forwardEulerbodyOutput/" + body.get_name() + ".txt");
        data << "x , y , z" << endl;
        data << body.get_pos()(0,0) << " , " << body.get_pos()(0,1) << " , " << body.get_pos()(0,2) << endl;
        if(bodyCount==0) {
            for(int i = 0; i < n-1; i++) { // Makes the sun static.. for now
                body.new_pos(body.get_pos().row(i), i+1);
                data << body.get_pos()(i+1,0) << " , " << body.get_pos()(i+1,1) << " , " << body.get_pos()(i+1,2) << endl;
            }
        }
        else {
            for(int i = 0; i < n-1; i++) {
                forwardEuler(body, system, i, dt);
                data << body.get_pos()(i+1,0) << " , " << body.get_pos()(i+1,1) << " , " << body.get_pos()(i+1,2) << endl;

            }
        }

    }

    /*
    data.open("../../forwardEulerbodyOutput/" + Earth.get_name() + ".txt");
    data << "x , y , z" << endl;
    data << Earth.get_pos()(0,0) << " , " << Earth.get_pos()(0,1) << " , " << Earth.get_pos()(0,2) << endl;

    // Forward Euler

    for(int i = 0; i < n-1; i++) {
        forwardEuler(Earth, system, i, dt);
        //cout << "Iteration: " << i << " complete!" << endl;
        data << Earth.get_pos()(i+1,0) << " , " << Earth.get_pos()(i+1,1) << " , " << Earth.get_pos()(i+1,2) << endl;

    }
    data.close();
    cout << "Forward Euler complete!" << endl;

    data.open("../../velocityVerletbodyOutput/" + Earth.get_name() + ".txt");
    data << "x , y , z" << endl;
    data << Earth.get_pos()(0,0) << " , " << Earth.get_pos()(0,1) << " , " << Earth.get_pos()(0,2) << endl;


    // Velocity Verlet

    for(int i = 0; i < n-1; i++) {
        velocityVerlet(Earth, system, i, dt, dt_pos, dt_vel);

        //cout << "Iteration: " << i << " complete!" << endl;
        data << Earth.get_pos()(i+1,0) << " , " << Earth.get_pos()(i+1,1) << " , " << Earth.get_pos()(i+1,2) << endl;

    }
    data.close();
    cout << "Velocity Verlet complete!" << endl;
    */

    //system("cd ..");
    //system("cd ..");
    //system("python ../../plot.py");

    return 0;
}
