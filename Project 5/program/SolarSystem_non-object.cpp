/*
    Last edited: 13.12.2019 13:47 by Erlend T. North
*/

// THIS CODE IS NOT FUNCTIONAL. IT IS A PLACEHOLDER FOR NON-ORIENTED STUFF, WHICH IS CURRENTLY COMMENTED OUT.
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

void deprecated_acceleration(Body body, Body *BodyList, int i) { // Sums over acceleration from all other bodies to specific body
    arma::rowvec vec;
    for(int j = 0; j < 2; j++) { // The number here is the amount of bodies in BodyList
        if(BodyList[j].get_name() != body.get_name()) {
            vec = BodyList[j].get_pos().row(i)-body.get_pos().row(i);
            body.new_acc(body.get_acc() - (GMs*BodyList[j].get_mass()) / (pow(mod(vec), 3))*vec);
        }
    }
}

void forwardEuler(Body body, arma::Row<Body> bodyList, int i, double dt) { // Performs 1 iteration of Forward Euler
    body.acceleration(bodyList, i);
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

    /*
    // Initial data 01.01.2020 00:00
    r_e(0,0) = -1.701443808233178E-01; // AU
    r_e(0,1) = 9.765603186945023E-01;
    r_e(0,2) = -1.822571331401604E-05;
    v_e(0,0) = -1.724754632489101E-02*365.242199; // AU/day -> AU/year
    v_e(0,1) = -2.983511998041463E-03*365.242199;
    v_e(0,2) = 6.558216388187520E-07*365.242199;
    */

    cout << "hey" << endl;

    Body Sun("Sun", n);
    Body Earth("Earth", n);
    arma::Row<Body> bodyList = {Sun, Earth};
    cout << bodyList[0].get_name() << endl;
    /*
    Body Mercury("Mercury", n);
    Body Venus("Venus", n);
    Body Mars("Mars", n);
    Body Jupiter("Jupiter", n);
    Body Saturn("Saturn", n);
    Body Uranus("Uranus", n);
    Body Neptune("Neptune", n);
    Body Pluto("Pluto", n);
    Body* BodyList[10]{&Sun, &Mercury, &Venus, &Earth, &Mars, &Jupiter, &Saturn, &Uranus, &Neptune, &Pluto};
    */
    //Body** List = new Body*[2];
    //vector<Body> BodyList{Sun, Earth};

    //Body BodyList[2]{Sun, Earth};


    /*
    Body Sun("Sun", 1.0, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Mercury("Mercury", 1.0/6060606.06, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Venus("Venus", 1.0/408163.27, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Earth("Earth", 1.0/333333.33, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Mars("Mars", 1.0/3030303.03, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Jupiter("Jupiter", 1.0/1052.63, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Saturn("Saturn", 1.0/3636.36, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Uranus("Uranus", 1.0/22727.27, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Neptune("Neptune", 1.0/19417.48, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    Body Pluto("Pluto", 1.0/152671755.7, arma::mat(n,3, arma::fill::zeros), arma::mat(n,3, arma::fill::zeros));
    */

    Earth.acceleration(bodyList, 0);
    Earth.get_acc().print();

    /*
    ofstream outfile;
    outfile.open("../../Euleroutput.txt");
    outfile << "t, rx, ry, rz" << endl;
    outfile << setprecision(8) << t[0] << ", " << r_e(0,0) << ", " << r_e(0,1) << ", " << r_e(0,2) << endl;
    //cout << v_e(0,0) << ", " << v_e(0,1) << ", " << v_e(0,2) << endl;

    double a;


    // Forward Euler

    for(int i = 1; i < n; i++) {
        a = -(4*pow(M_PI,2)) / pow(sqrt(pow((r_e(i-1,0)-r_s(0,0)),2) + pow((r_e(i-1,1)-r_s(0,1)),2) + pow((r_e(i-1,2)-r_s(0,2)),2)),3);
        v_e(i) = v_e(i-1) + a*(r_e(i-1)-r_s(0)) *dt;
        r_e(i) = r_e(i-1) + v_e(i-1)*dt;
        //cout << v_e(i,0) << ", " << v_e(i,1) << ", " << v_e(i,2) << endl;
        outfile << setprecision(8) << t[i] << ", " << r_e(i,0) << ", " << r_e(i,1) << ", " << r_e(i,2) << endl;
    }

    forwardEuler(Earth, bodyList, 0, dt);
    Earth.get_pos().print();
    outfile.close();

    cout << "Euler done!" << endl;

    */

    /*
    outfile.open("../../VVerletoutput.txt");
    outfile << "t, rx, ry, rz" << endl;
    outfile << setprecision(8) << t[0] << ", " << r_e(0,0) << ", " << r_e(0,1) << ", " << r_e(0,2) << endl;
    //cout << v_e(0,0) << ", " << v_e(0,1) << ", " << v_e(0,2) << endl;

    double a_scalar;
    arma::vec a_old(3);
    arma::vec a_new(3);

    arma::mat r_s(n, 3, arma::fill::zeros);
    arma::mat r_e(n, 3, arma::fill::zeros);
    arma::mat v_e(n, 3, arma::fill::zeros);

    // Velocity Verlet
    for(int i = 1; i < n; i++) {
        a_scalar = -(4*pow(M_PI,2)) / pow(sqrt(pow((r_e(i-1,0)-r_s(0,0)),2) + pow((r_e(i-1,1)-r_s(0,1)),2) + pow((r_e(i-1,2)-r_s(0,2)),2)),3);
        a_old = a_scalar*(r_e(i-1)-r_s(0));
        r_e(i,0) = r_e(i-1,0) + v_e(i-1,0)*dt + (dt*dt)/2*a_old(0);
        r_e(i,1) = r_e(i-1,1) + v_e(i-1,1)*dt + (dt*dt)/2*a_old(1);
        r_e(i,2) = r_e(i-1,2) + v_e(i-1,2)*dt + (dt*dt)/2*a_old(2);
        a_scalar = -(4*pow(M_PI,2)) / pow(sqrt(pow((r_e(i,0)-r_s(0,0)),2) + pow((r_e(i,1)-r_s(0,1)),2) + pow((r_e(i,2)-r_s(0,2)),2)),3);
        a_new = a_scalar*(r_e(i)-r_s(0));
        v_e(i,0) = v_e(i-1,0) + (a_new(0)+a_old(0))*(dt/2);
        v_e(i,1) = v_e(i-1,1) + (a_new(1)+a_old(1))*(dt/2);
        v_e(i,2) = v_e(i-1,2) + (a_new(2)+a_old(2))*(dt/2);


        v_e(i,0) = v_e(i-1,0) + a*(r_e(i-1,0)-r_s(0,0)) *dt;
        v_e(i,1) = v_e(i-1,1) + a*(r_e(i-1,1)-r_s(0,1)) *dt;
        v_e(i,2) = v_e(i-1,2) + a*(r_e(i-1,2)-r_s(0,2)) *dt;
        r_e(i,0) = r_e(i-1,0) + v_e(i-1,0)*dt;
        r_e(i,1) = r_e(i-1,1) + v_e(i-1,1)*dt;
        r_e(i,2) = r_e(i-1,2) + v_e(i-1,2)*dt;

        //cout << v_e(i,0) << ", " << v_e(i,1) << ", " << v_e(i,2) << endl;
        outfile << setprecision(8) << t[i] << ", " << r_e(i,0) << ", " << r_e(i,1) << ", " << r_e(i,2) << endl;
    }
    outfile.close();
    */


    //v_e.print("v");
    //r_e.print("r");

    //system("cd ..");
    //system("cd ..");
    //system("python ../../plot.py");

    return 0;
}
