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

const double GMs = 4*M_PI*M_PI;

using namespace std;

double acceleration(Body body1, Body body2, double acc = 0) {
    acc += (GMs*body1.get_mass()) / (pow( abs(body1.get_pos()(0) - body2.get_pos()(0)) ,2));
    return acc;
}

void ForwardEuler(Body* BodyList, Body body, int n) {
    for(int i = 0; i < n-1; i++) {
        //vel(i+1) = vel(i) + acceleration();
    }
}

int main() {

    ifstream infile("../../input.txt");

    int n, yr;
    double dt;

    string trash;
    if(infile.is_open()) {
        infile >> trash >> trash >> trash;
        cout << trash;
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

    arma::mat r_s(n, 3, arma::fill::zeros);
    arma::mat r_e(n, 3, arma::fill::zeros);
    arma::mat v_e(n, 3, arma::fill::zeros);

    // Initial data 01.01.2020 00:00

    /*
    r_e(0,0) = -1.701443808233178E-01; // AU
    r_e(0,1) = 9.765603186945023E-01;
    r_e(0,2) = -1.822571331401604E-05;
    v_e(0,0) = -1.724754632489101E-02*365.242199; // AU/day -> AU/year
    v_e(0,1) = -2.983511998041463E-03*365.242199;
    v_e(0,2) = 6.558216388187520E-07*365.242199;
    */
    r_e(0,0) = 1; // AU
    r_e(0,1) = 0;
    r_e(0,2) = 0;
    v_e(0,0) = 0; // AU/day -> AU/year
    v_e(0,1) = 2*M_PI;
    v_e(0,2) = 0;

    cout << "hey" << endl;

    Body Sun("Sun", n);
    Body Earth("Earth", n);
    Body Mercury("Mercury", n);
    Body Venus("Venus", n);
    Body Mars("Mars", n);
    Body Jupiter("Jupiter", n);
    Body Saturn("Saturn", n);
    Body Uranus("Uranus", n);
    Body Neptune("Neptune", n);
    Body Pluto("Pluto", n);


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

    Body* BodyList[10]{&Sun, &Mercury, &Earth, &Mars, &Jupiter, &Saturn, &Uranus, &Neptune, &Pluto};
    //Earth.init(r_e, v_e);

    cout << acceleration(Sun, Earth) << endl;

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
    outfile.close();

    cout << "Euler done!" << endl;

    outfile.open("../../VVerletoutput.txt");
    outfile << "t, rx, ry, rz" << endl;
    outfile << setprecision(8) << t[0] << ", " << r_e(0,0) << ", " << r_e(0,1) << ", " << r_e(0,2) << endl;
    //cout << v_e(0,0) << ", " << v_e(0,1) << ", " << v_e(0,2) << endl;

    double a_scalar;
    arma::vec a_old(3);
    arma::vec a_new(3);

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



    //v_e.print("v");
    //r_e.print("r");

    //system("cd ..");
    //system("cd ..");
    //system("python ../../plot.py");

    return 0;
}
