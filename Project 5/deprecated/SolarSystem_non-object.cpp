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

const double GMs = 4*M_PI*M_PI;

using namespace std;


// Functions


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

    arma::mat r_s(n, 3, arma::fill::zeros);
    arma::mat r_e(n, 3, arma::fill::zeros);
    arma::mat v_e(n, 3, arma::fill::zeros);

    /*
    // Initial data 01.01.2020 00:00
    r_e(0,0) = -1.701443808233178E-01; // AU
    r_e(0,1) = 9.765603186945023E-01;
    r_e(0,2) = -1.822571331401604E-05;
    v_e(0,0) = -1.724754632489101E-02*365.242199; // AU/day -> AU/year
    v_e(0,1) = -2.983511998041463E-03*365.242199;
    v_e(0,2) = 6.558216388187520E-07*365.242199;
    */


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

        /*
        This is Euler
        v_e(i,0) = v_e(i-1,0) + a*(r_e(i-1,0)-r_s(0,0)) *dt;
        v_e(i,1) = v_e(i-1,1) + a*(r_e(i-1,1)-r_s(0,1)) *dt;
        v_e(i,2) = v_e(i-1,2) + a*(r_e(i-1,2)-r_s(0,2)) *dt;
        r_e(i,0) = r_e(i-1,0) + v_e(i-1,0)*dt;
        r_e(i,1) = r_e(i-1,1) + v_e(i-1,1)*dt;
        r_e(i,2) = r_e(i-1,2) + v_e(i-1,2)*dt;
        */

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
