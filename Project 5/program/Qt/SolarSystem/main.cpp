/*
    Last edited: 11.12.2019 14:22 by Erlend T. North
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <armadillo>
#include <cmath>

using namespace std;

double force(double m1, double m2, double r1, double r2) {
    double g = 1;
    return (g*m1*m2)/(fabs(r2-r1));
}


int main() {
    ifstream infile("../../input.txt");

    int n, years;
    double dt;

    string line;
    if(infile.is_open()) {
        infile >> n;
        infile >> years;
    }
    else {
        cout << "Unable to open file";
        n = 100;
        years = 100;
    }

    dt = (1.0*(double)years)/(double)n;

    arma::vec t(n); t[0]=0;
    for (int i = 1; i < n; i++) {
        t[i] = i*dt;
    }

    arma::mat r_s(1, 3, arma::fill::zeros);
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



    ofstream outfile;
    outfile.open("../../Euleroutput.txt");
    outfile << "t, rx, ry, rz" << endl;
    outfile << setprecision(8) << t[0] << ", " << r_e(0,0) << ", " << r_e(0,1) << ", " << r_e(0,2) << endl;
    //cout << v_e(0,0) << ", " << v_e(0,1) << ", " << v_e(0,2) << endl;

    double a;

    // Forward Euler
    for(int i = 1; i < n; i++) {
        a = -(4*pow(M_PI,2)) / pow(sqrt(pow((r_e(i-1,0)-r_s(0,0)),2) + pow((r_e(i-1,1)-r_s(0,1)),2) + pow((r_e(i-1,2)-r_s(0,2)),2)),3);
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



    outfile.open("../../VVerletoutput.txt");
    outfile << "t, rx, ry, rz" << endl;
    outfile << setprecision(8) << t[0] << ", " << r_e(0,0) << ", " << r_e(0,1) << ", " << r_e(0,2) << endl;
    //cout << v_e(0,0) << ", " << v_e(0,1) << ", " << v_e(0,2) << endl;


    // Velocity Verlet
    for(int i = 1; i < n; i++) {
        a = -(4*pow(M_PI,2)) / pow(sqrt(pow((r_e(i-1,0)-r_s(0,0)),2) + pow((r_e(i-1,1)-r_s(0,1)),2) + pow((r_e(i-1,2)-r_s(0,2)),2)),3);
        r_e(i,0) = r_e(i-1,0) + v_e(i-1,0)*dt + (dt*dt)/2;

        /*
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
