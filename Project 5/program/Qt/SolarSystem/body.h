/*
    Last edited: 12.12.2019 14:05 by Erlend T. North
*/

#include <iostream>
#include <armadillo>
#include <string>

using namespace std;

class Body {
    // Attributes
    string name;
    double mass;
    arma::mat pos;
    arma::rowvec vel;
    arma::rowvec acc;

public:
    // Methods
    Body(string new_name, int n) {
        name = new_name;
        pos = arma::mat(n,3,arma::fill::zeros);
        vel = arma::rowvec(3,arma::fill::zeros);
        acc = arma::rowvec(3,arma::fill::zeros);
        ifstream infile("../../BodyData/" + name + ".txt");

        string trash;
        if(infile.is_open()) {
            infile >> name >> trash >> mass;
            infile >> pos(0,0) >> trash >> pos(0,1) >> trash >> pos(0,2);
            infile >> vel(0) >> trash >> vel(1) >> trash >> vel(2);
        }
        else {
            cout << "Unable to open file" << endl;
        }
    }

    string get_name() {
        return name;
    }

    double get_mass() {
        return mass;
    }

    arma::mat get_pos() {
        return pos;
    }

    arma::mat get_vel() {
        return vel;
    }

    arma::mat get_acc() {
        return acc;
    }

    void new_pos(arma::rowvec new_pos, int i) {
        pos.row(i) = new_pos;
    }

    void new_vel(arma::rowvec new_vel) {
        vel = new_vel;
    }

    void new_acc(arma::rowvec new_acc) {
        acc = new_acc;
    }

};
