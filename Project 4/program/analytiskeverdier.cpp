// Sist endret: 20.11.2019 12:03 av Alexandra Jahr Kolstad

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <cmath>

/*
std::random_device rd;
std::mt19937_64 mt_gen1(rd());
std::uniform_real_distribution<double> rand_frac1(0,1);

long long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count(); //random_device did not work (always made the same sequence)
std::mt19937_64 mt_gen2(seed);
std::uniform_real_distribution<double> rand_frac2(0,1);
*/

double beta(int kB, int T){
    return 1 / (kB * T) ;
}

double Z(int J, int kB, int T){
    return 4 * (cosh(8 * beta(kB, T) * J ) + 3 ) ;
}

double energy1(int J, int kB, int T) {
    return - (32*J/Z(J, kB, T)) * sinh(8 * beta(kB, T) * J) ;
}

double energy2 (int J, int kB, int T) {
    return (256*J*J / Z(J, kB, T)) * cosh(8 * beta(kB, T) * J) ;
}

double absmagnetization1 (int J, int kB, int T){
    return (8 / Z(J, kB, T)) * (exp(8 * beta(kB, T) * J) + 2) ;
}

double magnetization2(int J, int kB, int T) {
    return 32 / (Z(J, kB, T)) * (exp(8 * beta(kB, T) * J) + 1) ;
}

double heatcapacity(int J, int kB, int T) {
    return 1/(kB * T*T) * (energy2(J, kB, T) - energy1(J, kB, T) * energy1(J, kB, T)) ;
}

double susceptibility(int J, int kB, int T, int magnetization1) {
    return 1/(kB * T) * (magnetization2(J, kB, T) - magnetization1*magnetization1) ;
}

int main(int argc, char const *argv[]) {

    //std::cout << rand_frac2(mt_gen2(seed)) << "\n" ;

    //int random = rand_frac2(mt_gen2(seed)) ;

    int magnetization1 = 0 ;

    int kB = atoi(argv[1]) ;
    int T = atoi(argv[2]) ;
    int J = atoi(argv[3]) ;

    std::cout << "Energy1 = " << energy1(J, kB, T)/4 << "\n" ;
    std::cout << "Energy2 = " << energy2(J, kB, T)/4 << "\n" ;
    std::cout << "Heatcapacity = " << heatcapacity(J, kB, T)/4 << "\n" ;
    std::cout << "Magnetization1 = " << magnetization1/4 << "\n" ;
    std::cout << "Absmagnetization1 = " << absmagnetization1(J, kB, T)/4 << "\n" ;
    std::cout << "Magnetization2 = " << magnetization2(J, kB, T)/4 << "\n" ;
    std::cout << "Susceptibility = " << susceptibility(J, kB, T, magnetization1)/4 << "\n" ;




    /*
    char pointer[4] =  "Hei" ;
    char* p = pointer ;
    std::cout << "Navn:" << p << std::endl ;
    */
    /*
    std::ofstream out(argv[1]) ;
    out << argv[1] ;
    std::cout << out << std::endl ;
    */

    /*
    std::string outFile = "";
    if( argc == 10 ) {
        outFile = argv[1];
    }
    else {
      std::cout << "Usage: ./cppfile InputFile OutputFile\n";
      return 1;
    }

    std::cout << "Ny fil: " << outFile << std::endl ;
    */

    return 0;
}
