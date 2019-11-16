#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>

std::random_device rd;
std::mt19937_64 mt_gen1(rd());
std::uniform_real_distribution<double> rand_frac1(0,1);

long long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count(); //random_device did not work (always made the same sequence)
std::mt19937_64 mt_gen2(seed);
std::uniform_real_distribution<double> rand_frac2(0,1);

int main(int argc, char const *argv[]) {

    //std::cout << rand_frac2(mt_gen2(seed)) << "\n" ;

    //rand_frac2(mt_gen2(seed)) ;










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
