/*
Last changed: 08.09.2019 16:05 by Erlend
*/

#include <iostream> // input and output from command-line
#include <cmath> // math
#include <chrono> // timer
#include <string> // string
#include <fstream> // file
#include <algorithm> // std::max
#include <armadillo> // LU-decomposition

namespace ch = std::chrono;



//create dynamic array and fill with num
double* array(int n, double num=0) {
  double* arr = new double[n];
  for(int i = 0; i < n; i++) {
    arr[i] = num;
  }
  return arr;
}

//f(x)=100e^(-10x)
double* func(double* x, int n) {
  double* f = new double[n];
  double h = (x[n-1] - x[0]) / (n-1);
  for(int i = 0; i < n; i++) {
    f[i] = pow(h,2)*100*exp(-10*x[i]);
  }
  return f;
}

//closed-form solution, u(x)
double* cf_sol(double* x, int n) {
  double* u = new double[n];
  for(int i = 0; i < n; i++) {
    u[i] = 1 - (1 - exp(-10))*x[i] - exp(-10*x[i]);
  }
  return u;
}

//linspace for dynamic array
double* linspace(double start,double stop, int n) {
  double h = (stop - start)/(n-1);
  double* arr = new double[n];
  for(int i = 0; i < n; i++) {
    arr[i] = start + stop*i*h;
  }
  return arr;
}

int main() {

  //define n
  int n;
  std::cout << "Choose your n: ";
  std::cin >> n;

  //x array
  double* x = linspace(0,1,n);

  //step length
  double h = (x[n-1] - x[0]) / (n-1);
  std::cout << "h = " << h << std::endl;

  //defining array
  double *a = array(n-1, -1);
  double *d = array(n, 2);
  double *c = array(n-1, -1);
  double *v_gen = array(n);
  double *d_new = array(n); d_new[0] = d[0];
  double *b_tld = func(x, n);
  double *b_tld_new = array(n); b_tld_new[0] = b_tld[0];
  double *ad = array(n); //We can save 1 FLOP by pre-calculating a[i-1]/d_new[i-1]


  ch::steady_clock::time_point start = ch::steady_clock::now();

  // forward substitution
  for(int i = 1; i < n; i++) {
    ad[i-1] = a[i-1]/d_new[i-1]; // measure once, cut twice
    d_new[i] = d[i] - ad[i-1]*c[i-1];
    b_tld_new[i] = b_tld[i] - b_tld_new[i-1]*ad[i-1];
  }

  v_gen[n-1] = b_tld_new[n-1]/d_new[n-1]; // initial condition for v[-1] given by the backward subsitution

  // backward substitution
  for(int i = n-2; i > 0; i--) { //by setting i >= 0, we can fix the general-algorithms discontinuity in the 1st point, but we KNOW the value so can ignore it
    v_gen[i] = (b_tld_new[i] - c[i]*v_gen[i+1])/d_new[i];
  }

  ch::steady_clock::time_point stop = ch::steady_clock::now();

  ch::duration<double> time_span_general = ch::duration_cast<ch::nanoseconds>(stop - start);
  std::cout << "Time used for general algorithm = " << time_span_general.count()  << "s" << std::endl;

  /* Debugging snippet for seeing the difference of d_new from the general and special algorithms
  std::cout << "index, value\n";
  for(int i = 0; i < n; i++){
      std::cout << "    " << i << ", " << d_new[i] << std::endl;
  }*/

  //clean-up before starting the specialized algorithm
  delete[] d_new;
  delete[] b_tld_new;
  delete[] ad;



  // initialize arrays for specialized algorithm
  double *v_spl = array(n);
  d_new = array(n,2); d_new[0] = d[0];
  b_tld_new = array(n); b_tld_new[0] = b_tld[0];

  // filling d_new
  double j = 0;
  for(int i = 1; i < n; i++) {
    j = i;
    d_new[i] = (j + 1)/j;
  }

  /* Debugging snippet for seeing the difference of d_new from the general and special algorithms
  std::cout << "index, value\n";
  for(int i = 0; i < n; i++){
      std::cout << "    " << i << ", " << d_new[i] << std::endl;
  }*/

  //start clock
  start = ch::steady_clock::now();

  // forward substitution
  for(int i = 1; i < n; i++) {
    b_tld_new[i] = b_tld[i] + (b_tld_new[i-1])/(d_new[i-1]);
  }

  v_spl[n-1] = b_tld_new[n-1]/d_new[n-1]; // initial condition for v[-1] given by the backward subsitution

  //backward substitution
  for(int i = n-2; i > 0; i--) {
    v_spl[i] = (b_tld_new[i] + v_spl[i+1])/d_new[i];
  }

  stop = ch::steady_clock::now();

  ch::duration<double> time_span_special = ch::duration_cast<ch::nanoseconds>(stop - start);
  std::cout << "Time used for specialized algorithm = " << time_span_special.count()  << "s" << std::endl;

  //clean-up before the LU-decomposition
  delete[] a;
  delete[] d;
  delete[] c;
  delete[] d_new;
  delete[] b_tld_new;

  //exact solution
  double* u = cf_sol(x, n);

  //computing error of general algorithm
  double *eps = array(n);
  double max_eps = 0;
  for(int i = 0; i < n; i++) {
    eps[i] = std::log10( std::abs( (v_gen[i] - u[i]) /u[i] ) );
    if(i==0) {
      max_eps = eps[i];
    }
    else {
      max_eps = std::max(eps[i], max_eps);
    }

  }

  std::cout << "\n" "Printing eps: " "\n" ;
  for ( int i = 0 ; i < n ; i++) {
    std::cout << eps[i] << "\n" ;
  }

  delete[] eps; //eps is no longer used

  //initializing variables for file-writing
  std::string filename = "";
  std::fstream outfile;

  if(n<=1000){ //if to save storage and to avoid "out of memory"-errors
    //A-matrix for LU-decomposition
    arma::mat A(n,n,arma::fill::zeros);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(i == j) {
          A(i,j) = 2;
        }
        else if(i == j+1) {
            A(i,j) = -1;
        }
        else if(i == j-1) {
            A(i,j) = -1;
        }
    }
    }

    //make a vector of the pointer-array b_tld
    arma::vec b_tld_vec = arma::zeros<arma::vec>(n);
    for(int i = 0; i < n; i++) {
      b_tld_vec[i] = b_tld[i];
    }

    delete[] b_tld; //b_tld is no longer used

    start = ch::steady_clock::now();


    //LU-decomposition
    arma::mat L, U, P;
    arma::lu(L, U, P, A);

    //Solver
    arma::vec z = arma::solve(L,b_tld_vec); //solves Lz=y
    arma::vec v_LU = arma::solve(U,z); //solves Ux=z

    stop = ch::steady_clock::now();

    ch::duration<double> time_span_LU = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used for LU-decomposition and solving = " << time_span_LU.count()  << "s" << std::endl;


    std::cout << "Calculations complete!" << std::endl << "Writing files" << std::endl;

    //data-file for algorithms
    filename = "../../data" + std::to_string(n) + ".txt";
    outfile.open(filename, std::fstream::out | std::ofstream::trunc);
    outfile << "x, v_gen, v_spl, v_LU, u" << std::endl;

    for(int i = 0; i < n; i++) {
      outfile << x[i] << ", ";
      outfile << v_gen[i] << ", ";
      outfile << v_spl[i] << ", ";
      outfile << v_LU[i] << ", ";
      outfile << u[i] << std::endl;
    }
    outfile.close();
  }
  else {
      std::cout << "Calculations complete!" << std::endl << "n is larger than 1000, only special files will be made." << std::endl;
  }

  //deletion of straggling arrays
  delete[] x;
  delete[] u;
  delete[] v_gen;
  delete[] v_spl;


  //data-file for timings
  if(n==1000000){ //if to save storage
    filename = "../../timings" + std::to_string(n) + ".txt";
    outfile.open(filename, std::fstream::out | std::ofstream::app);
    outfile << time_span_general.count() << ", ";
    outfile << time_span_special.count() << std::endl;
    outfile.close();
  }



  //file for error of the standard algorithm (SNIPPET IS COMMENTED OUT TO AVOID DUPLICATE VALUES)
  /*
  filename = "../../error.txt";
  outfile.open(filename, std::fstream::out | std::ofstream::app);
  outfile << std::log10(h) << ", " << max_eps << std::endl;
  outfile.close();
  */


  //end of program
  std::cout << "\nProgram complete!\n";
  return 0;
}
