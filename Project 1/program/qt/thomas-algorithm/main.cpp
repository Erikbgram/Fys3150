/*
Last changed: 06.09.2019 17:05 by Erlend
*/

#include <iostream> // input and output
#include <cmath> // math
#include <chrono> // timer
#include <string> // string
#include <fstream> // file

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
  double *a = array(n,-1);
  double *d = array(n, 2);
  double *c = array(n,-1);
  double *v = array(n);
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

  v[n-1] = b_tld_new[n-1]/d_new[n-1]; // initial condition for v[-1]
  // backward substitution
  for(int i = n-2; i > 0; i--) {
    v[i] = (b_tld_new[i] - c[i]*v[i+1])/d_new[i];
  }

  ch::steady_clock::time_point stop = ch::steady_clock::now();

  ch::duration<double> time_span_general = ch::duration_cast<ch::nanoseconds>(stop - start);
  std::cout << "Time used for general algorithm = " << time_span_general.count()  << "s" << std::endl;

  double* u = cf_sol(x, n);

  //data-file for general algorithm
  std::string filename = "../../data_general" + std::to_string(n) + ".txt";
  std::fstream outfile;
  outfile.open(filename, std::fstream::out | std::ofstream::trunc);
  outfile << "x    , v    , u    " << std::endl;

  for(int i = 0; i < n; i++) {
    outfile << x[i] << ", ";
    outfile << v[i] << ", ";
    outfile << u[i] << std::endl;
  }
  outfile.close();
  delete[] v;
  delete[] d_new;
  delete[] b_tld_new;
  delete[] ad;


  // initialize arrays for specialized algorithm
  v = array(n);
  d_new = array(n); d_new[0] = d[0];
  b_tld_new = array(n); b_tld_new[0] = b_tld[0];

  // filling d_new
  for(int i = 1; i < n; i++) {
    d_new[i] = (i+1)/i;
  }


  start = ch::steady_clock::now();

  // forward substitution
  for(int i = 1; i < n; i++) {
    b_tld_new[i] = b_tld[i] + (b_tld_new[i-1])/(d_new[i-1]);
  }

  v[n-1] = b_tld_new[n-1]/d_new[n-1]; // initial condition for v[-1]

  //backward substitution
  for(int i = n-2; i > 0; i--) {
    v[i] = (b_tld_new[i] + v[i+1])/d_new[i];
  }

  stop = ch::steady_clock::now();

  ch::duration<double> time_span_special = ch::duration_cast<ch::nanoseconds>(stop - start);
  std::cout << "Time used for specialized algorithm = " << time_span_special.count()  << "s" << std::endl;

  filename = "../../data_special" + std::to_string(n) + ".txt";
  outfile.open(filename, std::fstream::out | std::ofstream::trunc);
  outfile << "x, v, u" << std::endl;

  //data-file for special algorithm
  for(int i = 0; i < n; i++) {
    outfile << x[i] << ", ";
    outfile << v[i] << ", ";
    outfile << u[i] << std::endl;
  }
  outfile.close();
  delete[] a;
  delete[] d;
  delete[] c;
  delete[] v;
  delete[] d_new;
  delete[] b_tld;
  delete[] b_tld_new;

  //data-file for timings
  filename = "../../timings" + std::to_string(n) + ".txt";
  outfile.open(filename, std::fstream::out | std::ofstream::app);
  outfile << time_span_general.count() << ", " << time_span_special.count() << std::endl;

  //end of program
  std::cout << "\nProgram complete!\n";
  return 0;
}
