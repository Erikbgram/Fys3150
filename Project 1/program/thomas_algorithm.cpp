#include <iostream> // input and output
#include <cmath> // math
#include <ctime> // time
#include <iomanip> // percision
#include <string> // string
#include <sstream> // string contactinating
#include <fstream> // file

//create dynamic array and fill with num
double* array(int n, double num=0) {
  double* arr = new double[n];
  for(int i = 0; i < n; i++) {
    arr[i] = num;
  }
  return arr;
}

//delete dynamic array
void del_array(double* arr) {
  delete arr;
}

//f(x)=100e^(-10x)
double* func(double* x, int n) {
  double* f = new double[n];
  double h = (x[n-1] - x[0]) / (n-1);
  std::cout << "h = " << h << std::endl;
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

  //defining array
  double* a = array(n,-1);
  double* d = array(n, 2);
  double* c = array(n,-1);
  double* v = array(n);
  double* d_new = array(n); d_new[0] = d[0];
  double* b_tld = func(x, n);
  double* b_tld_new = array(n);
  double* ad = array(n); //We can save 1 FLOP by pre-calculating a[i-1]/d_new[i-1]

  //clock
  clock_t start, stop;
  start = clock();

  //forward substitution
  for(int i = 1; i < n; i++) {
    ad[i-1] = a[i-1]/d_new[i-1]; // measure once, cut twice
    d_new[i] = d[i] - ad[i-1]*c[i-1];
    b_tld_new[i] = b_tld[i] - b_tld_new[i-1]*ad[i-1];
  }

  v[n-1] = b_tld_new[n-1]/d_new[n-1]; //initial condition for v[-1]

  for(int i = n-2; i > 0; i--) {
    v[i] = (b_tld_new[i] - c[i]*v[i+1])/d_new[i];
  }

  stop = clock();

  double timeused = (double) (stop - start)/(CLOCKS_PER_SEC );
  std::cout << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  std::cout << std::setprecision(10) << std::setw(20) << "Time used for algorithm = " << timeused  << "s" << std::endl;

  double* u = cf_sol(x, n);

  std::stringstream filename;
  filename << "data" << n << ".txt";
  std::fstream outfile;
  outfile.open(filename.str(), std::fstream::out | std::ofstream::trunc);
  outfile << "x    , v    , u    " << std::endl;

  for(int i = 0; i < n; i++) {
    outfile << std::setprecision(5) << std::setw(5) << x[i] << ", " << std::setprecision(5) << std::setw(7) << v[i] << ", " << std::setprecision(5) << std::setw(7) << u[i] << std::endl;
  }

  outfile.close();

  //end of program
  std::cout << "\nProgram complete!\n";
  return 0;
}
