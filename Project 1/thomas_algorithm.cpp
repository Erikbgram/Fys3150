#include <iostream>
#include <cmath>
#include <ctime>


//create dynamic array
double* dyn_array(int n) {
  double* arr = new double[n];
  for(int i = 0; i < n; i++) {
    arr[i] = 0;
  }
  return arr;
}
//delete dynamic array
void del_dyn_array(double* arr) {
  delete arr;
}

//f(x)=100e^(-10x)
double* func(double* x, int n) {
  double* f = new double[n];
  for(int i = 0; i < n; i++) {
    f[i] = 100*exp(10*x[i]);
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

//linspace for dyn_array
double* linspace(double start,double stop, int n) {
  double h = (stop - start)/(n-1);
  double* x = new double[n];
  for(int i = 0; i < n; i++) {
    x[i] = start + stop*i*h;
  }
  return x;
}


int main() {
  //define
  int n;
  std::cout << "Choose your n: ";
  std::cin >> n;

  //x array
  double* x = linspace(0,1,n);

  //step length
  double h = (x[n-1] - x[0]) / (n-1);
  std::cout << "h = " << h << std::endl;




  //end of program
  std::cout << "\nProgram complete!\n";
  return 0;
}
