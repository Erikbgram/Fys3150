/*
  Last edited by Erlend T. North 10:52 17/10/2019
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;

// Functions
double psi(double x1, double y1, double z1, double x2, double y2, double z2, double alpha = 2) { // This function defines the function to integrate
    double value = exp(-2*alpha*(sqrt(x1 * x1 + y1 * y1 + z1 * z1) + sqrt(x2 * x2 + y2 * y2 + z2 * z2)));
    double length = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
      if(length < ZERO) {
          return 0;
        }
    return value / length ;
}

double psi_sphere(double r1, double r2, double t1, double t2, double p1, double p2, double alpha = 2) { // This function defines the function to integrate
    double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
    double value = exp(-3*alpha*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2);
    double length = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosb);
      if(length < ZERO) {
          return 0;
        }
    return (value) / length ;
}

void gauleg(double x1, double x2, double x[], double w[], int n) {
    /*
    This  function takes the lower and upper limits of integration x1, x2
    calculates and returns the abcissas in x[0, ..., n-1] and the weights
    in w[0, ..., n-1] of length n of the Gauss-Legendre n-point quadrature formulae.
    */
    int         m,j,i;
    double      z1,z,xm,xl,pp,p3,p2,p1;
    double      const  pi = 3.14159265359;
    double      *x_low, *x_high, *w_low, *w_high;

    m  = (n + 1)/2; // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    x_low  = x; // pointer initialization
    x_high = x + n - 1;
    w_low  = w;
    w_high = w + n - 1;

    for(i = 1; i <= m; i++) { // loops over desired roots
        z = cos(pi * (i - 0.25)/(n + 0.5));
        /*
        Starting with the above approximation to the ith root
        we enter the mani loop of refinement bt Newtons method.
        */

        do {
            p1 =1.0;
            p2 =0.0;

            // Loop up recurrence relation to get the Legendre polynomial evaluated at x

            for(j = 1; j <= n; j++) {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
        }

            /*
    p1 is now the desired Legrendre polynomial. Next compute
    ppp its derivative by standard relation involving also p2,
    polynomial of one lower order.
    */

            pp = n * (z * p1 - p2)/(z * z - 1.0);
            z1 = z;
            z  = z1 - p1/pp;                   // Newton's method
        }

        while(fabs(z - z1) > ZERO);
        /*
        Scale the root to the desired interval and put in its symmetric
        counterpart. Compute the weight and its symmetric counterpart
        */

        *(x_low++)  = xm - xl * z;
        *(x_high--) = xm + xl * z;
        *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
        *(w_high--) = *(w_low++);
    }
}

double gammln( double xx) {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void gauss_laguerre(double *x, double *w, int n, double alf) {
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            double z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        }
        else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        }
        else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}

double* linspace(double start,double stop, int n) {
  double h = (stop - start)/(n-1);
  double* arr = new double[n];
  for(int i = 0; i < n; i++) {
    arr[i] = start + stop*i*h;
  }
  return arr;
}

int main(int argc, char *argv[]) {
    int n = atoi(argv[1]);
    double la = atof(argv[2]);

    //Legendre
    double *x = new double[n];
    double *w = new double[n];

    // Set up the mesh points and weights
    gauleg(-la, la, x, w, n);

    // Evaluate the integral with the Gauss-Legendre method
    // Note that we initialize the sum. Here brute force gauss-legendre
    double legendre_sum = 0.0;
    /*for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
                for(int l = 0; l < n; l++) {
                    for(int o = 0; o < n; o++) {
                        for(int p = 0; p < n; p++) {
                          legendre_sum += (w[i] * w[j] * w[k] * w[l] * w[o] * w[p]) * psi(x[i], x[j], x[k], x[l], x[o], x[p]);
                        }
                    }
                }
            }
        }
    }*/

    //Laguerre
    double *r = new double[n];
    double *the = new double[n];
    double *phi = new double[n];
    double *wr = new double[n];
    double *wthe = new double[n];
    double *wphi = new double[n];

    gauss_laguerre(r, wr, n, 0);
    gauleg(0, 2*M_PI, the, wthe, n);
    gauleg(0, M_PI, phi, wphi, n);

    double laguerre_sum = 0.0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
                for(int l = 0; l < n; l++) {
                    for(int o = 0; o < n; o++) {
                        for(int p = 0; p < n; p++) {
                          laguerre_sum += (wr[i] * wr[j] * wthe[k] * wthe[l] * wphi[o] * wphi[p]) * psi_sphere(r[i], r[j], the[k], the[l], phi[o], phi[p]);
                        }
                    }
                }
            }
        }
    }



    double exact = (5*M_PI*M_PI)/(16*16);

    // Final output
    cout  << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Gaussian-Legendre quad = " << setw(20) << setprecision(15)  << legendre_sum << endl;
    cout << "Exact answer = " << setw(20) << setprecision(15) << exact << endl;
    cout << "E-ror = " << setw(20) << setprecision(15) << fabs(exact-legendre_sum) << endl;
    cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15)  << laguerre_sum << endl;
    cout << "Exact answer = " << setw(20) << setprecision(15) << exact << endl;
    cout << "E-ror = " << setw(20) << setprecision(15) << fabs(exact-laguerre_sum) << endl;
    delete [] x;
    delete [] w;

    fstream outfile;
    outfile.open("../../error.txt", std::fstream::out | std::ofstream::app);
    outfile << n << " , " << la << " , " << fabs(exact-legendre_sum) << endl;
    outfile.close();

  return 0;
}