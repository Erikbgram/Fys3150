#include <iostream>
#include <armadillo>
#include <chrono>

using namespace std;
namespace ch = std::chrono;

//  the offdiag function, using Armadillo
void offdiag(arma::mat A, int *p, int *q, int n); {
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              max = aij;  p = i; q = j;
           }
       }
   }
}


int main()
{
    cout << "Hello World!" << endl;
    return 0;
}
