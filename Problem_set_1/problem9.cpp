#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>


void solver();

int N = 1000; // Defines how many x-values we want to calculate u for, including boundaries
int n = N-2; // Number of steps, since we do not want to calculate the boundary conditions
arma::vec x = arma::linspace(0,1,N);
double h = x(1)-x(0);
arma::vec u = 1-(1-exp(-10))*x-exp(-10*x);
int width = 12;
int prec  = 4;

int main()
{
    solver();
    return 0;
}
//7n FLOPs
void solver()
{   
    // Defining vectors
    arma::vec c_new = arma::zeros(n-1);
    arma::vec g = exp(-10*x.subvec(1,n))*100*pow(h,2);
    arma::vec g_new = arma::zeros(n);
    c_new(0) = -1/2;
    g_new(0) = g(0)/2;
    //Forwards substitution
    for (int i = 1; i < n-1; i++){
    c_new(i) = -1/(2+c_new(i-1));    
    }

    for (int i = 1; i < n; i++){
    g_new(i) = (g(i)+g_new(i-1))/(2+c_new(i-1));
    }

    arma::vec v = arma::zeros(N); // Length of v is N since we want to include boundary conditions. n is N-2
    v(n) = g_new(n-1);
    //Backwards substitution
    for (int i = n-2; i > 0; i--){
        v(i) = g_new(i)-c_new(i)*v(i+1);
    }
    
    //Writing to file
    std::string filename2 = "x_v9.txt";
    std::ofstream ofile2;
    ofile2.open(filename2);
    
    for (int i = 0; i < N; i++){
    ofile2 << std::setw(width) << std::setprecision(prec) << std::scientific << x(i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << v(i)
          << std::endl;
    }
    
    ofile2.close();
}    

