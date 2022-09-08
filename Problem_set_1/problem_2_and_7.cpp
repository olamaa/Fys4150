//g++ problem_2_and_7.cpp -o problem_2_and_7 -std=c++11 -O2 -larmadillo
// ./problem_2_and_7

#include <iostream>
#include <armadillo>
#include <string>
#include <chrono>
double poisson(int n);
// Declaring variables for future exercises
int N = 1e1;
int N1 = 1e2;
int N2 = 1e3;
int N3 = 1e4;
int N4 = 1e5;
int N5 = 1e6;
int N6 = 1e7;


int main(int argc, const char * argv[]){
    //  Problem 2
    std::string filename = "exercise_2.txt";
    int n = 1000;
    // Defining vectors x and u
    arma::vec x = arma::linspace(0, 1, n);
    arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
    // Writing to txt file
    std::ofstream ofile;
    ofile.open(filename);
    for (int i = 0; i < n; i++){
        ofile << x(i) << std::setw(12) << std::setprecision(4) << std::scientific << u(i) <<std::endl;}
    ofile.close();
    
    poisson(N);
    poisson(N1);
    poisson(N2);
    poisson(N3);
    poisson(N4);
    poisson(N5);
    poisson(N6);
    
    return 0;
}


double poisson(int N){
    arma::vec x = arma::linspace(0, 1, N);
    arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
    int n = N - 2;
    // Beginning of exercise 7
    std::string newfilename = "exercise_7_" + std::to_string(N) + ".txt";
    arma::vec a = -1*arma::ones(n);
    arma::vec b = 2*arma::ones(n);
    arma::vec c = -1*arma::ones(n-1);
    double h = x(1) - x(0);
    arma::vec g = pow(h, 2)*100*exp(-10*x.subvec(1, n));
    a(0) = 0;
    arma::vec c_prime = arma::zeros(n-1);
    arma::vec g_prime = arma::ones(n);
    arma::vec v = arma::zeros(N);
    // Using the "initial conditions"
    c_prime(0) = c(0)/b(0);
    g_prime(0) = g(0)/b(0);
    
    // General Thomas algorithm, exercise 7a)
    for (int i = 1; i < n-1; i++){
       c_prime(i) = c(i)/(b(i) - a(i)*c_prime(i-1));}
    for (int i = 1; i < n; i++){
        g_prime(i) = (g(i)-a(i)*g_prime(i-1))/(b(i) - a(i)*c_prime(i-1));}
    v(n) = g_prime(n-1);
    for(int i = n-2; i > 0; i--){
        v(i) = g_prime(i) - c_prime(i)*v(i+1);}
    
    std::ofstream ofile;
    ofile.open(newfilename);
    for (int i = 0; i < N; i++){
        // Should have added u here to make u and v the same precision.
        ofile << x(i) << std::setw(12) << std::setprecision(4) << std::scientific << v(i) << std::setw(12) << n << std::endl;}
    ofile.close();
    return 0;
    
}
