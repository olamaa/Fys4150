//g++ problem_10.cpp -o problem_10 -std=c++11 -O2 -larmadillo
// ./problem_10

#include <iostream>
#include <armadillo>
#include <string>
#include <chrono>
double general(int n);
double specialized(int n);
// Declaring variables
int N = 1e1;
int N1 = 1e2;
int N2 = 1e3;
int N3 = 1e4;
int N4 = 1e5;
int N5 = 1e6;
int N6 = 1e7;


int main(int argc, const char * argv[]){
    // Problem 10
    // Declaring vectors
    int h = 100;
    arma::vec time_avg_N_general = arma::zeros(h);
    arma::vec time_avg_N_general1 = arma::zeros(h);
    arma::vec time_avg_N_general2 = arma::zeros(h);
    arma::vec time_avg_N_general3 = arma::zeros(h);
    arma::vec time_avg_N_general4 = arma::zeros(h);
    arma::vec time_avg_N_general5 = arma::zeros(h);
    arma::vec time_avg_N_general6 = arma::zeros(h);
    arma::vec time_avg_N_specialized = arma::zeros(h);
    arma::vec time_avg_N_specialized1 = arma::zeros(h);
    arma::vec time_avg_N_specialized2 = arma::zeros(h);
    arma::vec time_avg_N_specialized3 = arma::zeros(h);
    arma::vec time_avg_N_specialized4 = arma::zeros(h);
    arma::vec time_avg_N_specialized5 = arma::zeros(h);
    arma::vec time_avg_N_specialized6 = arma::zeros(h);
    
    for (int i = 0; i < h; i++){
        // Calling on the general function
        time_avg_N_general(i) = general(N);
        time_avg_N_general1(i) = general(N1);
        time_avg_N_general2(i) = general(N2);
        time_avg_N_general3(i) = general(N3);
        time_avg_N_general4(i) = general(N4);
        time_avg_N_general5(i) = general(N5);
        time_avg_N_general6(i) = general(N6);
        time_avg_N_specialized(i) = specialized(N);
        time_avg_N_specialized1(i) = specialized(N1);
        time_avg_N_specialized2(i) = specialized(N2);
        time_avg_N_specialized3(i) = specialized(N3);
        time_avg_N_specialized4(i) = specialized(N4);
        time_avg_N_specialized5(i) = specialized(N5);
        time_avg_N_specialized6(i) = specialized(N6);}
    // Printing the average time it takes to complete the function for different number of steps
std::cout << "First we present with " + std::to_string(h) + " rounds, the average amount of time for the generalized Thomas algortihm" << std::endl;
    std::cout << "For N = " + std::to_string(N)+ "        we have time: " + std::to_string(arma::sum(time_avg_N_general)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N1)+ "       we have time: " + std::to_string(arma::sum(time_avg_N_general1)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N2)+ "      we have time: " + std::to_string(arma::sum(time_avg_N_general2)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N3)+ "     we have time: " + std::to_string(arma::sum(time_avg_N_general3)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N4)+ "    we have time: " + std::to_string(arma::sum(time_avg_N_general4)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N5)+ "   we have time: " + std::to_string(arma::sum(time_avg_N_general5)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N6)+ "  we have time: " + std::to_string(arma::sum(time_avg_N_general6)/h)+" seconds" <<std::endl;

    std::cout << "Secondly we present with " + std::to_string(h) + " rounds, the average amount of time for the specialized Thomas algortihm" << std::endl;
    std::cout << "For N = " + std::to_string(N)+ "        we have time: " + std::to_string(arma::sum(time_avg_N_specialized)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N1)+ "       we have time: " + std::to_string(arma::sum(time_avg_N_specialized1)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N2)+ "      we have time: " + std::to_string(arma::sum(time_avg_N_specialized2)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N3)+ "     we have time: " + std::to_string(arma::sum(time_avg_N_specialized3)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N4)+ "    we have time: " + std::to_string(arma::sum(time_avg_N_specialized4)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N5)+ "   we have time: " + std::to_string(arma::sum(time_avg_N_specialized5)/h)+" seconds" <<std::endl;
    std::cout << "For N = " + std::to_string(N6)+ "  we have time: " + std::to_string(arma::sum(time_avg_N_specialized6)/h)+" seconds" <<std::endl;
    
    return 0;
}


double specialized(int N){
    // Hitting the stopwatch, declaring vectors and variables
    auto t0 = std::chrono::high_resolution_clock::now();
    arma::vec x = arma::linspace(0, 1, N);
    arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
    int n = N - 2;

    double h = x(1) - x(0);
    arma::vec g = pow(h, 2)*100*exp(-10*x.subvec(1, n));
    arma::vec c_prime = arma::zeros(n-1);
    arma::vec g_prime = arma::ones(n);
    arma::vec v = arma::zeros(N);
    // Using the "initial conditions"
    c_prime(0) = -1/2;
    g_prime(0) = g(0)/2;
    
    // General Thomas algorithm, exercise 7a)
    for (int i = 1; i < n-1; i++){
       c_prime(i) = -1/(2 + c_prime(i-1));}
    for (int i = 1; i < n; i++){
        g_prime(i) = (g(i) + g_prime(i-1))/(2 + c_prime(i-1));}
    v(n) = g_prime(n-1);
    for(int i = n-2; i > 0; i--){
        v(i) = g_prime(i) - c_prime(i)*v(i+1);}
    // Hitting the stopwatch
    auto t1 = std::chrono::high_resolution_clock::now();
    double duration_seconds_N = std::chrono::duration<double>(t1 - t0).count();
    return duration_seconds_N;
    
}


double general(int N){
    // Hitting the stopwatch, declaring vectors and variables
    auto t0 = std::chrono::high_resolution_clock::now();
    arma::vec x = arma::linspace(0, 1, N);
    arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
    int n = N - 2;
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
    
    // Hitting the stopwatch
    auto t1 = std::chrono::high_resolution_clock::now();
    double duration_seconds_N = std::chrono::duration<double>(t1 - t0).count();
    return duration_seconds_N;
    
}
