////g++ problem_2_and_7.cpp -o problem_2_and_7 -std=c++11 -O2 -larmadillo
//// ./problem_2_and_7
//
//#include <iostream>
//#include <armadillo>
//#include <string>
//#include <chrono>
//double poisson(int n);
//// Declaring variables for future exercises
//int N = 1e1;
//int N1 = 1e2;
//int N2 = 1e3;
//int N3 = 1e4;
//int N4 = 1e5;
//int N5 = 1e6;
//int N6 = 1e7;
//
//
//int main(int argc, const char * argv[]){
//    //  Problem 2
//    std::string filename = "exercise_2.txt";
//    int n = 1000;
//    // Defining vectors x and u
//    arma::vec x = arma::linspace(0, 1, n);
//    arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
//    // Writing to txt file
//    std::ofstream ofile;
//    ofile.open(filename);
//    for (int i = 0; i < n; i++){
//        ofile << x(i) << std::setw(12) << std::setprecision(4) << std::scientific << u(i) <<std::endl;}
//    ofile.close();
//
//    poisson(N);
//    poisson(N1);
//    poisson(N2);
//    poisson(N3);
//    poisson(N4);
//    poisson(N5);
//    poisson(N6);
//
//    return 0;
//}
//
//
//double poisson(int N){
//    arma::vec x = arma::linspace(0, 1, N);
//    arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
//    int n = N - 2;
//    // Beginning of exercise 7
//    std::string newfilename = "exercise_7_" + std::to_string(N) + ".txt";
//    arma::vec a = -1*arma::ones(n-1);
//    arma::vec b = 2*arma::ones(n);
//    arma::vec c = -1*arma::ones(n-1);
//    double h = x(1) - x(0);
//    arma::vec g = pow(h, 2)*100*exp(-10*x.subvec(1, n));
//    //std::cout << g.size() << std::endl;
//    //a(0) = 0;
////    arma::vec c_prime = arma::zeros(n-1);
////    arma::vec g_prime = arma::ones(n);
//    arma::vec v = arma::zeros(N);
//    // Using the "initial conditions"
////    c_prime(0) = c(0)/b(0);
////    g_prime(0) = g(0)/b(0);
//
//    // General Thomas algorithm, exercise 7a)
//
//    for(int i=0; i<n-1;i++){
//        double cb_ratio = c(i)/b(i);
//        b(i+1) -= a(i) * cb_ratio;
//        g(i+1) -= g(i) * cb_ratio;
//
//    }
//    v(n) = g(n-1)/b(n-1);
//    for (int i=n-2; i>0; i--){
//        v(i+1) = (g(i) - a(i)*v(i+2))/b(i);
//        //std::cout << v(i+1) << std::endl;
//    }
//
////    for (int i = 1; i < n-1; i++){
////       c_prime(i) = c(i)/(b(i) - a(i)*c_prime(i-1));}
////    for (int i = 1; i < n; i++){
////        g_prime(i) = (g(i)-a(i)*g_prime(i-1))/(b(i) - a(i)*c_prime(i-1));}
////    v(n) = g_prime(n-1);
////    for(int i = n-2; i > 0; i--){
////        v(i+1) = g_prime(i) - c_prime(i)*v(i+2);}
//
//    std::ofstream ofile;
//    ofile.open(newfilename);
//    for (int i = 0; i < N; i++){
//        // Should have added u here to make u and v the same precision.
//        ofile << std::setw(24) << std::setprecision(8) << std::scientific << u(i) << std::setw(24) << x(i) << std::setw(24) << v(i) << std::endl;}
//    ofile.close();
//    return 0;
//
//}



#include <iostream>
#include <armadillo>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
void poisson(int N);
// Declaring variables for future exercises
int N = 1e1 + 2;
int N1 = 1e2 + 2;
int N2 = 1e3 + 2;
int N3 = 1e4 + 2;
int N4 = 1e5 + 2;
int N5 = 1e6 + 2;
int N6 = 1e7 + 2;


int main(){
    //  Problem 2
    std::string filename = "exercise_2.txt";
    int n_2 = 1000;
    // Defining vectors x and u
    arma::vec x_2 = arma::linspace(0, 1, n_2);
    arma::vec u_2 = 1 - (1 - exp(-10))*x_2 - exp(-10*x_2);
    // Writing to txt file
    std::ofstream ofile;
    ofile.open(filename);
    for (int i = 0; i < n_2; i++){
        ofile << x_2(i) << std::setw(12) << std::setprecision(4) << std::scientific << u_2(i) <<std::endl;}
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


void poisson(int N){
    arma::vec x = arma::linspace(0,1,N);
    double h = x(1)-x(0);
    arma::vec u = 1-(1-exp(-10))*x-exp(-10*x);
    int width = 40;
    int prec  = 16;
    int n = N - 2;
    arma::vec a = arma::ones(n)*-1;
    arma::vec b = arma::ones(n)*2;
    arma::vec c = arma::ones(n-1)*-1;
    a(0) = 0;
    arma::vec c_new = arma::zeros(n-1);
    arma::vec g = exp(-10*x.subvec(1,n))*100*pow(h,2);
    arma::vec g_new = arma::zeros(n);
    c_new(0) = c(0)/b(0);
    g_new(0) = g(0)/b(0);
    for (int i = 1; i < n-1; i++){
    c_new(i) = c(i)/(b(i)-a(i)*c_new(i-1));
    }
    for (int i = 1; i < n; i++){
    g_new(i) = (g(i)-a(i)*g_new(i-1))/(b(i)-a(i)*c_new(i-1));
    }

    arma::vec v = arma::zeros(N);
    v(n) = g_new(n-1);
    for (int i = n-2; i >= 0; i--){
        v(i+1) = g_new(i)-c_new(i)*v(i+2);
    }

    // Beginning of exercise 7
    std::string newfilename = "exercise_7_" + std::to_string(n) + ".txt";
    std::ofstream ofile_2;
    ofile_2.open(newfilename);
    for (int i = 0; i < N; i++){
        // Should have added u here to make u and v the same precision.
        //ofile << std::setw(24) << std::setprecision(8) << std::scientific << u(i) << std::setw(24) << x(i) << std::setw(24) << v(i) << std::endl;}
        ofile_2 << std::setw(width) << std::setprecision(prec) << std::scientific << x(i)
              << std::setw(width) << std::setprecision(prec) << std::scientific << v(i)
              << std::endl;}
    ofile_2.close();


}
