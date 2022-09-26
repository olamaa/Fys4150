//g++ problem2.cpp -o problem2 -std=c++11 -O2 -larmadillo
// ./problem2


#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

arma::mat A_matrix(int &n)
{   int N = n-1;
    arma::mat A = arma::zeros(N,N);
    arma::vec x = arma::linspace(0,1,N+2);
    double h = 1./n;
    double a = -1./pow(h,2);
    double d = 2./pow(h,2);
    for (int row = 0; row < N-1; row++){
        for (int column = 0; column < N-1; column++){
            if (row == column){
                A(row,column) = d;
                A(row,column+1) = a;
                A(row+1,column) =a;
            }
        }
    }
    A(N-1,N-1) = d;
    return A;
}

void test(int &n){
    double h = 1./n;
    double a = -1./pow(h,2);
    double d = 2./pow(h,2);
    arma::vec lambda = arma::zeros(n-1);
    arma::mat eigenvecs = arma::zeros(n-1,n-1);
    for (int i =1; i<n;i++){
        lambda(i-1) = d + 2.*a*cos(i*M_PI/n);
        for (int element =1;element<n;element++){
            eigenvecs(i-1,element-1) = sin(element*i*M_PI/n);
        }
    }
    std::cout << lambda << std::endl;
    std::cout << arma::normalise(arma::trans(eigenvecs)) << std::endl;
}

double find_max(arma::mat& A, int &k, int &l){ //largest off-diagonal values (3a)
//    Symmetric implies that it is sufficient to look at the upper triangle
    arma::mat upper_triang = arma::trimatu(A, 1);
    arma::mat non_diag = arma::abs(upper_triang);
    arma::uword x = non_diag.index_max();
    double max_val = non_diag(x);
    k = x % A.n_rows;
    l = x / A.n_cols;
    return max_val;
}
    
    
void all(int &k, int &l, double &t, double &c, double &s, double &tau,arma::mat &R,arma::mat &S, arma::mat &A ){//4a)
    int N = A.n_cols;
    int iterations = 0;
    while (find_max(A,k,l) > pow(10,-8)){
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0){
            t = 1/(tau + sqrt(1+pow(tau, 2)));
        }
        else {
            t = -1/(-tau + sqrt(1+pow(tau,2)));
        }
        c = 1/sqrt(1+pow(t, 2));
        s = t*c;
        S(k, k) = c;
        S(l, l) = c;
        S(k, l) = -s;
        S(l, k) = s;
        
        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
        a_kk = A(k,k);
        a_ll = A(l,l);
        // changing the matrix elements with indices k and l
        A(k,k) = a_kk*pow(c,2) - 2.0*A(k,l)*c*s + a_ll*pow(s,2);
        A(l,l) = a_kk*pow(s,2) + 2.0*A(k,l)*c*s + a_ll*pow(c,2);
        A(k,l) = 0.0; // hard-coding of the zeros
        A(l,k) = 0.0;
        // and then we change the remaining elements
        for ( int i = 0; i < N; i++ ) {
            if ( i != k && i != l ){
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k) = a_ik*c - a_il*s;
                A(k,i) = A(i,k);
                A(i,l) = a_il*c + a_ik*s;
                A(l,i) = A(i,l);
            }
            // Finally, we compute the new eigenvectors
            r_ik = R(i,k);
            r_il = R(i,l);
            R(i,k) = r_ik*c - r_il*s;
            R(i,l) = r_il*c + r_ik*s;
            iterations += 1;
        }
    }
    arma::vec A_diag = A.diag();
    arma::vec &A_diag_change = A_diag;
    arma::vec saved_indicies = arma::zeros(3); //the three lowest eigenvalues
    int lambda_val = 0;
    int &lambda = lambda_val;
    for (int i = 0; i < 3; i++){
        lambda = A_diag_change.index_min();
        saved_indicies(i) = lambda;
        A_diag_change(lambda) = arma::datum::inf;
    }
    int width = 40;
    int prec  = 16;
    std::string newfilename = "exercise_6.txt";
    std::ofstream ofile_1;
    arma::vec x = arma::linspace(0, 1, N+2);
    arma::vec x_sub = x.subvec(1, N);
    arma::mat R_trans = arma::trans(R);
    ofile_1.open(newfilename);
    ofile_1 << std::setw(width) << 0
            << std::setw(width) << 0
            << std::setw(width) << 0
            << std::setw(width) << 0
            << std::endl;
    for (int i = 0; i < N; i++){
        ofile_1 << std::setw(width) << x_sub(i)
                << std::setw(width) << R_trans(saved_indicies(0), i)
                << std::setw(width) << R_trans(saved_indicies(1), i)
                << std::setw(width) << R_trans(saved_indicies(2), i)
                << std::endl;}
    
    ofile_1 << std::setw(width) << 1
            << std::setw(width) << 0
            << std::setw(width) << 0
            << std::setw(width) << 0
            << std::endl;
    ofile_1.close();
}
                      

void test_matrix(){ //Test code (3b)
    arma::mat test_m = arma::eye(4, 4);
    test_m(3, 0) = 0.5;
    test_m(2, 1) = -0.7;
    test_m(1, 2) = -0.7;
    test_m(0, 3) = 0.5;
    int k = 0;
    int l = 0;
    find_max(test_m, k, l);
}

int main()
{   int n = 100;
    int N = n-1; // (4b, N = 6)
    arma::mat A = A_matrix(n);
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    //std::cout << eigval << std::endl;
    //std::cout << arma::normalise(-1*eigvec) << std::endl;
    //test_matrix();
    
    double t = 0;
    double c = 0;
    double s = 0;
    double tau = 0;
    arma::mat R = arma::eye(N,N);
    arma::mat S = arma::eye(N,N);
    int k = 0;
    int l = 0;
    double test_val = find_max(A, k, l);
    all(k,l,t,c,s,tau,R,S,A);
    //test(n);
    
    
    
    // Generate random N*N matrix
    //arma::mat B = arma::mat(N, N).randn();
    // Symmetrize the matrix by reflecting the upper triangle to lower triangle
    //B = arma::symmatu(B);
    //all(k,l,t,c,s,tau,R,S,B);
    return 0;
}
