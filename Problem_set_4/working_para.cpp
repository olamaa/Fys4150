#include<armadillo>
#include<iostream>
#include<vector>
#include<cstdlib>
#include<iomanip>
#include "omp.h"
#include<chrono>


inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}
arma::mat initialize(int N, double temp,arma::vec& w,double& beta,double&E_tot,double&M_tot,std::string& rand){
    beta = 1/temp;
    E_tot = 0;
    M_tot = 0;
    arma::mat spin_matrix = arma::ones(N,N);
    if (rand=="True"){
        rand = "random";
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<int> dis(1,2);
        for (int j = 0; j < N; j++){
            for (int i = 0; i < N; i++){
                int random_number = dis(gen);
                if(random_number == 2){
                    spin_matrix(j,i) *= -1;

                }

            }
        }
    }
    else{  
        rand = "not_random";
    }

    for (int j = 0;j<N;j++){
        for(int i = 0;i<N;i++){
            E_tot -= spin_matrix(j,i)*spin_matrix(j,periodic(i,N,1));
            E_tot -= spin_matrix(i,j)*spin_matrix(periodic(i,N,1),j);
            M_tot += spin_matrix(j,i);
        }
    }
    for (int i = 1;i<3;i++){
        w(i-1) = exp(-beta*(4*i));
    }
    return spin_matrix;
}

void metropolis(int N,arma::mat& spin_matrix,arma::vec w,double beta,double&E_tot,double&M_tot){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(0.,1.);
    std::uniform_int_distribution<int> dist_int(0,N-1);
    int x;
    int y;
    int delta_E;
    for (int j = 0; j < N; j++){
        for (int i = 0; i < N; i++){
        
            x = int(dist_int(gen));
            y = int(dist_int(gen));

            delta_E = 2*spin_matrix(y,x)*(spin_matrix(periodic(y,N,1),x)+spin_matrix(periodic(y,N,-1),x) + spin_matrix(y,periodic(x,N,1)) + spin_matrix(y,periodic(x,N,-1)));
            if (delta_E <= 0){
                spin_matrix(y,x) *=-1;
                E_tot += delta_E;
                M_tot += 2*spin_matrix(y,x);
            }    
            else{
                if (delta_E == 8){
                    if (dist(gen)<= w(1)){
                        spin_matrix(y,x) *=-1;
                        E_tot += delta_E;
                        M_tot += 2*spin_matrix(y,x);
                    }
                    else{
                    }
                }
                else{
                    if (dist(gen)<= w(0)){
                        spin_matrix(y,x) *=-1;
                        E_tot += delta_E;
                        M_tot += 2*spin_matrix(y,x);
                    }
                    else{
                    }
                }   
            }         
        } 
    }
}

void operations(int N, double temp,int number_of_cycles,std::string rand, std::string write_to_file_1,std::string write_to_file_2){   
    arma::vec w = {0,0};
    double n_spins = N*N;
    double beta;
    double E_tot;
    double M_tot;
    double tot = number_of_cycles;
    arma::mat averages = arma::zeros(4,number_of_cycles);
    arma::mat spin_matrix = initialize(N,temp,w,beta,E_tot,M_tot,rand);
    for (int cycle=0;cycle<number_of_cycles;cycle++){
        metropolis(N,spin_matrix,w,beta,E_tot,M_tot);
        averages(0,cycle) = E_tot;
        averages(1,cycle) = E_tot*E_tot;
        averages(2,cycle) = abs(M_tot);
        averages(3,cycle) = M_tot*M_tot;
        }
    if (write_to_file_1=="True"){
    int width = 20;
    int prec = 12; 
    std::ofstream ofile;
    std::string filename = "L_" + std::to_string(N) + "_T_"  + std::to_string(temp) + "_" + rand + ".txt";
    ofile.open(filename);
    ofile<< std::setw(width) << "i" <<
                std::setw(width) << "<eps>"<< 
                std::setw(width) << "<|m|>"<<
                std::setw(width) << "eps" <<
                std::endl; 
    double avg1=0;
    double avg2=0;
    double avg_1;
    double avg_2;

    for (int i = 1;i<=number_of_cycles;i++){
        avg1 += averages(0,i-1)/(n_spins);
        avg2 += averages(2,i-1)/(n_spins);
        ofile<< std::setw(width)<<std::setprecision(prec) << i-1 <<
                std::setw(width)<<std::setprecision(prec) << avg1/i<< 
                std::setw(width)<<std::setprecision(prec) << avg2/i<<
                std::setw(width)<<std::setprecision(prec) << averages(0,i-1)/n_spins <<
                std::endl;       
    }
    ofile.close();
    }
    else{
    }
    double average_E=0;
    double average_eps=0;
    double average_E_squared=0;
    double average_eps_squared = 0;
    double average_abs_M = 0;
    double average_abs_m = 0;
    double average_M_squared = 0;
    double average_m_squared = 0;
    for (int i =15000;i<number_of_cycles;i++){
        average_E += averages(0,i);
        average_E_squared += averages(1,i);
        average_abs_M += averages(2,i);
        average_M_squared += averages(3,i);
    }
    double effective_number_of_cycles = number_of_cycles-15000;
    average_E = average_E/effective_number_of_cycles;
    average_E_squared = average_E_squared/effective_number_of_cycles;
    average_M_squared = average_M_squared/effective_number_of_cycles;
    average_abs_M = average_abs_M/effective_number_of_cycles;
    average_eps = average_E/n_spins;
    average_eps_squared = average_E_squared/pow(n_spins,2);
    average_abs_m = average_abs_M/n_spins;
    average_m_squared = average_M_squared/pow(n_spins,2);
    double C_v_num = 1./n_spins*(beta/temp)*(average_E_squared-pow(average_E,2));
    double Chi_num = 1./n_spins*beta*(average_M_squared-pow(average_abs_M,2));
    if (write_to_file_2=="True"){
    int width = 20;
    int prec = 12; 
    std::ofstream ofile;
    std::string filename = "8a_beehive.txt";
    ofile.open(filename, std::ofstream::app);
        ofile<< std::setw(width)<<std::setprecision(prec) << temp <<
                std::setw(width)<<std::setprecision(prec) << N <<
                std::setw(width)<<std::setprecision(prec) << average_eps<< 
                std::setw(width)<<std::setprecision(prec) << average_abs_m <<
                std::setw(width)<<std::setprecision(prec) << C_v_num <<
                std::setw(width)<<std::setprecision(prec) << Chi_num <<
                std::endl;   
        ofile.close();
        }
    else{
    }
    if (N ==2){
    double average_E=0;
    double average_eps=0;
    double average_E_squared=0;
    double average_eps_squared = 0;
    double average_abs_M = 0;
    double average_abs_m = 0;
    double average_M_squared = 0;
    double average_m_squared = 0;
    for (int i =0;i<number_of_cycles;i++){
        average_E += averages(0,i);
        average_E_squared += averages(1,i);
        average_abs_M += averages(2,i);
        average_M_squared += averages(3,i);
    }

    average_E = average_E/number_of_cycles;
    average_E_squared = average_E_squared/number_of_cycles;
    average_M_squared = average_M_squared/number_of_cycles;
    average_abs_M = average_abs_M/number_of_cycles;
    average_eps = average_E/n_spins;
    average_eps_squared = average_E_squared/pow(n_spins,2);
    average_abs_m = average_abs_M/n_spins;
    average_m_squared = average_M_squared/pow(n_spins,2);

    double C_v_num = 1./n_spins*(beta/temp)*(average_E_squared-pow(average_E,2));
    double chi_num = 1./n_spins*beta*(average_M_squared-pow(average_abs_M,2));
    
    double partition_function = 2*(exp(8*beta)+exp(-8*beta)+6);
    double eps_avg = -4/partition_function*(exp(beta*8)-exp(-beta*8));
    double abs_m_avg = 2/partition_function*(exp(8*beta)+1);
    double C_v_avg = 1/4.*beta/temp*(4*4*8/partition_function*(exp(beta*8)+exp(-beta*8))-4*4/(partition_function*partition_function)*pow((4*exp(beta*8)-4*exp(-beta*8)),2));
    double chi_avg = beta/4.*(4.*4.*2./partition_function*(exp(beta*8.)+1.)-4.*4.*2.*2./pow(partition_function,2)*pow(exp(beta*8.)+1.,2));
    std::cout << "Property                 " << "Numerical    " << "Analytical"<< std::endl;
    std::cout << "Average energy per spin  " << average_eps << "    "<< eps_avg << std::endl;
    std::cout <<"Average |m|               " << average_abs_m << "    "<< abs_m_avg << std::endl;
    std::cout << "Specific heat capacity   " << C_v_num << "    " << C_v_avg << std::endl;
    std::cout << "The susceptibility       " << chi_num << "    " << chi_avg << std::endl;
    }

}

int main(){
    //task 4
    operations(2,1,1e5,"True","False","False");
    
    //task 5 / 6
    operations(20,1,1e5,"True","True","False");
    operations(20,2.4,1e5,"True","True","False");
    operations(20,1,1e5,"False","True","False");
    operations(20,2.4,1e5,"False","True","False");
    
    
    //task 8
    int num_t = 100;
    arma::vec T = arma::linspace(2.1,2.4,num_t);
    arma::vec L = {40, 60, 80, 100};
    int width = 20;
    
    std::ofstream ofile;
    std::string filename = "8a_beehive.txt";
    ofile.open(filename);
    ofile<< std::setw(width) << "T" <<
            std::setw(width) << "L"<<
            std::setw(width) << "<eps>"<<
            std::setw(width) << "<|m|>"<<
            std::setw(width) << "C_v"<<
            std::setw(width) << "chi"<<
            std::endl;
    ofile.close();
    double count = 0.;
    double tot_num = arma::size(T)(0)*arma::size(L)(0);
    #pragma omp parallel for collapse(2)
    //for loop to parallelize
    for (int i = 0;i < arma::size(L)(0);i++){
        for (int j = 0;j< num_t;j++){
            operations(L(i),T(j),1e6,"True","False","True");
            count +=1.;
            std::cout << (count/tot_num)*100. << "%" << std::endl;

        }
    }
    
    
    //task 7
    int N_cores = 40;
    L = arma::ones(N_cores)*3;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0;i<N_cores;i++){
        operations(L(i),2.4,1e5,"True","False","False");
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0;i<N_cores;i++){
        operations(L(i),2.4,1e5,"True","False","False");
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t1-t0).count();
    double time1 = std::chrono::duration<double>(t3-t2).count();
    std::cout <<"Timing test: without parallell;" << time << ", " << "with parallell;"<< time1 << ", with a speed up factor of "<< time/time1 << std::endl;
} 
