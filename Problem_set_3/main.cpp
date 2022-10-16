#include<particle.hpp>
#include<armadillo>
#include<iostream>
#include<penningtrap.hpp>



int main()
{   // Magnetic field strength (T) and electric potential (V)
    double k_e = 1.38935333*pow(10,5);
    double T = 9.64852558*pow(10,1);
    double V = 9.64852558*pow(10,7);

// Penning trap configuration constants
    double B0 = 1.*T;
    double V0 = 25.*pow(10,-3)*V;//*pow(10, -3)*V; //10*V (Original definition. Changed after message)
    double d = 500.; //micrometer
 
   double atomic_mass = 1.66053906660*pow(10,-27);
   double mass_Ca_I = 40.078;
   double c = 299792458.;
   double elementary_charge = 1.60217663*1e-19;
   int q = 1.;
   double mass_Ca_II = (11.87172*elementary_charge/(pow(c,2)*atomic_mass) + mass_Ca_I)-5.4858*pow(10,-4);
   double omega_0_squared =  pow(q,2)*pow(B0,2)/pow(mass_Ca_II,2);
   double omega_z_squared = 2*q*V0/(mass_Ca_II*pow(d,2));
   
   

   
   std::cout << mass_Ca_II << std::endl;
   arma::vec position1 = {20.,0.,20.};
   arma::vec velocity1 = {0.,25.,0.};
   arma::vec position2 = {25.,25.,0.};
   arma::vec velocity2 ={0.,40.,5.};
   Particle test_particle_1 = Particle(q,mass_Ca_II,position1,velocity1);
   Particle test_particle_2 = Particle(q,mass_Ca_II,position2,velocity2);
   PenningTrap test_trap = PenningTrap(B0,V0,d);
   test_trap.add_particle(test_particle_1);
   test_trap.add_particle(test_particle_2);


   int time = 50;
   test_trap.RK4(time,4000,"yes"); //"yes" for interactions


   return 0;

}