#include<particle.hpp>
#include<armadillo>
#include<iostream>
#include<vector>
#include<penningtrap_beehive.hpp>
#include<cmath>


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
   arma::vec initial_position1 = {20.,0.,20.};
   arma::vec initial_velocity1 = {0.,25.,0.};
   arma::vec initial_position2 = {25.,25.,0.};
   arma::vec initial_velocity2 ={0.,40.,5.};
   arma::vec position1 = {20.,0.,20.};
   arma::vec velocity1 = {0.,25.,0.};
   arma::vec position2 = {25.,25.,0.};
   arma::vec velocity2 ={0.,40.,5.};
   //Particle test_particle_1 = Particle(q,mass_Ca_II,position1,velocity1);
   //Particle test_particle_2 = Particle(q,mass_Ca_II,position2,velocity2);
   PenningTrap test_trap = PenningTrap(B0,V0,d);
   //test_trap.add_particle(test_particle_1);
   //test_trap.add_particle(test_particle_2);

    
   //int time = 50;

//   for (int i=2;i<6;i++){
//    test_trap.RK4(time,1000*pow(2,i),"yes","yes");
//    test_trap.remove_all_particles();
//
//    test_trap.add_particle(test_particle_1);
//    test_trap.add_particle(test_particle_2);
//    
//    test_trap.RK4(time,1000*pow(2,i),"no","yes");
//    test_trap.remove_all_particles();
//
//    test_trap.add_particle(test_particle_1);
//    test_trap.add_particle(test_particle_2);
//    
//    test_trap.forward_Euler(time,1000*pow(2,i),"yes");
//    test_trap.remove_all_particles();
//   
//    test_trap.add_particle(test_particle_1);
//    test_trap.add_particle(test_particle_2);
//
//    test_trap.forward_Euler(time,1000*pow(2,i),"no");
//    
//   
//   }

   //test_trap.task9(d,q,mass_Ca_II,{0.1,0.4,0.7},arma::linspace(0.8,2,(2-0.8)/0.002),500,40000,"without");
   //test_trap.task9(d,q,mass_Ca_II,{0.1,0.4,0.7},500,40000,"with",arma::linspace(1,2,1/0.02),arma::linspace(0.2,2,(1.6-1.2)/0.04),arma::linspace(1.1,1.7,(1.7-1.1)/0.06));
   test_trap.task9(d,q,mass_Ca_II,0.1,500,40000,"with",arma::linspace(1,2,1./0.005));



    //test_trap.task9(d,q,mass_Ca_II,0.4,500,40000,"with",arma::linspace(0.8,2.,(2.-0.8)/0.006));
    //test_trap.task9(d,q,mass_Ca_II,0.7,500,40000,"with",arma::linspace(0.8,2.5,(2.5-0.8)/0.0085));
   return 0;

}