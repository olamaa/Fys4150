#ifndef __penningtrap_hpp__
#define __penningtrap_hpp__
#include<iostream>
#include<particle.hpp>
#include<armadillo>
#include<fstream>
#include<vector>

class PenningTrap
{ 
public:
  // Constructor
  double B0_,V0_, d_;
  double prefactor_e_field_;
  arma::vec r;
  arma::vec f_;
  double f_value_;
  double omega_V_value_;
  arma::vec omega_V_;
  arma::vec sum;
  arma::vec coordinate_dependent_vector;
  arma::vec b_field_unit_vector_;

  std::vector<Particle> particles_;
  arma::vec pos_current_;
  arma::vec pos_not_current_; 

  PenningTrap(double B0, double V0, double d);

  // Add a particle to the trap
  void add_particle(Particle new_particle);
  

  //remove all particles
  void remove_all_particles();

  //External electric field at paricle i
  arma::vec external_E_field_on_current_particle(int current_particle,double time);


  // External magnetic field at particle i
  arma::vec external_B_field_on_current_particle(int current_particle);  

  // Force on particle_i from other particles
  arma::vec internal_E_field_on_current_particle(int current_particle);//,std::vector<Particle>& copy_of_particles);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force_with_interactions(int current_particle,double time);
  
  // without interactions
  arma::vec total_force_without_interactions(int current_particle,double time);

 
  //simulate with forward Euler
  void forward_Euler(double time,int time_steps,std::string with_or_without_interactions);
  //simulate with Runge-Kutta 4
  void RK4(double time,int time_steps,std::string with_or_without_interactions,std::string make_files);
  //task 9
  void task9(double d,int charge, double mass,arma::vec f,double time,int time_steps,std::string with_or_without_interactions,arma::vec w1,arma::vec w2,arma::vec w3);
};

#endif