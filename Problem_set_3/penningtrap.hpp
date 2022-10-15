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
  arma::vec sum;
  arma::vec coordinate_dependent_vector;
  arma::vec b_field_unit_vector_;
  //std::vector<Particle> particles_to_be_referenced_;
  //std::vector<Particle>& particles_ = particles_to_be_referenced_;
  std::vector<Particle> particles_;
  arma::vec pos_current_;
  arma::vec pos_not_current_; 
  //arma::vec sum_to_be_pointed_to;
  //arma::vec *sum_ = &sum_to_be_pointed_to;
  //std::vector<arma::vec> positions_to_be_referenced_; 
  //std::vector<arma::vec>& positions_ = positions_to_be_referenced_; 

  //std::vector<std::ofstream> ofile_;
  //std::string filename_to_be_pointed_to;
  //std::string *filename_ = &filename_to_be_pointed_to;

  //arma::vec a_to_be_pointed_to;
  //arma::vec *a = &a_to_be_pointed_to;



  //std::vector<std::string> filenames_;
  
 //std::ofstream ofile_to_be_pointed_to;
 //std::ofstream *ofile_ = &ofile_to_be_pointed_to;
  
  
  
  //std::vector<Particle> particles_to_be_pointed_to;
  //std::vector<Particle> *particles_ = &particles_to_be_pointed_to ;
  //std::vector<Particle> copy_of_particles_;

  //particles_ = std::vector<Particle>;
  PenningTrap(double B0, double V0, double d);

  //// Add a particle to the trap
  void add_particle(Particle& new_particle);
//
  //External electric field at point r=(x,y,z)
  arma::vec external_E_field_on_current_particle(int current_particle);

//
  //// External magnetic field at point r=(x,y,z)
  arma::vec external_B_field_on_current_particle();  
//
  //// Force on particle_i from particle_j
  arma::vec internal_E_field_on_current_particle(int current_particle);//,std::vector<Particle>& copy_of_particles);
//
  //// The total force on particle_i from the external fields
  //arma::vec total_force_external(int i);
//
  //// The total force on particle_i from the other particles
  //arma::vec total_force_particles(int i);
//
  //// The total force on particle_i from both external fields and other particles
  arma::vec total_force_with_interactions(int current_particle);

  arma::vec total_force_without_interactions(int current_particle);
//
  //// Evolve the system one time step (dt) using Runge-Kutta 4th order
  //void evolve_RK4(double dt);
//
  //// Evolve the system one time step (dt) using Forward Euler
  void forward_Euler(double time,int time_steps,std::string with_or_without_interactions);

  void RK4(double time,int time_steps,std::string with_or_without_interactions);
};

#endif