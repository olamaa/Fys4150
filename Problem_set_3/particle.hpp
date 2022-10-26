#ifndef __particle_hpp__
#define __particle_hpp__

#include <iostream>
#include <armadillo>

class Particle 
{
public: //defining variables
   double charge_, mass_;
   arma::vec velocity_,position_;

   Particle(double charge,double mass, arma::vec position,arma::vec velocity);
   //functions to get the member variables
   double get_mass();

   double get_charge();

   arma::vec get_position();

   arma::vec get_velocity();

   //to change the positions and velocities

   void change_position(arma::vec new_position);

   void change_velocity(arma::vec new_velocity);
};

#endif
