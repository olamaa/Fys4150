#ifndef __particle_hpp__
#define __particle_hpp__

#include <iostream>
#include <armadillo>

class Particle 
{
public:
   double charge_, mass_;
   arma::vec velocity_,position_;

    //type1 my_variable = some_value;
    //type1 *pointer_to_my_variable = &my_variable;

    //arma::vec position_to_be_pointed_to_;
    //arma::vec *position_ = &position_to_be_pointed_to_;

    //arma::vec velocity_to_be_pointed_to_;
    //arma::vec *velocity_ = &velocity_to_be_pointed_to_;


   //arma::vec position_to_be_referenced;
   //arma::vec& position_ = position_to_be_referenced;
   //
   //arma::vec velocity_to_be_referenced;
   //arma::vec& velocity_ = velocity_to_be_referenced;

   Particle(double charge,double mass, arma::vec position,arma::vec velocity);

   double get_mass();

   double get_charge();

   //void set_position();

   arma::vec get_position();

   arma::vec get_velocity();

   void change_position(arma::vec new_position);

   void change_velocity(arma::vec new_velocity);
};

#endif
