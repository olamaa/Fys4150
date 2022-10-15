#include <particle.hpp>

Particle::Particle(double charge,double mass, arma::vec& position,arma::vec& velocity)
{



    charge_ = charge;
    mass_ = mass;
    position_ = position;
    velocity_ = velocity;
    //position_to_be_pointed_to_ = position;
    //velocity_to_be_pointed_to_ = velocity;
    //position_ *= &position_to_be_pointed_to_;
    //type1 my_variable = some_value;
    //type1 *pointer_to_my_variable = &my_variable;
    //position_ = position_to_be_pointed_to_;
    //position_ = {0,0,0};
    
    //arma::vec position_ = position;
    //arma::vec *pointer_to_position_ = &position_;
    //arma::ved posit


    //arma::vec velocity_ = velocity;
}


//void set_position(arma::vec initial_position){
//    position_ = initial_position;
//}

double Particle::get_mass(){
    return mass_;
}

double Particle::get_charge(){
    return charge_;
}

arma::vec Particle::get_position(){
    return position_;
}

arma::vec Particle::get_velocity(){
    return velocity_;
}

void Particle::change_position(arma::vec new_position){
    position_ = new_position;
    
}

void Particle::change_velocity(arma::vec new_velocity){
    velocity_ = new_velocity;
}