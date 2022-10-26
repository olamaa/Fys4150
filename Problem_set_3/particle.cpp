#include <particle.hpp>

Particle::Particle(double charge,double mass, arma::vec position,arma::vec velocity)
{   //define vriables
    charge_ = charge;
    mass_ = mass;
    position_ = position;
    velocity_ = velocity;
}

//unnecessary functions to get the member variabels...
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
//to change the positions and velocities
void Particle::change_position(arma::vec new_position){
    position_ = new_position;
    
}

void Particle::change_velocity(arma::vec new_velocity){
    velocity_ = new_velocity;
}