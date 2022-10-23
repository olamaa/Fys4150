#include<penningtrap_beehive.hpp>
#include<iostream>
#include<particle.hpp>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iterator>
#include<string>
#include<fstream>
#include<iomanip>
 


//if test to check if there is only one particle is important. Can then compare with analyticl solution... maybe inside the time solver
  
PenningTrap::PenningTrap(double B0, double V0, double d)
{
    B0_ = B0;
    V0_ = V0;
    d_ = d;
    prefactor_e_field_= V0_/(pow(d_,2));
    b_field_unit_vector_ = {0,0,1};
    //f_ = 1;
 
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle new_particle){
    
    particles_.push_back(new_particle);
}

void PenningTrap::remove_all_particles(){
    
    particles_ = {};
}


//// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field_on_current_particle(int current_particle,double time){
    r = particles_[current_particle].get_position();
    coordinate_dependent_vector = {r(0),r(1),-2*r(2)};
    //if (f_.size() == 0){
    //r = particles_[current_particle].get_position();
    //coordinate_dependent_vector = {r(0),r(1),-2*r(2)};
      //  return prefactor_e_field_*coordinate_dependent_vector;
    
    
    //else{
        if (arma::norm(r,2) > d_){
            return {0,0,0};
        }
        else{
            return prefactor_e_field_*(1+f_value_*cos(omega_V_value_*time))*coordinate_dependent_vector;
        }
        
    }



arma::vec PenningTrap::external_B_field_on_current_particle(int current_particle){
    r = particles_[current_particle].get_position();
    
    if (arma::norm(r,2) > d_){
            return {0,0,0};
        }
    else{
        return b_field_unit_vector_*B0_;
        }
}  

arma::vec PenningTrap::internal_E_field_on_current_particle(int current_particle){
    pos_current_ = (particles_)[current_particle].get_position();
    sum = {0,0,0};
    for (int not_current_particle = 0; not_current_particle < (particles_).size();not_current_particle++){
        
        
        if (not_current_particle != current_particle){
            //pos_current_ = (particles_)[current_particle].get_position();
            pos_not_current_ = (particles_)[not_current_particle].get_position();
            
            sum += (particles_)[not_current_particle].get_charge()*(pos_current_-pos_not_current_)/pow(arma::norm((pos_current_-pos_not_current_),2),3);
        }
        
    } 

    return 1.38935333*pow(10,5)*sum;
}

arma::vec PenningTrap::total_force_with_interactions(int current_particle, double time){
return (particles_)[current_particle].get_charge()*(external_E_field_on_current_particle(current_particle,time) + internal_E_field_on_current_particle(current_particle))
    + arma::cross((particles_)[current_particle].get_charge()*(particles_)[current_particle].get_velocity(),external_B_field_on_current_particle(current_particle));
}

arma::vec PenningTrap::total_force_without_interactions(int current_particle,double time){
return (particles_)[current_particle].get_charge()*(external_E_field_on_current_particle(current_particle,time))
    + arma::cross((particles_)[current_particle].get_charge()*(particles_)[current_particle].get_velocity(),external_B_field_on_current_particle(current_particle));
}

//euler-cromer
void PenningTrap::forward_Euler(double time,int time_steps,std::string with_or_without_interactions){

    double dt = time/time_steps;
    int width = 40;
    int prec = 12;
    double current_time;

    if (with_or_without_interactions == "yes"){
        std::string filename; 
        std::ofstream ofile_;
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";        
            ofile_.open(filename);
            ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << 0 << std::endl; 
            ofile_.close();     
        }

 
        arma::vec a;
        
        for (int time_iteration = 0;time_iteration < time_steps;time_iteration++){
            current_time = dt*time_iteration;
            for (int current_particle_uniq = 0;current_particle_uniq < (particles_).size();current_particle_uniq++){
                a = total_force_with_interactions(current_particle_uniq,current_time)/((particles_)[current_particle_uniq].get_mass());
                (particles_)[current_particle_uniq].change_position(((particles_)[current_particle_uniq].get_position()+((particles_)[current_particle_uniq].get_velocity()*dt)));
                (particles_)[current_particle_uniq].change_velocity((particles_)[current_particle_uniq].get_velocity()+a*dt);
                //filename = "particle_" + std::to_string(current_particle_uniq_1)+ "_with_interactions.txt";             
                filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << current_time + dt << std::endl; 
                ofile_.close();     
                }
            }
  
        }
    else{

        //std::cout << "after_if" << std::endl;
   
        std::string filename; 
        std::ofstream ofile_;
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            //filename = "particle_" + std::to_string(current_particle_uniq) + "_without_interactions.txt";           
            filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
            ofile_.open(filename);
            ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << 0 << std::endl; 
            ofile_.close();     
        }

        //std::cout << "line_205" << std::endl;
        arma::vec a;
        
        for (int time_iteration = 0;time_iteration< time_steps;time_iteration++){
            current_time = dt*time_iteration;
            //std::cout << "time_loop" << std::endl;
            for (int current_particle_uniq_1 = 0;current_particle_uniq_1 < (particles_).size();current_particle_uniq_1++){
                //std::cout << "round 2??" << std::endl;
                //std::cout << current_particle_uniq_1 << std::endl;
                a = total_force_without_interactions(current_particle_uniq_1,time_iteration*dt)/((particles_)[current_particle_uniq_1].get_mass());
                //std::cout << a << std::endl;
                //std::cout << "line_212" << std::endl;
                //std::cout << "velocity_before" << std::endl;
                //std::cout << (particles_)[current_particle_uniq_1].get_velocity() << std::endl;
                (particles_)[current_particle_uniq_1].change_velocity((particles_)[current_particle_uniq_1].get_velocity()+a*dt);
                //std::cout << "velocity_after" << std::endl;
                //std::cout << (particles_)[current_particle_uniq_1].get_velocity() << std::endl;
                //std::cout << "line_254" << std::endl;
                (particles_)[current_particle_uniq_1].change_position(((particles_)[current_particle_uniq_1].get_position()+((particles_)[current_particle_uniq_1].get_velocity()*dt)));
                //filename = "particle_" + std::to_string(current_particle_uniq_1)+ "_without_interactions.txt";
                filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq_1) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq_1].get_position()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq_1].get_position()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq_1].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq_1].get_velocity()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq_1].get_velocity()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq_1].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << current_time +dt << std::endl; 
                ofile_.close();     
                }
            }

        }

}

//runge-kutta 4
void PenningTrap::RK4(double time,int time_steps,std::string with_or_without_interactions,std::string make_files){
    double dt = time/time_steps;
    int width = 40;
    int prec = 12;
    std::vector<arma::mat> k_all_a;
    std::vector<arma::mat> k_all_v;
    int k_n = 4;
    arma::mat k_to_be_pushed_back = arma::zeros(k_n,3);
    double current_time;

    
    for (int i=0;i<particles_.size();i++){
    k_all_a.push_back(k_to_be_pushed_back);
    k_all_v.push_back(k_to_be_pushed_back);
    }


    //k_all.push_back(k2);

    //arma::mat k_a = arma::zeros(particles_.size(),4);
    //arma::mat k_v = arma::zeros(particles_.size(),4);

    

    if (with_or_without_interactions == "yes"){
        std::string filename; 
        std::ofstream ofile_;
        if (make_files == "yes"){
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            //filename = "particle_" + std::to_string(current_particle_uniq) + "_with_interactions_RK4.txt";        
            filename = "RK4_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
            ofile_.open(filename);
            ofile_<< 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << 0 << std::endl; 
            ofile_.close();        
            }
        }   
        arma::mat pos = arma::zeros(particles_.size(),3);
        arma::mat vel = arma::zeros(particles_.size(),3);
        arma::mat update_pos;
        arma::mat update_vel;
        for (int time_iteration = 0;time_iteration< time_steps;time_iteration++){
            current_time = dt*time_iteration;


            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](0,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,current_time)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](0,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                pos(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_position());
                vel(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](0,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](0,arma::span::all)));
            }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](1,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,time_iteration*dt+dt/2.)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](1,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //changes back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](1,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](1,arma::span::all)));        
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](2,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,time_iteration*dt+dt/2.)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](2,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //changes back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt*k_all_a[current_particle](2,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt*k_all_v[current_particle](2,arma::span::all)));        
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](3,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,time_iteration*dt+dt)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](3,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //changes back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
            //finally updates the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                update_vel = arma::trans(particles_[current_particle].get_velocity())+dt/6*(k_all_a[current_particle](0,arma::span::all)+
                2*k_all_a[current_particle](1,arma::span::all)+2*k_all_a[current_particle](2,arma::span::all)+k_all_a[current_particle](3,arma::span::all));

                update_pos = arma::trans(particles_[current_particle].get_position())+dt/6*(k_all_v[current_particle](0,arma::span::all)+
                2*k_all_v[current_particle](1,arma::span::all)+2*k_all_v[current_particle](2,arma::span::all)+k_all_v[current_particle](3,arma::span::all));

                particles_[current_particle].change_velocity(arma::trans(update_vel));
                particles_[current_particle].change_position(arma::trans(update_pos));        
                //filename = "particle_" + std::to_string(current_particle) + "_with_interactions_RK4.txt";           
                if (make_files == "yes"){
                filename = "RK4_particle_" + std::to_string(current_particle) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_position()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_position()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_velocity()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_velocity()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << current_time+dt << std::endl; 
                ofile_.close();
                }       
                }
    
            }
        }

    else{
        std::string filename; 
        std::ofstream ofile_;
        if (make_files == "yes"){
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            //filename = "particle_" + std::to_string(current_particle_uniq) + "_without_interactions_RK4.txt";        
            filename = "RK4_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
            ofile_.open(filename);
            ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << 0 << std::endl; 
            ofile_.close();         
            }
        }

        arma::mat pos = arma::zeros(particles_.size(),3);
        arma::mat vel = arma::zeros(particles_.size(),3);
        arma::mat update_pos;
        arma::mat update_vel;
        for (int time_iteration = 0;time_iteration<= time_steps;time_iteration++){
            current_time = dt*time_iteration;

            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](0,arma::span::all) = arma::trans(total_force_without_interactions(current_particle,time_iteration*dt)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](0,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                pos(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_position());
                vel(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](0,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](0,arma::span::all)));
            }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](1,arma::span::all) = arma::trans(total_force_without_interactions(current_particle,time_iteration*dt+dt/2.)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](1,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //changes back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](1,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](1,arma::span::all)));        
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](2,arma::span::all) = arma::trans(total_force_without_interactions(current_particle,time_iteration*dt+dt/2.)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](2,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //changes back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt*k_all_a[current_particle](2,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt*k_all_v[current_particle](2,arma::span::all)));        
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](3,arma::span::all) = arma::trans(total_force_without_interactions(current_particle,time_iteration*dt+dt)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](3,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //changes back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
            //finally updates the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                update_vel = arma::trans(particles_[current_particle].get_velocity())+dt/6*(k_all_a[current_particle](0,arma::span::all)+
                2*k_all_a[current_particle](1,arma::span::all)+2*k_all_a[current_particle](2,arma::span::all)+k_all_a[current_particle](3,arma::span::all));

                update_pos = arma::trans(particles_[current_particle].get_position())+dt/6*(k_all_v[current_particle](0,arma::span::all)+
                2*k_all_v[current_particle](1,arma::span::all)+2*k_all_v[current_particle](2,arma::span::all)+k_all_v[current_particle](3,arma::span::all));

                particles_[current_particle].change_velocity(arma::trans(update_vel));
                particles_[current_particle].change_position(arma::trans(update_pos));        
                //filename = "particle_" + std::to_string(current_particle) + "_without_interactions_RK4.txt";           
                if (make_files == "yes"){
                filename = "RK4_particle_" + std::to_string(current_particle) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_position()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_position()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_velocity()(0) <<
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_velocity()(1) << 
                std::setw(width) << std::setprecision(prec) << (particles_)[current_particle].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << current_time+dt << std::endl; 
                ofile_.close();         
                }
                }   
    
            }
        }
    
}


void PenningTrap::task9(double d,int charge, double mass, double f,double time,int time_steps,std::string with_or_without_interactions,arma::vec w1){
    f_ = f;
    //omega_V_ = omega_V;
    int number_of_particles_inside_trap;
    std::ofstream ofile;
    std::string filename;
    int width = 40;
    arma::vec rand_pos;
    arma::vec rand_vel;
    arma::vec widths = w1;
    //widths.push_back(w2);
    //widths.push_back(w3);
    
    //for (int amplitude = 0;amplitude < f_.size();amplitude++){
            f_value_ = f_;
            std::cout <<  "first loop" << std::endl;
            filename = std::to_string(f_value_)+"_"+with_or_without_interactions + ".txt";
            ofile.open(filename);
            ofile.close();
            ofile.open(filename, std::ofstream::app);
            omega_V_ = widths;
        for (int omega = 0;omega<omega_V_.size();omega++){
            remove_all_particles();
            
            //Particle particle_x = Particle(q,mass_Ca_II,rand_pos,rand_vel);
            //std::vector<Particle> all_particles;
            for (int i = 0;i<100;i++){
                //std::cout <<  "441" << std::endl;
               //all_particles.push_back();
               arma::arma_rng::set_seed_random();
               //std::cout <<  "444" << std::endl;
               rand_pos.randn(3)*0.1*d;
               //std::cout <<  "446" << std::endl;
               rand_vel.randn(3)*0.1*d;
               //std::cout <<  "448" << std::endl;
               Particle particle_x = Particle(charge,mass,rand_pos,rand_vel);
               //std::cout <<  "450" << std::endl;
               add_particle(particle_x);
               //std::cout <<  "452" << std::endl;
            }
            number_of_particles_inside_trap = 0;
            omega_V_value_ = omega_V_(omega);
            std::cout << omega_V_value_  << std::endl;
            RK4(time,time_steps,"yes","don't make all the files");
            
            
            for (int i = 0;i<particles_.size();i++){
                if (arma::norm(particles_[i].get_position(),2) <=d_){
                    number_of_particles_inside_trap += 1;
                    }
                else{}
            }

            //for loop to change initial values of all particles back to original
            //filename = std::to_string(f_value_) + ".txt";
            //ofile.open(filename);
            //ofile << std::endl;
            //ofile.close();
            //ofile.open(filename, std::ofstream::app);
            //for each amplitude: write to file  omega number_of_particles_inside_trap
            ofile << std::setw(40)<< omega_V_value_ << std::setw(40) << number_of_particles_inside_trap << std::endl;
        }
        ofile.close();
            

        //}

    //}
}
