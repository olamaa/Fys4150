#include<penningtrap.hpp>
#include<iostream>
#include<particle.hpp>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iterator>
#include<string>
#include<fstream>
#include<iomanip>
 


  
PenningTrap::PenningTrap(double B0, double V0, double d)
{
    //defining variables
    B0_ = B0;
    V0_ = V0;
    d_ = d;
    prefactor_e_field_= V0_/(pow(d_,2));
    b_field_unit_vector_ = {0,0,1};
    
 
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle new_particle){
    
    particles_.push_back(new_particle);
}

//remove all particles
void PenningTrap::remove_all_particles(){
    
    particles_ = {};
}


//// External electric field on particle i
arma::vec PenningTrap::external_E_field_on_current_particle(int current_particle,double time){
    r = particles_[current_particle].get_position();
    coordinate_dependent_vector = {r(0),r(1),-2*r(2)};
    //checking if amplitudes exist
    if (f_.size() == 0){
        return prefactor_e_field_*coordinate_dependent_vector;
    }
    
    else{
        //checking if the particles are outside the trap
        if (arma::norm(r,2) > d_){
            return {0,0,0};
        }
        else{
            return prefactor_e_field_*(1+f_value_*cos(omega_V_value_*time))*coordinate_dependent_vector;
        }
        
    }
}

//external magnetinc field on particles i
arma::vec PenningTrap::external_B_field_on_current_particle(int current_particle){
    r = particles_[current_particle].get_position();
    //checking if the particles are outside the trap
    if (arma::norm(r,2) > d_){
            return {0,0,0};
        }
    else{
        return b_field_unit_vector_*B0_;
        }
}  

//internal electric field, i.e. particle to particle interactions
arma::vec PenningTrap::internal_E_field_on_current_particle(int current_particle){
    pos_current_ = (particles_)[current_particle].get_position();
    sum = {0,0,0};
    //running through all particles
    for (int not_current_particle = 0; not_current_particle < (particles_).size();not_current_particle++){
        //checking that we do not add the electric field from particle i
        if (not_current_particle != current_particle){
            pos_not_current_ = (particles_)[not_current_particle].get_position();
            
            sum += (particles_)[not_current_particle].get_charge()*(pos_current_-pos_not_current_)/pow(arma::norm((pos_current_-pos_not_current_),2),3);
        }
        
    } 

    return 1.38935333*pow(10,5)*sum;
}

//total force on particle i with interactions
arma::vec PenningTrap::total_force_with_interactions(int current_particle, double time){
return (particles_)[current_particle].get_charge()*(external_E_field_on_current_particle(current_particle,time) + internal_E_field_on_current_particle(current_particle))
    + arma::cross((particles_)[current_particle].get_charge()*(particles_)[current_particle].get_velocity(),external_B_field_on_current_particle(current_particle));
}
//total force on particle i without interactions
arma::vec PenningTrap::total_force_without_interactions(int current_particle,double time){
return (particles_)[current_particle].get_charge()*(external_E_field_on_current_particle(current_particle,time))
    + arma::cross((particles_)[current_particle].get_charge()*(particles_)[current_particle].get_velocity(),external_B_field_on_current_particle(current_particle));
}

//forward Euler method
void PenningTrap::forward_Euler(double time,int time_steps,std::string with_or_without_interactions){

    double dt = time/time_steps;
    int width = 40;
    int prec = 12;
    double current_time;
    //writing the initial conditions to file
    if (with_or_without_interactions == "yes"){
        std::string filename; 
        std::ofstream ofile_;
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            filename = "Forward_Euler_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
            ofile_.open(filename);
            ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << 0 << std::endl; 
            ofile_.close();     
        }

 
        arma::vec a;
        //calculating and writing the updated velocities and positions of the particles 
        for (int time_iteration = 0;time_iteration < time_steps;time_iteration++){
            current_time = dt*time_iteration;
            for (int current_particle_uniq = 0;current_particle_uniq < (particles_).size();current_particle_uniq++){
                a = total_force_with_interactions(current_particle_uniq,current_time)/((particles_)[current_particle_uniq].get_mass());
                (particles_)[current_particle_uniq].change_position(((particles_)[current_particle_uniq].get_position()+((particles_)[current_particle_uniq].get_velocity()*dt)));
                (particles_)[current_particle_uniq].change_velocity((particles_)[current_particle_uniq].get_velocity()+a*dt);
                filename = "Forward_Euler_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
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
        //the same as above only without interactions   
        std::string filename; 
        std::ofstream ofile_;
        
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            filename = "Forward_Euler_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
            ofile_.open(filename);
            ofile_<< std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << std::setprecision(prec) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::setprecision(prec) << 0 << std::endl; 
            ofile_.close();     
        }

        arma::vec a;
        
        for (int time_iteration = 0;time_iteration< time_steps;time_iteration++){
            current_time = dt*time_iteration;
            for (int current_particle_uniq_1 = 0;current_particle_uniq_1 < (particles_).size();current_particle_uniq_1++){
                a = total_force_without_interactions(current_particle_uniq_1,time_iteration*dt)/((particles_)[current_particle_uniq_1].get_mass());
                (particles_)[current_particle_uniq_1].change_velocity((particles_)[current_particle_uniq_1].get_velocity()+a*dt);
                (particles_)[current_particle_uniq_1].change_position(((particles_)[current_particle_uniq_1].get_position()+((particles_)[current_particle_uniq_1].get_velocity()*dt)));
                filename = "Forward_Euler_particle_" + std::to_string(current_particle_uniq_1) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
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

    //creating the k variables
    for (int i=0;i<particles_.size();i++){
    k_all_a.push_back(k_to_be_pushed_back);
    k_all_v.push_back(k_to_be_pushed_back);
    }
    
    //writing initial values to file
    if (with_or_without_interactions == "yes"){
        std::string filename; 
        std::ofstream ofile_;
        if (make_files == "yes"){
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
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
        //running through all times
        for (int time_iteration = 0;time_iteration< time_steps;time_iteration++){
            current_time = dt*time_iteration;

            //getting k_0
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](0,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,current_time)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](0,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }

            //evolving system to step i+1/2    
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                pos(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_position());
                vel(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](0,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](0,arma::span::all)));
            }

            //getting k_1
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](1,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,time_iteration*dt+dt/2.)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](1,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }

            //reverting the system, i.e. changing back the positions and the velocities to step i
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
                //evolving system
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](1,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](1,arma::span::all)));        
                }
            //getting k_2
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](2,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,time_iteration*dt+dt/2.)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](2,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //change back the positions and the velocities
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){
                particles_[current_particle].change_position(arma::trans(pos(current_particle,arma::span::all)));
                particles_[current_particle].change_velocity(arma::trans(vel(current_particle,arma::span::all)));
                }
                //evolve system
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt*k_all_a[current_particle](2,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt*k_all_v[current_particle](2,arma::span::all)));        
                }
                //get k_3
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](3,arma::span::all) = arma::trans(total_force_with_interactions(current_particle,time_iteration*dt+dt)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](3,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            //change back the positions and the velocities
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
                //writes to file
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
        //same as above without interactions
        std::string filename; 
        std::ofstream ofile_;
        if (make_files == "yes"){
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
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
        for (int time_iteration = 0;time_iteration< time_steps;time_iteration++){
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

//task 9
void PenningTrap::task9(double d,int charge, double mass, arma::vec f,double time,int time_steps,std::string with_or_without_interactions,arma::vec w1,arma::vec w2,arma::vec w3){
    f_ = f;
    int number_of_particles_inside_trap;
    std::ofstream ofile;
    std::string filename;
    int width = 40;
    arma::vec rand_pos;
    arma::vec rand_vel;
    //adding different ranges for the frequency analysis
    std::vector<arma::vec> widths;
    widths.push_back(w1);
    widths.push_back(w2);
    widths.push_back(w3);

    //running through all frequencies for all amplitudes
    for (int amplitude = 0;amplitude < f_.size();amplitude++){
            f_value_ = f_(amplitude);
            std::cout <<  "first loop" << std::endl;
            filename = std::to_string(f_value_)+"_"+with_or_without_interactions + "_test.txt";
            ofile.open(filename);
            ofile.close();
            ofile.open(filename, std::ofstream::app);
            omega_V_ = widths[amplitude];
        for (int omega = 0;omega<omega_V_.size();omega++){
            //removes particles i the trap before we make the new particles
            remove_all_particles();
            //creating and adding 100 particles to the trap
            for (int i = 0;i<100;i++){
               arma::arma_rng::set_seed_random();
               rand_pos.randn(3)*0.1*d;
               rand_vel.randn(3)*0.1*d;
               Particle particle_x = Particle(charge,mass,rand_pos,rand_vel);
               add_particle(particle_x);
            }
            number_of_particles_inside_trap = 0;
            omega_V_value_ = omega_V_(omega);
            //running the simulation
            if (with_or_without_interactions=="with"){
                RK4(time,time_steps,"yes","don't make all the files");
            }
            else{
                RK4(time,time_steps,"no","don't make all the files");
            }
            //counting particles
            for (int i = 0;i<particles_.size();i++){
                if (arma::norm(particles_[i].get_position(),2) <=d_){
                    number_of_particles_inside_trap += 1;
                    }
                else{}
            }
            //writing trapped particles to file
            ofile << std::setw(40)<< omega_V_value_ << std::setw(40) << number_of_particles_inside_trap << std::endl;
        }
        ofile.close();    
    }
}
