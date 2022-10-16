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
 


//if test to check if there is only one particle is important. Can then compare with analyticl solution... maybe inside the time solver
  
PenningTrap::PenningTrap(double B0, double V0, double d)
{
    B0_ = B0;
    V0_ = V0;
    d_ = d;
    prefactor_e_field_= V0_/(pow(d_,2));
    b_field_unit_vector_ = {0,0,1};
 
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle& new_particle){
    
    particles_.push_back(new_particle);
    //arma::mat k1 = {{1,0,0},{2,0,0},{3,0,0},{4,0,0}};
    //arma::mat k2 = {{5,0,0},{6,0,0},{7,0,0},{8,0,0}};
    ////std::vector<arma::mat> k_s;
    //k_all.push_back(k1);
    //k_all.push_back(k2);


}


//// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field_on_current_particle(int current_particle){
    //std::cout << "external_e_field" << std::endl;
    r = particles_[current_particle].get_position();
    //std::cout << "position_in_external_E_field" << std::endl;
    //std::cout << r << std::endl;
    coordinate_dependent_vector = {r(0),r(1),-2*r(2)};
    //std::cout << "external_e_field_end" << std::endl;
    return prefactor_e_field_*coordinate_dependent_vector;
}



arma::vec PenningTrap::external_B_field_on_current_particle(){
    return b_field_unit_vector_*B0_;
}
  

arma::vec PenningTrap::internal_E_field_on_current_particle(int current_particle){

    sum = {0,0,0};
    for (int not_current_particle = 0; not_current_particle < (particles_).size();not_current_particle++){
        
        
        if (not_current_particle != current_particle){
            pos_current_ = (particles_)[current_particle].get_position();
            pos_not_current_ = (particles_)[not_current_particle].get_position();
            
            sum += (particles_)[not_current_particle].get_charge()*(pos_current_-pos_not_current_)/pow(arma::norm((pos_current_-pos_not_current_),2),3);
        }
        
    } 

    return 1.38935333*pow(10,5)*sum;
}

arma::vec PenningTrap::total_force_with_interactions(int current_particle){
return (particles_)[current_particle].get_charge()*(external_E_field_on_current_particle(current_particle) + internal_E_field_on_current_particle(current_particle))
    + arma::cross((particles_)[current_particle].get_charge()*(particles_)[current_particle].get_velocity(),external_B_field_on_current_particle());
}

arma::vec PenningTrap::total_force_without_interactions(int current_particle){
return (particles_)[current_particle].get_charge()*(external_E_field_on_current_particle(current_particle))
    + arma::cross((particles_)[current_particle].get_charge()*(particles_)[current_particle].get_velocity(),external_B_field_on_current_particle());
}

//euler-cromer
void PenningTrap::forward_Euler(double time,int time_steps,std::string with_or_without_interactions){

    double dt = time/time_steps;
    int width = 40;

    if (with_or_without_interactions == "yes"){
        std::string filename; 
        std::ofstream ofile_;
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";        
            ofile_.open(filename);
            ofile_<< std::setw(width) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << 0 << std::endl; 
            ofile_.close();     
        }

 
        arma::vec a;
        
        for (int time_iteration = 1;time_iteration<= time_steps;time_iteration++){
            for (int current_particle_uniq = 0;current_particle_uniq < (particles_).size();current_particle_uniq++){
                a = total_force_with_interactions(current_particle_uniq)/((particles_)[current_particle_uniq].get_mass());
                (particles_)[current_particle_uniq].change_velocity((particles_)[current_particle_uniq].get_velocity()+a*dt);
                (particles_)[current_particle_uniq].change_position(((particles_)[current_particle_uniq].get_position()+((particles_)[current_particle_uniq].get_velocity()*dt)));
                //filename = "particle_" + std::to_string(current_particle_uniq_1)+ "_with_interactions.txt";             
                filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << (particles_)[current_particle_uniq].get_position()(0) <<
                std::setw(width) << (particles_)[current_particle_uniq].get_position()(1) << 
                std::setw(width) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(0) <<
                std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(1) << 
                std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << std::to_string(dt*time_iteration) << std::endl; 
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
            ofile_<< std::setw(width) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << 0 << std::endl; 
            ofile_.close();     
        }

        //std::cout << "line_205" << std::endl;
        arma::vec a;
        
        for (int time_iteration = 1;time_iteration<= time_steps;time_iteration++){
            //std::cout << "time_loop" << std::endl;
            for (int current_particle_uniq_1 = 0;current_particle_uniq_1 < (particles_).size();current_particle_uniq_1++){
                //std::cout << "round 2??" << std::endl;
                //std::cout << current_particle_uniq_1 << std::endl;
                a = total_force_without_interactions(current_particle_uniq_1)/((particles_)[current_particle_uniq_1].get_mass());
                //std::cout << a << std::endl;
                //std::cout << "line_212" << std::endl;
                //std::cout << "velocity_before" << std::endl;
                //std::cout << (particles_)[current_particle_uniq_1].get_velocity() << std::endl;
                (particles_)[current_particle_uniq_1].change_velocity((particles_)[current_particle_uniq_1].get_velocity()+a*dt);
                //std::cout << "velocity_after" << std::endl;
                //std::cout << (particles_)[current_particle_uniq_1].get_velocity() << std::endl;
                //std::cout << "line_254" << std::endl;
                (particles_)[current_particle_uniq_1].change_position(((particles_)[current_particle_uniq_1].get_position()+((particles_)[current_particle_uniq_1].get_velocity()*dt)));
                std::cout << "line_257" << std::endl;
                //filename = "particle_" + std::to_string(current_particle_uniq_1)+ "_without_interactions.txt";
                filename = "Euler_Cromer_particle_" + std::to_string(current_particle_uniq_1) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << (particles_)[current_particle_uniq_1].get_position()(0) <<
                std::setw(width) << (particles_)[current_particle_uniq_1].get_position()(1) << 
                std::setw(width) << (particles_)[current_particle_uniq_1].get_position()(2) << std::setw(width) << (particles_)[current_particle_uniq_1].get_velocity()(0) <<
                std::setw(width) << (particles_)[current_particle_uniq_1].get_velocity()(1) << 
                std::setw(width) << (particles_)[current_particle_uniq_1].get_velocity()(2) << std::setw(width) << std::to_string(dt*time_iteration) << std::endl; 
                ofile_.close();     
                }
            }

        }

}

//runge-kutta 4
void PenningTrap::RK4(double time,int time_steps,std::string with_or_without_interactions){
    std::cout << "line 173" << std::endl;
    double dt = time/time_steps;
    int width = 40;
    std::vector<arma::mat> k_all_a;
    std::vector<arma::mat> k_all_v;
    int k_n = 4;
    arma::mat k_to_be_pushed_back = arma::zeros(k_n,3);
    std::cout << "line 180" << std::endl;


    
    for (int i=0;i<particles_.size();i++){
    k_all_a.push_back(k_to_be_pushed_back);
    k_all_v.push_back(k_to_be_pushed_back);
    }
    std::cout << "line 188" << std::endl;


    //k_all.push_back(k2);

    //arma::mat k_a = arma::zeros(particles_.size(),4);
    //arma::mat k_v = arma::zeros(particles_.size(),4);

    

    if (with_or_without_interactions == "yes"){
        std::string filename; 
        std::ofstream ofile_;
        std::cout << "line 201" << std::endl;
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            //filename = "particle_" + std::to_string(current_particle_uniq) + "_with_interactions_RK4.txt";        
            filename = "RK4_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
            ofile_.open(filename);
            ofile_<< 
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(2) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << 0 << std::endl; 
            ofile_.close();        
            }

        arma::mat pos = arma::zeros(particles_.size(),3);
        arma::mat vel = arma::zeros(particles_.size(),3);
        arma::mat update_pos;
        arma::mat update_vel;
        for (int time_iteration = 1;time_iteration<= time_steps;time_iteration++){


            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](0,arma::span::all) = arma::trans(total_force_with_interactions(current_particle)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](0,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                pos(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_position());
                vel(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](0,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](0,arma::span::all)));
            }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](1,arma::span::all) = arma::trans(total_force_with_interactions(current_particle)/((particles_)[current_particle].get_mass()));
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
                k_all_a[current_particle](2,arma::span::all) = arma::trans(total_force_with_interactions(current_particle)/((particles_)[current_particle].get_mass()));
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
                k_all_a[current_particle](3,arma::span::all) = arma::trans(total_force_with_interactions(current_particle)/((particles_)[current_particle].get_mass()));
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
                filename = "RK4_particle_" + std::to_string(current_particle) + "_n_steps_" + std::to_string(time_steps)+ "_with_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << (particles_)[current_particle].get_position()(0) <<
                std::setw(width) << (particles_)[current_particle].get_position()(1) << 
                std::setw(width) << (particles_)[current_particle].get_position()(2) << std::setw(width) << (particles_)[current_particle].get_velocity()(0) <<
                std::setw(width) << (particles_)[current_particle].get_velocity()(1) << 
                std::setw(width) << (particles_)[current_particle].get_velocity()(2) << std::setw(width) << std::to_string(dt*time_iteration) << std::endl; 
                ofile_.close();       
                }
    
            }
        }

    else{
        std::string filename; 
        std::ofstream ofile_;
        for (int current_particle_uniq = 0; current_particle_uniq < (particles_).size(); current_particle_uniq++){
            //filename = "particle_" + std::to_string(current_particle_uniq) + "_without_interactions_RK4.txt";        
            filename = "RK4_particle_" + std::to_string(current_particle_uniq) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
            ofile_.open(filename);
            ofile_<< std::setw(width) << (particles_)[current_particle_uniq].get_position()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_position()(2) << std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(0) <<
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(1) << 
            std::setw(width) << (particles_)[current_particle_uniq].get_velocity()(2) << std::setw(width) << 0 << std::endl; 
            ofile_.close();         
            }

        arma::mat pos = arma::zeros(particles_.size(),3);
        arma::mat vel = arma::zeros(particles_.size(),3);
        arma::mat update_pos;
        arma::mat update_vel;
        for (int time_iteration = 1;time_iteration<= time_steps;time_iteration++){


            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                std::cout << "inside loop" << std::endl;
                k_all_a[current_particle](0,arma::span::all) = arma::trans(total_force_without_interactions(current_particle)/((particles_)[current_particle].get_mass()));
                k_all_v[current_particle](0,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                pos(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_position());
                vel(current_particle,arma::span::all) = arma::trans(particles_[current_particle].get_velocity());
                particles_[current_particle].change_velocity(particles_[current_particle].get_velocity()+arma::trans(dt/2*k_all_a[current_particle](0,arma::span::all)));
                particles_[current_particle].change_position(particles_[current_particle].get_position()+arma::trans(dt/2*k_all_v[current_particle](0,arma::span::all)));
            }
            for (int current_particle = 0;current_particle < particles_.size(); current_particle++){ 
                k_all_a[current_particle](1,arma::span::all) = arma::trans(total_force_without_interactions(current_particle)/((particles_)[current_particle].get_mass()));
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
                k_all_a[current_particle](2,arma::span::all) = arma::trans(total_force_without_interactions(current_particle)/((particles_)[current_particle].get_mass()));
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
                k_all_a[current_particle](3,arma::span::all) = arma::trans(total_force_without_interactions(current_particle)/((particles_)[current_particle].get_mass()));
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
                filename = "RK4_particle_" + std::to_string(current_particle) + "_n_steps_" + std::to_string(time_steps)+ "_without_interactions.txt";
                ofile_.open(filename, std::ofstream::app);
                ofile_<< std::setw(width) << (particles_)[current_particle].get_position()(0) <<
                std::setw(width) << (particles_)[current_particle].get_position()(1) << 
                std::setw(width) << (particles_)[current_particle].get_position()(2) << std::setw(width) << (particles_)[current_particle].get_velocity()(0) <<
                std::setw(width) << (particles_)[current_particle].get_velocity()(1) << 
                std::setw(width) << (particles_)[current_particle].get_velocity()(2) << std::setw(width) << std::to_string(dt*time_iteration) << std::endl; 
                ofile_.close();         
                }
    
            }
        }
    
}