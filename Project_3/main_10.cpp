#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <ctime>


using namespace std;
using namespace arma;

int main (){ 
   
   	
	ofstream myfile_omega;
	myfile_omega.open ("omega.txt");
	
	
//arma_rng::set_seed_random();

double n_step = 1.0e3;
double t_tot = 500.;
double dt = t_tot/n_step;
double n_particles=100.;
double t=0;
vector<int> p_inside;

vec position_0 (3);
vec velocity_0 (3);
mat  r_step (3,n_particles);
mat v_step (3,n_particles);
mat pos_0 (3, n_particles);
mat vel_0 (3, n_particles);
double V_in=0.0025*9.64852558*1.0e7;


//PenningTrap my_trap_inter = PenningTrap(9.65e1, V_in, 500.,1);
PenningTrap my_trap_nointer = PenningTrap(9.65e1, V_in, 500.,0,0);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

for (int k=0; k<n_particles;k++) { 
 	position_0= vec(3).randn() * 0.1 * my_trap_nointer.d_;
 	velocity_0= vec(3).randn() * 0.1 * my_trap_nointer.d_;

	my_particle.position_=position_0;
	my_particle.velocity_=velocity_0;
	//my_trap_inter.add_particle(my_particle);
	my_trap_nointer.add_particle(my_particle);
	
	pos_0.col(k)=position_0;
	vel_0.col(k)=velocity_0;

	//cout<<"initial position particle "<< k+1<<": "<< my_particle.position()<<endl;
}



// no interactions

for (double om =0.2; om<=2.5; om += 0.02){				//loop over omega_v
	my_trap_nointer.omega_v_= 2;
	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_nointer.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				r_step.col(j)=pos_0.col(j);
				v_step.col(j)=vel_0.col(j);
			}
			else {
				t=dt*(i+1);
				my_trap_nointer.evolve_RK4(dt,j, r_step, v_step,t);
			}
		}
	}
	p_inside.push_back(my_trap_nointer.count_particles(r_step)/n_particles);
	myfile_omega<<om<<endl;
}
		

	myfile_omega.close();

	
	ofstream myfile1;
	myfile1.open ("p_inside_0.1.txt");
	for (int i = 0; i < p_inside.size(); i++) {
	myfile1<<p_inside.at(i)<<endl;
	}
	myfile1.close();
	

return 0;

}

/*
if (n_particles>1){

// with interactions 

	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_inter.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				R.slice(j).col(i)=my_trap_inter.particle_collection[j].position();
				V.slice(j).col(i)=my_trap_inter.particle_collection[j].velocity();
				r_step.col(j)=my_trap_inter.particle_collection[j].position();
				v_step.col(j)=my_trap_inter.particle_collection[j].velocity();
			}
			else {
				my_trap_inter.evolve_RK4(dt,j, r_step, v_step);
				//my_trap_inter.evolve_forward_Euler(dt,j, r_step, v_step);
				R.slice(j).col(i)=r_step.col(j);
				V.slice(j).col(i)=v_step.col(j);
			}
		}
		}
	
	ofstream myfile1_inter;
	myfile1_inter.open ("r1_inter.txt");
	myfile1_inter<<R.slice(0).t();
	myfile1_inter.close();

	ofstream myfile2_inter;
	myfile2_inter.open ("r2_inter.txt");
	myfile2_inter<<R.slice(1).t();
	myfile2_inter.close();
	
	ofstream myfile3_inter;
	myfile3_inter.open ("v1_inter.txt");
	myfile3_inter<<V.slice(0).t();
	myfile3_inter.close();
	
	ofstream myfile4_inter;
	myfile4_inter.open ("v2_inter.txt");
	myfile4_inter<<V.slice(1).t();
	myfile4_inter.close();
	
	cout<< "interaction "<< my_trap_inter.interaction_<<endl;	
	cout<< "force inter total "<< my_trap_inter.total_force(0,0)<<endl;
//	cout<< "force inter particle "<< my_trap_inter.force_particle(0,1)<<endl;
}
*/