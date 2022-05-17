#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	int L = atoi(argv[1]);		// set dimension of the LxL matrix
	double T = atof(argv[2]);	// set temperature
	arma_rng::set_seed_random();
	
	// to store data in a file
	ofstream myfile;
	// set ordered or random initial state
	string state = argv[3];
	imat spin_matrix;
	if (state == "ordered"){
		spin_matrix = ones<imat>(L,L);
		myfile.open ("5a_ord_"+std::to_string((int)T)+"_"+std::to_string(L)+".txt");
	}
	else{
		spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;
		myfile.open ("5a_rand_"+std::to_string((int)T)+"_"+std::to_string(L)+".txt");
	}
	
	int max_cycles = 1000000;		// number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages	
	
	// initial energy and magnetization
	double E = energy(spin_matrix);
	double M = mag(spin_matrix);
	
	double X;	 // susceptibility
	double Cv; 	 // specific heat capacity 
	
	double X_an, Cv_an, E_an, M_an, E2_an, M2_an, Z_an;	 // analytical values
	
	// set up array for the Boltzmann factor
	vec Bf (17, fill::zeros);
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-1.*i/T);
	}
	
	for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){	// do it for <max_cycles> cycles
		// Monte Carlo 
		Metropolis(spin_matrix, E, M, Bf);
		// update averages
		average(0) += E; 
		average(1) += E*E; 
		average(2) += fabs(M);
		average(3) += M*M;
	
		myfile<<n_cycles<<" "<<setw(15) << setprecision(8) <<average(0)/(double)(n_cycles*L*L)<<" "<<average(2)/(double)(n_cycles*L*L)<<endl;
	}
	
	myfile.close();
	
	if (L==2 && T==1){

		average(0) = average(0)*(1./(max_cycles*L*L));
		average(2) = average(2)*(1./(max_cycles*L*L));
		average(1) = average(1)*(1./(max_cycles*L*L));
		average(3) = average(3)*(1./(max_cycles*L*L));

		// Compute specific heat capacity and susceptibility
		Cv = (average(1) - pow(average(0),2)*(L*L))/(T*T);
		X = (average(3) - pow(average(2),2)*(L*L))/T;

		// Compute analytical values for T=1 and a lattice with L=2
		Z_an = 2*(exp(8)+exp(-8)+6);
		E_an = 4/Z_an*(exp(-8)-exp(8));
		M_an = 2/Z_an*(exp(8)+2);
		E2_an = 32/Z_an*(exp(8)+exp(-8));
		M2_an = 8/Z_an*(exp(8)+1);

		Cv_an = E2_an - pow(E_an, 2)*4;
		X_an = M2_an - pow(M_an,2)*4;

		cout << "energy per spin:" << average(0) << endl;
		cout << "magnetization per spin:" << average(2) << endl;
		cout << "specific heat capacity per spin:" << Cv << endl;
		cout << "susceptibility per spin:" << X << endl;

		cout << "E an-num:" << E_an-average(0) << endl;
		cout << "M an-num:" << M_an-average(2) << endl;
		cout << "Cv an-num:" << Cv_an-Cv << endl;
		cout << "X an-num:" << X_an-X << endl;		

		cout << "E an per spin:" << E_an << endl;
		cout << "M an per spin:" << M_an << endl;
		cout << "Cv an per spin:" << Cv_an << endl;
		cout << "X an per spin:" << X_an<< endl;	
	}

	return 0;
}
