#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 

	// read simulation parameters as command-line input
	double h = atof(argv[1]);
	double dt = atof(argv[2]);
	double T = atof(argv[3]);
	double x_c = atof(argv[4]);
	double sigma_x = atof(argv[5]);
	double p_x = atof(argv[6]);
	double y_c = atof(argv[7]);
	double sigma_y = atof(argv[8]);
	double p_y = atof(argv[9]);
	double v_0 = atof(argv[10]);

	int M = 1/h;
	int M_2 = M - 2;
	int slit = 2; //number of slits
	mat V = set_pot (M_2, h, slit, v_0);	//set potential matrix
	cx_mat U_0 (M_2, M_2, fill::zeros);
	Gauss_pack (M_2, x_c, y_c, sigma_x, sigma_y, p_x, p_y, U_0);

	sp_mat A_im (M_2*M_2, M_2*M_2);
	sp_mat B_im (M_2*M_2, M_2*M_2);

	//imaginary part of A and B
	fill_AB (M_2, h, dt, V, A_im, B_im);

	//put imaginary and real part together
	//sp_cx_mat A (eye<sp_mat>(M_2*M_2), A_im);
	mat identity = eye(M_2*M_2, M_2*M_2);
	sp_mat identity_sparse = sp_mat(identity);
	sp_cx_mat A (identity_sparse, A_im);
	sp_cx_mat B (identity_sparse, B_im);
	cx_cube U (M_2, M_2, T/dt + 1);
	cx_vec U_n (M_2*M_2);
	cx_vec b (M_2*M_2);

	cx_vec U_0_flat = U_0.as_col();
	U_0_flat = U_0_flat/norm(U_0_flat);
	b = B*U_0_flat;
	U_0 = reshape(U_0_flat, M_2, M_2);
	U.slice(0) = U_0;
	mat probability (T/dt + 1, 2);
	probability(0,0) = real(cdot(U_0_flat, U_0_flat));
	probability(0,1) = 0;

	for (int j = 1; j <= T/dt; j ++) {
		//solve the matrix equation
		U_n = spsolve(A, b);
		b = B*U_n;
		U.slice(j) = reshape(U_n, M_2, M_2);	//fill the cube
		probability(j,0) = real(cdot(U_n, U_n));
		probability(j,1) = j*dt;
	}

	if (v_0 == 0){ 
		probability.save("probability.txt", arma_ascii);
	}
	else{
		probability.save("probability_ds.txt", arma_ascii);
	}
	
	// save as arma_binary format
	U.save("U.bin");

	// save in raw_ascii format
//	U.save("U.txt", raw_ascii);

	//A.print();
	//cout<<endl;
	//B.print();

}