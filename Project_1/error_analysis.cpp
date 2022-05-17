// This program computes an approximation to the second derivative of
// u(x)  using different stepsizes h

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

struct solution {
             float x;
             float v;
            };
            
            
// Count rows
int flcount(ifstream& in){
                  string s;
                  int k=0;
                  getline(in,s);
                  while (!in.eof()) {
                      k++;
                      getline(in,s);
                  }

                  in.clear();
                  in.seekg(0);
return k;
                }     

int main()
{

ifstream in;
        in.open("solutions_7_10^2pts.txt");

         if ( in.fail() ) {
            cout << "Attention! problems with file opening "<< endl;
         return 1; 
         }

// Caricamento array e stampa risultati

            int k = flcount(in);
            solution *S = new solution [k];
        
            for(int i=0; i < k; i++){ 
                in>>S[i].x;
                in>>S[i].v;
             }
                                     
             //cout<<"I dati caricati nell'array sono: "<<endl;

             //for(int i=0; i<k; i++){
             //   cout<<S[i].x<<"\t"<<S[i].v<<endl;
             //}
             //cout<<endl;


  int n = 100;  // n number of points in the grid

  // Parameters for output formatting
  int width = 18;
  int prec  = 10;

  // Stepsize
  //double h = 0.090909;

  // Output a header 
  std::cout << "#" << std::setw(width-1) << "function u"
            //<< std::setw(width) << "d2u_approx"
            //<< std::setw(width) << "d2u_exact"
            << std::setw(width) << "abs_error"
            << std::setw(width) << "rel_error"
            << std::setw(width) << "log10(abs_error)"
            << std::setw(width) << "log10(rel_error)"
            << std::endl;

// Defining empty arrays
  double u[n];
  double abs_err[n];
  double rel_err[n];
  double log_abs_err[n];
  double log_rel_err[n];


	for (int i=0; i<n; i++){
// Function: u(x) = 1-(1-exp(-10))*x-exp(-10x)
		u[i] = 1-(1-exp(-10))*S[i].x-exp(-10.*S[i].x); 

// Absolute error
		abs_err[i] = fabs((u[i] - S[i].v));

// Relative error
		rel_err[i] = fabs((u[i] - S[i].v) / S[i].v);

// log10 absolute error
		log_abs_err[i] = log10(fabs(u[i] - S[i].v));

// log10 relative error
		log_rel_err[i] = log10(fabs((u[i] - S[i].v) / S[i].v));


// Output to screen
	std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << u[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << abs_err[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << rel_err[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << log_abs_err[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << log_rel_err[i]
              << std::endl;
              
              
// Writing the solution on a file
//      ofstream myfile;
//      myfile.open ("err_analysis_10.txt");
 
//    if (!myfile ) { // file couldn't be opened
//      cerr << "Error: file could not be opened" << endl;
//      exit(1);
//   }

//	for (int i=0; i<n; i++) {
//		myfile << abs_err[i]<<" "<< rel_err[i] <<" "<< log_rel_err[i] << "\n"; 
//	}
 
//      myfile.close();
	}

  return 0;
}
