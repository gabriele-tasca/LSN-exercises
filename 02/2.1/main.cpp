/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   ///////////////////////

   int n_blocks = 10000;
   int n_steps = 100;

   ofstream outfile;
   outfile.open("results.out", ios::trunc);
   outfile << n_blocks << " # number of blocks" << endl;
   outfile << n_steps << " # number of throws in each block" << endl;
   outfile << "# number of steps # spread # uncertainty " << endl;

   ofstream posfile;
   posfile.open("pos.out", ios::trunc);

   vector<double> spread_at_step(n_steps, 0);

   vector<double> outer_sum(n_steps, 0);
   vector<double> outer_squares_sum(n_steps, 0);
   
   // vector<double> average_so_far(n_steps, 0);
   // vector<double> average_of_squares_so_far(n_steps, 0);
   // vector<double> uncertainty_so_far(n_steps, 0);
   
   for(int i=0; i<n_blocks; i++){

      fill(spread_at_step.begin(), spread_at_step.end(), 0); //very idiomatic! no memset here
      
      double x = 0;
      double y = 0;
      double z = 0;
      for(int j=0; j<n_steps; j++){
         double vers = 1.;
         if(rnd.Rannyu() > 0.5){
            vers = -1.;
         }
         double direct = rnd.Rannyu();
         if(direct < 0.3333333333333333){
            x += vers;
         }
         else if(direct < 0.6666666666666666){
            y += vers;
         }
         else if(direct < 0.9999999999999999){
            z += vers;
         }
         else{
            continue; // just in case
         }
         spread_at_step[j] =  pow(x,2) + pow(y,2) + pow(z,2) ;
         
         //end step loop
      }      

      for(int kappa = 0; kappa < n_steps; kappa++){
         spread_at_step[kappa] = sqrt(spread_at_step[kappa]);
         outer_sum[kappa] += spread_at_step[kappa];
         outer_squares_sum[kappa] += pow(spread_at_step[kappa],2); 

      }
      //end block loop
   }

      for(int kappa = 0; kappa < n_steps; kappa++){

         double final_average = outer_sum[kappa]/n_blocks;
         double final_average_of_squares = outer_squares_sum[kappa]/n_blocks;
         double final_uncertainty = sqrt(  1/double(n_blocks-1) * ( final_average_of_squares -  pow(final_average,2) ) );

         outfile << kappa+1 << " " << final_average << " "  << final_uncertainty << endl;
      }


   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
