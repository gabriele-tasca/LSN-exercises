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
   int n_throws = 10000; 

   ofstream outfile;
   outfile.open("results.out", ios::trunc);
   outfile << n_blocks << " # number of blocks" << endl;
   outfile << n_throws << " # number of throws in each block" << endl;
   outfile << "# progressive mean of integral # progressive uncertainty " << endl;

   double xmin = 0.;
   double xmax = 1.;

   double sigma = 0.4;
   double mean = 0.;

   // double magic_normalization = 0.313308687321;
   double magic_normalization = 0.4951;

   auto function = [](double x) { return 0.5*M_PI * cos( 0.5*M_PI*x );};

   double outer_sum = 0;
   double outer_squares_sum = 0;
   for(int i=0; i<n_blocks; i++){
      
      double f_average = 0;
      for(int j=0; j<n_throws; j++){
         double x1;
         do{
            x1 = fabs(rnd.Gauss(mean, sigma));
         }while( x1 > 1);

         f_average += function(x1)/( 1./magic_normalization *exp( -0.5*pow((x1 - mean)/(sigma) ,2))   );
         
      }
      f_average = f_average / n_throws;
      double integral = f_average * (xmax - xmin);
      
      outer_sum += integral;
      outer_squares_sum += pow(integral , 2);

      double average_so_far = outer_sum/(i+1);
      double average_of_squares_so_far = outer_squares_sum/(i+1);
      double uncertainty_so_far = sqrt(  (1/(float(i))) * ( average_of_squares_so_far -  pow(average_so_far,2)) );

      outfile << average_so_far << " " << uncertainty_so_far << endl;
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
