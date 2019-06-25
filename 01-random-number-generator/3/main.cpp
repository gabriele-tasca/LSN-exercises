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

   int n_blocks = 1000;
   int n_throws = 50000; 

   ofstream outfile;
   outfile.open("results.out", ios::trunc);
   outfile << n_blocks << " # number of blocks" << endl;
   outfile << n_throws << " # number of throws in each block" << endl;
   outfile << "# progressive mean of result # progressive uncertainty " << endl;

   double box = 100;
   double d = 1; // don't change this
   double l = 0.2;


   double outer_sum = 0;
   double outer_squares_sum = 0;
   for(int i=0; i<n_blocks; i++){

      double count = 0;
      for(int j=0; j<n_throws; j++){
         double x1 = rnd.Rannyu()*box - box/2.;
         //double y1 = rnd.Rannyu()*box - box/2.; // not needed

         double dx, dy, dist;
         do{
         dx = (rnd.Rannyu()*2* - 1)*l;
         dy = (rnd.Rannyu()*2* - 1)*l;
         dist = sqrt(  pow(dx,2) + pow(dy,2)  );
         }while( dist > l );
         dx = dx/dist * l;
         
         if( floor(x1 + dx) != floor(x1) ){ // needs l = 1 to work
            count++;
         }
      }
      double obs_pi = (2.* l * n_throws)/(count * d);

      
      outer_sum += obs_pi;
      outer_squares_sum += pow(obs_pi , 2);

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
