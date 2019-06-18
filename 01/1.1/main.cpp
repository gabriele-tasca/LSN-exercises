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


   int M = 10000000;      //total random points
   int N = 1000;         //number of blocks
   int L = int(M/N);    //throws in each block


   double inner_sum = 0;

   double mean_values_sum = 0;
   double mean_values_squared_sum = 0;


   ofstream out;
   out.open("results.out");
   out << N << endl;
   out << L << endl;
   out << M << endl;
   


   for(int i=0; i<N; i++){

      inner_sum = 0;
      for(int j=0; j<L; j++){
         double temp = rnd.Rannyu();
         inner_sum += temp;
      }

      double mean_value = inner_sum/L;
      mean_values_sum += mean_value;
      mean_values_squared_sum += pow(mean_value,2);

      double average_so_far = mean_values_sum/(i+1);
      double average_of_squares_so_far = mean_values_squared_sum/(i+1);
      double uncertainty_so_far = sqrt(  (1/(float(i))) * ( average_of_squares_so_far -  pow(average_so_far,2)) );

      out << i+1 << " " << average_so_far << " " << uncertainty_so_far << endl;
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
