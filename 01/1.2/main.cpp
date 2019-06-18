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
   int N = 10000;         //number of blocks
   int L = int(M/N);    //throws in each block


   double inner_stdev = 0;
   double variances_sum = 0;
   double variances_squared_sum = 0;


   ofstream out;
   out.open("results.out");
   out << N << endl;
   out << L << endl;
   out << M << endl;
   


   for(int i=0; i<N; i++){

      inner_stdev = 0;
      for(int j=0; j<L; j++){
         double temp = rnd.Rannyu();
         inner_stdev += (temp - 0.5);
      }

      double inner_variance =  pow(inner_stdev,2)/float(L);

      variances_sum += inner_variance;
      variances_squared_sum += pow(inner_variance,2);

      double average_so_far = variances_sum/(i+1);
      double average_of_squares_so_far = variances_squared_sum/(i+1);
      double uncertainty_so_far = sqrt(  (1/(float(i))) * ( average_of_squares_so_far -  pow(average_so_far,2)) );

      out << i+1 << " " << average_so_far << " " << uncertainty_so_far << endl;
   }




   rnd.SaveSeed();
   return 0;
}
