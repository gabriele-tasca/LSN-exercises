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

   int times = 100;
   int bins = 100;
   int throws_each_time = 10000;

   ofstream out;
   out.open("results.out");
   out << times << " #times" << endl;
   out << bins << " #bins" << endl;
   out << throws_each_time << " #throws each time" << endl;
   out << "# number # chi squared # progressive mean of chi squared # progressive uncertainty of chi squared" << endl;


   vector<int> bin_counter(bins,0);

   double outer_sum = 0;
   double outer_squares_sum = 0;

   for(int i=0; i<times; i++){
      for(int i=0; i<bins; i++){
         bin_counter[i] = 0;
      }
      for(int j=0; j<throws_each_time; j++)
      {
         double x = rnd.Rannyu();
         double scaled_x = x*bins;
         int bin_number = floor(scaled_x);
         bin_counter[bin_number] +=1;
      }
      
      double chisquared = 0;
      double temp_exv = throws_each_time/bins;
      for(int k=0; k<bins; k++){
         chisquared += (  pow(bin_counter[k] - temp_exv, 2)  )  /(temp_exv);
      }

      out << i+1 << " " << chisquared << " ";
      outer_sum += chisquared;
      outer_squares_sum += pow(chisquared,2);


      double average_so_far = outer_sum/(i+1);
      double average_of_squares_so_far = outer_squares_sum/(i+1);
      double uncertainty_so_far = sqrt(  (1/(float(i))) * ( average_of_squares_so_far -  pow(average_so_far,2)) );

      out << average_so_far << " " << uncertainty_so_far << endl;


   }

   rnd.SaveSeed();
   return 0;
}