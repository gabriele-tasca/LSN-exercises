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

   //PARAMETERS

   int ndim = 2;
   int ncities = 100;
   double L = 10;


   auto circle_random = [&](){
      //ndim is overridden and fixed as 2
      double theta = rnd.Rannyu()*2*M_PI;
      vector<double> pos(2);
      pos[0] = L*cos(theta);
      pos[1] = L*sin(theta);
      return pos;
   };

   auto box_random = [&](){
      vector<double> pos(ndim);
      for (int j=0; j<ndim; j++){
         pos[j] = rnd.Rannyu()*L;
      }
      return pos;
   };

   auto randomize = box_random;



   vector<vector<double>> positions(ncities);
   
   ofstream outf("city-pos.dat");

   for(int j=0; j<ncities; j++){
      auto temp = randomize();
      outf << temp[0] << " " << temp[1] << endl;
   }



   outf.close();

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
