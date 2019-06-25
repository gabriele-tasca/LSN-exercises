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
#include <string>

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

   int N_blocks = 10000;
   vector<int> L_throws_v({1,2,10,100});

   int len = L_throws_v.size();
   vector<ofstream> unif_outfiles_v(len); 
   vector<ofstream> exp_outfiles_v(len);
   vector<ofstream> lorentz_outfiles_v(len);
   for(int i=0; i<len; i++){  
      unif_outfiles_v[i].open("unif_"+to_string(L_throws_v[i])+".out" );
      unif_outfiles_v[i] << N_blocks << " # number of blocks" << endl;
      unif_outfiles_v[i] << L_throws_v[i] << " # number of throws per block" << endl;
      
      for (int k=0; k<N_blocks; k++){

         double unif_average = 0;
         for (int j=0; j<L_throws_v[i]; j++){
            unif_average +=  rnd.Rannyu() / L_throws_v[i];
         }
         unif_outfiles_v[i] << unif_average << endl;
      }

      exp_outfiles_v[i].open("exp_"+to_string(L_throws_v[i])+".out" ); 
      exp_outfiles_v[i] << N_blocks << " # number of blocks" << endl;
      exp_outfiles_v[i] << L_throws_v[i] << " # number of throws per block" << endl;
      
      for (int k=0; k<N_blocks; k++){

         double exp_average = 0;
         for (int j=0; j<L_throws_v[i]; j++){
            exp_average +=  rnd.Exp(1.) / L_throws_v[i];
         }
         exp_outfiles_v[i] << exp_average << endl;
      }

      lorentz_outfiles_v[i].open("lorentz_"+to_string(L_throws_v[i])+".out" );
      lorentz_outfiles_v[i] << N_blocks << " # number of blocks" << endl;
      lorentz_outfiles_v[i] << L_throws_v[i] << " # number of throws per block" << endl;

      for (int k=0; k<N_blocks; k++){

         double lorentz_average = 0;
         for (int j=0; j<L_throws_v[i]; j++){
            lorentz_average +=  rnd.Lorentz(0., 1.) / L_throws_v[i];
         }
         lorentz_outfiles_v[i] << lorentz_average << endl;
      }

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
