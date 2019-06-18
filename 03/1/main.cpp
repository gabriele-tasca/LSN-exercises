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

   double asset_price_0 = 100;
   double strike_asset_price = 100;
   double delivery_time = 1;
   double interest_rate = 0.1;
   double sigma = 0.25;

   double n_blocks = 10000;
   double n_time_steps = 100;
   double dt = delivery_time/n_time_steps;

   ofstream outfile;
   outfile.open("results.out");
   outfile << n_blocks << " # number of blocks" << endl;
   outfile << n_time_steps << " # number of throws in each block" << endl;
   outfile << "# progressive mean of result # progressive uncertainty " << endl;

   /// 4 files for the 4 options, copy-pasted by hand

   ofstream callfile;
   callfile.open("call.out");
   callfile << n_blocks << " # number of blocks" << endl;
   callfile << n_time_steps << " # number of throws in each block" << endl;
   callfile << "# progressive mean of result # progressive uncertainty " << endl;

   ofstream d_callfile;
   d_callfile.open("d_call.out");
   d_callfile << n_blocks << " # number of blocks" << endl;
   d_callfile << n_time_steps << " # number of throws in each block" << endl;
   d_callfile << "# progressive mean of result # progressive uncertainty " << endl;

   ofstream putfile;
   putfile.open("put.out");
   putfile << n_blocks << " # number of blocks" << endl;
   putfile << n_time_steps << " # number of throws in each block" << endl;
   putfile << "# progressive mean of result # progressive uncertainty " << endl;

   ofstream d_putfile;
   d_putfile.open("d_put.out");
   d_putfile << n_blocks << " # number of blocks" << endl;
   d_putfile << n_time_steps << " # number of throws in each block" << endl;
   d_putfile << "# progressive mean of result # progressive uncertainty " << endl;

   double asset_sum = 0;
   double asset_squares_sum = 0;

   double d_asset_sum = 0;
   double d_asset_squares_sum = 0;

   double call_sum = 0;
   double call_squares_sum = 0;

   double d_call_sum = 0;
   double d_call_squares_sum = 0;

   double put_sum = 0;
   double put_squares_sum = 0;

   double d_put_sum = 0;
   double d_put_squares_sum = 0;


   for(int i=0; i<n_blocks; i++){
      ////////////direct
      double asset_price_t = asset_price_0 * exp( (interest_rate - 0.5*pow(sigma,2))*delivery_time + sigma* rnd.Gauss(0,delivery_time) );
      asset_sum += asset_price_t;
      asset_squares_sum += pow(asset_price_t , 2);

      double call_profit = max( 0., asset_price_t - strike_asset_price );
      double call_option_price = call_profit * exp( - delivery_time*interest_rate);
      call_sum += call_option_price;
      call_squares_sum += pow(call_option_price , 2);

      double put_profit = max( 0., -asset_price_t + strike_asset_price );
      double put_option_price = put_profit * exp( - delivery_time*interest_rate);
      put_sum += put_option_price;
      put_squares_sum += pow(put_option_price , 2);

      ////////////discrete
      double d_asset_price_t = asset_price_0;
      for(int j=0; j<n_time_steps; j++){
         d_asset_price_t = d_asset_price_t*exp(  (interest_rate - 0.5*pow(sigma,2))*(dt) + sigma* rnd.Gauss(0,1)*sqrt(dt) );
      }
      d_asset_sum += d_asset_price_t;
      d_asset_squares_sum += pow(d_asset_price_t , 2);

      double d_call_profit = max( 0., d_asset_price_t - strike_asset_price );
      double d_call_option_price = d_call_profit * exp( - delivery_time*interest_rate);
      d_call_sum += d_call_option_price;
      d_call_squares_sum += pow(d_call_option_price , 2);

      double d_put_profit = max( 0., -d_asset_price_t + strike_asset_price );
      double d_put_option_price = d_put_profit * exp( - delivery_time*interest_rate);
      d_put_sum += d_put_option_price;
      d_put_squares_sum += pow(d_put_option_price , 2);



      double call_average_so_far = call_sum/(i+1);
      double call_average_of_squares_so_far = call_squares_sum/(i+1);
      double call_uncertainty_so_far = sqrt(  (1/(float(i))) * ( call_average_of_squares_so_far -  pow(call_average_so_far,2)) );

      callfile << call_average_so_far << " " << call_uncertainty_so_far << endl;


      double d_call_average_so_far = d_call_sum/(i+1);
      double d_call_average_of_squares_so_far = d_call_squares_sum/(i+1);
      double d_call_uncertainty_so_far = sqrt(  (1/(float(i))) * ( d_call_average_of_squares_so_far -  pow(d_call_average_so_far,2)) );

      d_callfile << d_call_average_so_far << " " << d_call_uncertainty_so_far << endl;


      double put_average_so_far = put_sum/(i+1);
      double put_average_of_squares_so_far = put_squares_sum/(i+1);
      double put_uncertainty_so_far = sqrt(  (1/(float(i))) * ( put_average_of_squares_so_far -  pow(put_average_so_far,2)) );

      putfile << put_average_so_far << " " << put_uncertainty_so_far << endl;


      double d_put_average_so_far = d_put_sum/(i+1);
      double d_put_average_of_squares_so_far = d_put_squares_sum/(i+1);
      double d_put_uncertainty_so_far = sqrt(  (1/(float(i))) * ( d_put_average_of_squares_so_far -  pow(d_put_average_so_far,2)) );

      d_putfile << d_put_average_so_far << " " << d_put_uncertainty_so_far << endl;

         // double call_final_average = call_sum/n_blocks;
         // double call_final_average_of_squares = call_squares_sum/n_blocks;
         // double call_uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( call_final_average_of_squares -  pow(call_final_average,2)) );

      

   }
   
   double final_average = asset_sum/n_blocks;
   double final_average_of_squares = asset_squares_sum/n_blocks;
   double uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( final_average_of_squares -  pow(final_average,2)) );

   outfile << final_average << " " << uncertainty_so_far << endl;


   double d_final_average = d_asset_sum/n_blocks;
   double d_final_average_of_squares = d_asset_squares_sum/n_blocks;
   double d_uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( d_final_average_of_squares -  pow(d_final_average,2)) );

   outfile << d_final_average << " " << d_uncertainty_so_far << endl;
   
//////// 4 options for call / put / discrete / direct

   double call_final_average = call_sum/n_blocks;
   double call_final_average_of_squares = call_squares_sum/n_blocks;
   double call_uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( call_final_average_of_squares -  pow(call_final_average,2)) );

   callfile << call_final_average << " " << call_uncertainty_so_far << endl;

   double d_call_final_average = d_call_sum/n_blocks;
   double d_call_final_average_of_squares = d_call_squares_sum/n_blocks;
   double d_call_uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( d_call_final_average_of_squares -  pow(d_call_final_average,2)) );

   outfile << d_call_final_average << " " << d_call_uncertainty_so_far << endl;

   double put_final_average = put_sum/n_blocks;
   double put_final_average_of_squares = put_squares_sum/n_blocks;
   double put_uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( put_final_average_of_squares -  pow(put_final_average,2)) );

   outfile << put_final_average << " " << put_uncertainty_so_far << endl;

   double d_put_final_average = d_put_sum/n_blocks;
   double d_put_final_average_of_squares = d_put_squares_sum/n_blocks;
   double d_put_uncertainty_so_far = sqrt(  (1/(float(n_blocks-1))) * ( d_put_final_average_of_squares -  pow(d_put_final_average,2)) );

   outfile << d_put_final_average << " " << d_put_uncertainty_so_far << endl;
   cout << d_put_final_average << " " << d_put_uncertainty_so_far << endl;


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
