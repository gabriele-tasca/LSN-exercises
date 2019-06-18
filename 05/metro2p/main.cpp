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
#include <functional>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open())
   {
      Primes >> p1 >> p2;
   }
   else
      cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open())
   {
      while (!input.eof())
      {
         input >> property;
         if (property == "RANDOMSEED")
         {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);
         }
      }
      input.close();
   }
   else
      cerr << "PROBLEM: Unable to open seed.in" << endl;

   ///////////////////////

   int n_blocks = 100;
   int n_throws = 2000;

   vector<double> x0 = {0, 0, 0};
   vector<double> x = x0;
   double sigma = 1.8;

   // auto mod_psi_1s = [](vector<double> x) {
   //    double r = 0;
   //    for (int i = 0; i < 3; i++)
   //    {
   //       r += pow(x[i], 2);
   //    }
   //    r = sqrt(r);
   //    return pow(M_E, -2. * r);
   // };

   auto mod_psi_2p = [](vector<double> x) { //squared modulus!
      double r = 0;
      for (int i = 0; i < 3; i++)
      {
         r += pow(x[i], 2);
      }
      r=sqrt(r);
      //double theta = acos(x[2]/r);
      auto hgh = r*r * exp(-r) * pow(x[2] / r, 2.);
      return hgh;
   };

   auto gauss = [&rnd](vector<double> mu, double sigma) {
      for (int i = 0; i < 3; i++)
      {
         mu[i] = rnd.Gauss(mu[i], sigma);
      }
      return mu;
   };

   ofstream hfile;
   hfile.open("h.out", ios::trunc);

   ofstream smplfile;
   smplfile.open("smpl.out", ios::trunc);

   double inner_sum = 0;
   double outer_sum = 0;
   double outer_sum_2 = 0;
   for (int i = 0; i < n_blocks; i++)
   {
      inner_sum = 0;
      int acc = 0;
      for (int j = 0; j < n_throws; j++)
      {

         vector<double> xnew = gauss(x, sigma);
         double alfa = min(1., mod_psi_2p(xnew) / mod_psi_2p(x));

         if (rnd.Rannyu() < alfa)
         {
            x = xnew;
            acc++;
         }
         smplfile << x[0] << " " << x[1] << " " << x[2] << endl;
         inner_sum += sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
      }
      cout << "accept: " << double(acc) / double(n_throws) << endl;
      inner_sum = inner_sum / n_throws;

      outer_sum += inner_sum;
      outer_sum_2 += pow(inner_sum, 2);

      double average_so_far = outer_sum / (i + 1);
      double average_of_squares_so_far = outer_sum_2 / (i + 1);
      double uncertainty_so_far = sqrt((1 / (float(i))) * (average_of_squares_so_far - pow(average_so_far, 2)));

      hfile << average_so_far << " " << uncertainty_so_far << endl;
   }

   smplfile.close();
   hfile.close();
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