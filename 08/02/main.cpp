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

//SIGN FUNCTION
template <typename T>
int sgn(T val)
{
   return (T(0) < val) - (val < T(0));
}


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

   //BASE FUNCTIONS

   auto psiT = [](double x, double mu, double sigma) {
      return exp(-(pow((x - mu), 2)) / (2 * pow(sigma, 2))) + exp(-(pow((x + mu), 2)) / (2 * pow(sigma, 2)));
   };
   // auto psiT = [](double x, double mu, double sigma) {
   //    if( fabs(x) > 1) {return 0;}
   //    else {return 1;}
   // };

   auto D2_psiT = [](double x, double mu, double sigma) {
      //return (1/(pow(sigma,4)))*exp( -( pow((x - mu),2) )/(2*sigma*sigma) )   *   ((1 +  exp( (2*mu*x)/(sigma*sigma) ) )*(mu*mu -sigma*sigma + x*x )  - 2*(-1+exp( (2*mu*x)/(sigma*sigma) )*mu*x )  );
      //return (1 / pow(sigma, 2)) * exp(-(pow(mu + x, 2)) / (2 * sigma * sigma)) * (exp(-pow(mu + x, 2) / sigma * sigma) * (mu - sigma - x) * (mu + sigma - x) + pow(mu + x, 2) - sigma * sigma);
      //return (-1. / pow(sigma, 2.)) * (exp(-pow(x - mu, 2.) / (2. * pow(sigma, 2.))) + -pow(x + mu, 2.) / (2. * pow(sigma, 2.))) + (1. / pow(sigma, 4)) * (pow(x - mu, 2.) * exp(-pow(x - mu, 2.) / (2. * pow(sigma, 2.))) + pow(x + mu, 2.) * exp(-pow(x + mu, 2.) / (2. * pow(sigma, 2.))));
      double s2 = pow(sigma,2);
      double s4 = pow(sigma,4);
      double xmu2plus = pow(x+ mu, 2);
      double xmu2min = pow(x- mu, 2);
      return -exp(-xmu2min/(2*s2))/s2 - exp(-xmu2plus/(2*s2))/s2 + (  xmu2min*exp(-xmu2min/(2*s2)) )/s4 + (  xmu2plus*exp(-xmu2plus/(2*s2)) )/s4;
   };

   auto V = [&psiT](double x) {
      return (pow(x, 4.) - 2.5 * pow(x, 2.));
   };  


   auto H_psiT = [&D2_psiT, &V, &psiT](double x, double mu, double sigma) {
      return -D2_psiT(x, mu, sigma) + V(x) * psiT(x, mu, sigma);
   };

   auto gauss = [&rnd](double mmu, double msigma) {
      mmu = rnd.Gauss(mmu, msigma);
      return mmu;
   };

   auto unif = [&rnd](double mean, double delt) {
      return rnd.Rannyu() * delt + mean - 0.5 * delt;
   };

   //METROPOLIS STEP PARAMS
   auto metro_T = unif;
   auto metro_f = [&psiT](double x, double mu, double sigma) {
      return pow(psiT(x, mu, sigma) ,2);
   };
   double metro_delta_zero = 1.2;

   //METROPOLIS STEP FUNCTION
   auto metropolis_step = [&metro_T, &metro_f, &rnd](double x, double mu, double sigma, int &a_count, double metro_delta) {
      double xnew = metro_T(x, metro_delta);
      double alfa = min(1., metro_f(xnew, mu, sigma) / metro_f(x, mu, sigma));
      if (rnd.Rannyu() < alfa)
      {
         x = xnew;
         a_count++;
      }
      if (rnd.Rannyu() < 0.5)  // manually enforces symmetry 
      {
         x = -x;
      }
      return x;
   };

   struct datapoint
   {
   public:
      double value;
      double error;
   };

   auto equilibrate = [&metropolis_step](int n_equilib, double x, double mu, double sigma) {
      int a_count = 0;
      for (int j = 0; j < n_equilib; j++)
      {
         x = metropolis_step(x, mu, sigma, a_count, sigma/2.);
      }
      cout << "equilibration phase completed. acceptance %:" << double(a_count) / double(n_equilib + 1) << endl;
      return x;
   };

   auto calc_energy = [&H_psiT, &D2_psiT, &psiT, &metropolis_step](double x, double mu, double sigma, int n_blocks, int n_throws){
      double glob_energy_sum = 0;
      double glob_energy_sum2 = 0;
      double glob_energy_error = 0;

      double current_e = 0;
      for (int jb = 0; jb < n_blocks; jb++)
      {
         int b_count = 0;
         double energy_sum = 0;
         double energy_sum2 = 0;
         for (int jt = 0; jt < n_throws; jt++)
         {
            x = metropolis_step(x, mu, sigma, b_count, sigma/2.);


            double H = H_psiT(x, mu, sigma)   / psiT(x, mu, sigma) ;

            //H = pow(x,3);

            //double H = D2_psiT(x, mu, sigma);
            // cout << H << "\n";
            energy_sum += H;
         }
         // cout << "acc metro int : " << double(b_count)/double(n_throws) << "\n";
         double energy_av = energy_sum / double(n_throws);
         //cerr << energy_av << "\n";

         // double energy_av2 = energy_sum2/double(n_throws);

         glob_energy_sum += energy_av;
         glob_energy_sum2 += energy_av * energy_av;
         glob_energy_error = sqrt( glob_energy_sum2/double(jb+1) - pow(glob_energy_sum/double(jb+1),2 )  )/sqrt(double(jb+1));


         current_e = glob_energy_sum / double(jb+1);
         //outfile << current_e << " " << glob_energy_error << endl;
         energy_sum = 0;
         energy_sum2 = 0;
         //cout << jb <<"th phase over. acceptance %:" << double(b_count)/double(n_throws +1) << endl;
      }

      datapoint energy;
      energy.value = current_e;
      energy.error = glob_energy_error;
      return energy;
   };

   auto sample = [&metropolis_step](double x, double mu, double sigma, int n) {
      int a_count = 0;
      vector<double> smpl(n, 0);
      for (int jt = 0; jt < n; jt++)
      {
         x = metropolis_step(x, mu, sigma, a_count, sigma/2.);
         smpl[jt] = x;
      }
      return smpl;
   };

   ////PARAMETERS

   int n_equilib = 200000;
   int n_throws = 10000;
   int n_blocks = 2*100;

   double x0 = 0;
   double sigmazero = 1.2;
   double muzero = 1.2;

   double last_sigma = 1.3;
   double last_mu = 1.7;

   int SA_n = 130;
   double SAsigmastep = 0.1;
   double SAmustep = 0.1;
   double SA_betazero = 100;
   double SA_betafac = 1.6;

   ofstream outfile;
   outfile.open("results.out", ios::trunc);
   outfile << n_throws << " # number of throws" << endl;

   ////////////////////////////////////////////////////////////////
   //////////////////////START/////////////////////////////////////
   ////////////////////////////////////////////////////////////////

   //EQUILIBRATION
   double x = x0;
   // x = equilibrate(n_equilib, x, muzero, sigmazero);


   //MINIMIZATION LOOP
   double sigma = sigmazero;
   double mu = muzero;
   double bbeta = 2;
   double lasten = calc_energy(x, muzero, sigmazero, n_blocks, n_throws).value;
   int accepted = 0;
   for (int jsa = 0; jsa < SA_n; jsa++)
   {

      double newsigma = rnd.Gauss(sigma, SAsigmastep);
      double newmu = rnd.Gauss(mu, SAmustep);

      double newen = calc_energy(x, newmu, newsigma, n_blocks, n_throws).value;

      double p = min(1., exp(-bbeta * (newen - lasten)));
      if (rnd.Rannyu() < p)
      {
         sigma = newsigma;
         mu = newmu;
         lasten = newen;
         accepted++;
         cout << "mu: " << mu << " | sigma: " << sigma << " | energy: " << newen << " | acceptance so far: " << double(accepted) / double(jsa+1) << endl;
         outfile << "mu: " << mu << " | sigma: " << sigma << " | energy: " << newen << " | acceptance so far: " << double(accepted) / double(jsa+1) << endl;
      }
      else
      {
         //cout << "        failed jump...." << endl;
      }
      if (double(accepted) / double(jsa+1) < 0.3 && jsa > 20)
      {
         SAmustep = SAmustep * 0.5;
         SAsigmastep = SAsigmastep * 0.3;
         cout << "trying to reduce step... acc: " << double(accepted) / double(jsa+1) << endl;
      }
      if (double(accepted) / double(jsa+1) < 0.2 && jsa > 20)  
      {
         break;
      }

      bbeta = bbeta * SA_betafac;
   }

   // cout <<"mu:"<< endl;
   // cin >> mu;   
   // cout <<"sigma:"<< endl;
   // cin >> sigma;
   // cout <<"en : "<< calc_energy(x, mu, sigma, n_blocks, 2*n_throws).value << endl;

   //PRINT FINAL PARAMETERS
   ofstream parfile;
   parfile.open("par.out", ios::trunc);

   parfile << "mu     |     sigma" << endl;
   parfile << mu << " " << sigma << endl;

   //SAMPLING OF FINAL WAVeFUNCTION

   int smplsz = 500000;
   ofstream samplefile;
   samplefile.open("sample.out", ios::trunc);

   auto ttt = sample(x0, mu, sigma, smplsz);
   for (int j = 0; j < smplsz; j++)
   {
      samplefile << ttt[j] << endl;
   }

   ofstream psifile;
   psifile.open("psi.out", ios::trunc);
   double wid = 5;
   for(double x=-wid; x<wid; x+=0.05)
   {
      psifile << x << " " << D2_psiT(x, mu, sigma) << endl;
   }
   psifile.close();
   // return 0;

   //CALCULATION OF FINAL ENERGY PRINTING EXPLICIT HISTORY

   ofstream hisfile;
   hisfile.open("his.out", ios::trunc);

   double glob_energy_sum = 0;
   double glob_energy_sum2 = 0;
   double glob_energy_error = 0;

   double current_e = 0;
   for (int jb = 0; jb < 150; jb++)
   {
      int a_count = 0;
      double energy_sum = 0;
      double energy_sum2 = 0;
      for (int jt = 0; jt < n_throws; jt++)
      {
         x = metropolis_step(x, mu, sigma, a_count, sigma/2.);
         // double H = 1.00000;
         // double H = psiT(x, mu, sigma)* psiT(x,mu,sigma);
         double H = H_psiT(x, mu, sigma) / psiT(x, mu, sigma) ;
         energy_sum += H;
      }
      double energy_av = energy_sum / double(n_throws);
      // double energy_av2 = energy_sum2/double(n_throws);

      glob_energy_sum += energy_av;
      glob_energy_sum2 += energy_av * energy_av;
      //glob_energy_error = sqrt((-glob_energy_sum2 / (double)jb + pow(glob_energy_sum / (double)jb, 2)) / (double)jb);
      glob_energy_error = sqrt( glob_energy_sum2/double(jb+1) - pow(glob_energy_sum/double(jb+1),2 )  )/sqrt(double(jb+1));

      current_e = glob_energy_sum / double(jb+1);
      hisfile << jb << " " << current_e << " " << glob_energy_error << endl;
      // cerr << jb<< " "<< " "<<  glob_energy_sum2/double(jb+1)   << " "<< pow(glob_energy_sum/double(jb+1),2 ) <<"->"<< glob_energy_error <<"\n";

      energy_sum = 0;
      energy_sum2 = 0;
      //cout << jb <<"th phase over. acceptance %:" << double(a_count)/double(n_throws +1) << endl;
   }

   hisfile.close();
   samplefile.close();
   outfile.close();
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
