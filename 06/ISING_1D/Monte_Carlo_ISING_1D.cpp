/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//aaa
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    // cerr << "completed reset" << endl;
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      // cerr << "completed move " << istep << endl;
      Measure();
      // cerr << "completed measure " << istep << endl;
      Accumulate(); //Update block averages
      // cerr << "completed accumulation " << istep << endl;
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- bbeta * H ), bbeta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1=0, p2=0;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  bbeta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  iu2 = 1; //Energy squared
  is = 2; //Spin
  is2 = 3; //Spin squared

  giu = 0;
  gic = 1;
  gim = 2;
  gix = 3;

 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o = 0;
  double p, energy_old, energy_new;

  for(int i=0; i<nspin; i++)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = floor(rnd.Rannyu()*double(nspin));
    // cerr << "o index " << o << endl;

    if(metro==1) //Metropolis
    {
      energy_old = Boltzmann(s[o], o);
      energy_new = Boltzmann(-s[o], o);
      p = min(1., exp(- bbeta* (energy_new - energy_old)) );
      if( rnd.Rannyu() < p){
        s[o] = -s[o];
        accepted++;
      }
      attempted++;
    }
    else //Gibbs sampling
    {
      p = 1/( 1 + exp(-bbeta*(Boltzmann(1, o) - Boltzmann(-1, o)) ) );
      if(rnd.Rannyu() < p){
        s[o] = -1;
      }else{
        s[o] = 1;
      }
    }

  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, u2 = 0.0, sp = 0.0, sp2 = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    sp += s[i];
  }
  sp2 = sp*sp;
  u2 += u*u;

  walker[iu] = u;
  walker[iu2] = u2;
  walker[is] = sp;
  walker[is2] = sp2;  
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
  
  ofstream Ene, Heat, Mag, Chi;
  const int wd=20;
  
  cout << "Block number " << iblk << endl;
  if(metro == 1){
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  }
  Ene.open("output.ene.0",ios::app);
  Heat.open("output.heat.0",ios::app);
  Mag.open("output.mag.0",ios::app);
  Chi.open("output.chi.0",ios::app);
  
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  stima_u2 = blk_av[iu2]/blk_norm/ (double)nspin /double(nspin);
  // stima_u2 = blk_av[iu2]/blk_norm/  (  (double)nspin*(double)nspin );

  stima_c = bbeta*bbeta*( stima_u2 - stima_u*stima_u  )*double(nspin);

  stima_s = blk_av[is]/blk_norm/(double)nspin; 
  stima_s2 = blk_av[is2]/blk_norm/(double)nspin/(double)nspin;
  stima_m = stima_s;
  stima_x = bbeta*stima_s2*double(nspin);

  glob_av[0]  += stima_u;
  glob_av2[0] += stima_u*stima_u;
  err_u=Error(glob_av[0],glob_av2[0],iblk);
  Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[0]/(double)iblk << setw(wd) << err_u << endl;
  Ene.close();

  glob_av[gic]  += stima_c;
  glob_av2[gic] += stima_c*stima_c;
  err_u=Error(glob_av[gic],glob_av2[gic],iblk);
  Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[gic]/(double)iblk << setw(wd) << err_u << endl;
  Heat.close();

  glob_av[gim]  += stima_m;
  glob_av2[gim] += stima_m*stima_m;
  err_u=Error(glob_av[gim],glob_av2[gim],iblk);
  Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[gim]/(double)iblk << setw(wd) << err_u << endl;
  Mag.close();

  glob_av[gix]  += stima_x;
  glob_av2[gix] += stima_x*stima_x;
  err_u=Error(glob_av[gix],glob_av2[gix],iblk);
  Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[gix]/(double)iblk << setw(wd) << err_u << endl;
  Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
