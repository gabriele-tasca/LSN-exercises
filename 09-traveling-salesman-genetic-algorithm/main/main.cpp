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

//PARAMETERS
#define NCITIES 30
#define NDIM 2
#define NBATCH 900
#define NSTEPS 1000000
#define OUTPUTSKIP 1

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

   //INTEGER PBC FUNCTION
   auto pbc = [](int i) {
      if (i >= NCITIES)
      {
         return i - NCITIES;
      }
      if (i < 0)
      {
         return NCITIES + i;
      }
      else
         return i;
   };

   //SETUP: POSITIONS
   vector<vector<double>> positions(NCITIES + 1);
   for (int jcit = 0; jcit < NCITIES; jcit++)
   {
      positions[jcit] = {0, 0};
   }

   //BASE PATH
   vector<int> path1(NCITIES);
   for (int jcit = 0; jcit < NCITIES; jcit++)
   {
      path1[jcit] = jcit;
   }

   //READ CITY CONFIG FROM FILE
   ifstream posfile("city-pos.dat");
   for (int jcit = 0; jcit < NCITIES; jcit++)
   {
      for (int kdim = 0; kdim < NDIM; kdim++)
      {
         double temp;
         posfile >> temp;
         // std::cout << jcit << " " << temp << " ";
         positions[jcit][kdim] = temp;
      }
      // std::cout << endl;
   }
   posfile.close();

   //DEFINITION OF LENGTH (COST) FUNCTION
   auto len = [&](vector<int> vec) {
      double tot = 0;
      double inner = 0;
      for (int i = 0; i < NCITIES - 1; i++)
      {
         for (int j = 0; j < NDIM; j++)
         {
            inner += pow(positions[vec[i]][j] - positions[vec[i + 1]][j], 2);
         }
         tot += sqrt(inner);
         inner = 0;
      }
      for (int j = 0; j < NDIM; j++)
      {
         inner += pow(positions[vec[NCITIES - 1]][j] - positions[vec[0]][j], 2);
      }
      tot += sqrt(inner);
      return tot;
   };

   //CHECK FUNCTION
   auto is_good = [&](vector<int> vec) {
      auto compare = path1;
      for (int z = 0; z < NCITIES; z++)
      {
         compare[vec[z]] = -99;
      }
      bool allgood = true;
      for (int z = 0; z < NCITIES; z++)
      {
         if (compare[z] != -99)
         {
            allgood = false;
         }
      }
      return allgood;
   };

   //MUTATION FUNCTIONS  // none of which modify the vectors. mutated vector are returned. the [&] is just for aesthetics (and pbc+RNG functions)
   auto swap = [&](vector<int> vec, int index1, int index2) {
      int temp = vec[index1];
      vec[index1] = vec[index2];
      vec[index2] = temp;
      return vec;
   };

   auto random_swap = [&](vector<int> vec) {
      int index1 = floor(rnd.Rannyu() * (NCITIES - 1));
      int index2 = floor(rnd.Rannyu() * (NCITIES - 1));
      return swap(vec, index1, index2);
   };

   auto shuffle = [&](vector<int> vec, int n = 2 * NCITIES) {
      for (int i = 0; i < n; i++)
      {
         random_swap(vec);
      }
      return vec;
   };

   auto full_shift = [&](vector<int> vec, int n) {
      auto copy = vec;
      for (int z = 0; z < NCITIES; z++)
      {
         copy[z] = vec[pbc(z + n)];
      }
      return copy;
   };

   auto invert = [&](vector<int> vec) {
      for (int z = 0; z < int(NCITIES / 2); z++)
      {
         int temp = vec[z];
         vec[z] = vec[NCITIES - z - 1];
         vec[NCITIES - z - 1] = temp;
      }
      return vec;
   };

   auto shift = [&](vector<int> vec, int first, int chunk, int shift) {
      auto copy = vec;
      for (int z = first; z < first + chunk; z++)
      {
         auto temp = vec[z];
         vec[z] = vec[pbc(z + shift)];
         vec[pbc(z + shift)] = temp;
      }
      return vec;
   };

   //CROSSOVER FUNCTION   //note that this doesn't modify the vector but returns copies of them, in vector
   auto cross = [&](vector<int> vec1, vector<int> vec2) {
      const int cut = int(rnd.Rannyu() * (NCITIES - 1));
      //const int cut = int(NCITIES/2);
      auto dummy1 = vec1;
      auto dummy2 = vec2;

      int count = cut;

      for (int i = 0; i < NCITIES; i++)
      {
         //loop in vec2, to add elements in the order in which they appear on vec2
         int now = vec2[i];

         //loop through the tail of the archived vec1 (dummy) to see if "now" needs to be added
         for (int j = cut; j < NCITIES; j++)
         {
            if (now == dummy1[j])
            {
               vec1[count] = now;
               count++;
               break;
            }
         }
         if (count == NCITIES)
            break;
      }
      //  only one son is returned. This creates another son, but whatever. We don't need him.
      //
      //
      // count = cut;
      // for(int i=0; i<NCITIES; i++){
      //    //loop in archived vec1, to add elements in the order in which they appear on vec2
      //    int now = dummy1[i];

      //    //loop through the tail of the archived vec1 (dummy) to see if "now" needs to be added
      //    for (int j = cut; j < NCITIES; j++)
      //    {
      //       if( now == dummy2[j]){
      //          vec2[count] = now;
      //          count++;
      //          break;
      //       }
      //    }
      //    if(count == NCITIES) break;

      // }
      // vector<vector<int>> a = {vec1, vec2};
      return vec1;
   };

   //GENERATE STARTING ELEMENTS
   vector<vector<int>> batch(NBATCH, path1);
   for (int kappa = 0; kappa < NBATCH; kappa++)
   {
      batch[kappa] = shuffle(batch[kappa]);
      // std::random_shuffle(batch[kappa].begin(),batch[kappa].end());
   };

   auto qscompare = [&](const vector<int> &a, const vector<int> &b) {
      if (len(a) < len(b))
      {
         return true;
      }
      else
         return false;
   };

   //FITNESS ORDERING FUNCTION
   auto order = [&]() {
      std::sort(batch.begin(), batch.end(), qscompare);
   };

   //SURVIVAL OF THE FITTEST FUNCTION
   //many versions
   auto live1 = [&]() {
      //random search
      for (int j = 1; j < NBATCH; j++)
      {
         //std::random_shuffle(batch[j].begin(),batch[j].end());
         batch[j] = shuffle(batch[j]);
      }
   };
   auto live2 = [&]() {
      //top ones remains untouched
      //middle has a chance of self-mutating (being replaced by mutated version)
      //bottom is discarded completely and replaced by mutations of top ones (who don't
      //   self-mutate to avoid destructive mutations) and crossovers
      double topfrac = 0.05;
      double midfrac = 0.5;
      int mid = ceil(NBATCH * topfrac);
      int bot = int(NBATCH * midfrac);
      for (int j = mid; j < bot; j++)
      {
         if (rnd.Rannyu() < 0.05)
         {
            batch[j] = random_swap(batch[j]);
         }
         else if (rnd.Rannyu() < 0.05)
         {
            batch[j] = invert(batch[j]);
         }
         else if (rnd.Rannyu() < 0.05)
         {
            batch[j] = shift(batch[j], int(rnd.Rannyu() * NCITIES / 3) + 2, 0, int(rnd.Rannyu() * NCITIES / 3) + 2);
         }
         else if (rnd.Rannyu() < 0.1)
         {
            batch[j] = full_shift(batch[j], int(rnd.Rannyu() * NCITIES / 3) + 1);
         }
         else if (rnd.Rannyu() < 0.05)
         {
            batch[j] = shuffle(batch[j]); //sucks to be you lol
         }
      }
      for (int j = bot; j < NBATCH; j++)
      {
         //same mutations with double rate, and using a top element as a model
         int kgood1 = int(rnd.Rannyu() * NBATCH * topfrac);
         if (rnd.Rannyu() < 0.05)
         {
            batch[j] = random_swap(batch[kgood1]);
         }
         else if (rnd.Rannyu() < 0.05)
         {
            batch[j] = invert(batch[kgood1]);
         }
         else if (rnd.Rannyu() < 0.05)
         {
            batch[j] = shift(batch[kgood1], int(rnd.Rannyu() * NCITIES / 3) + 2, 0, int(rnd.Rannyu() * NCITIES / 3) + 2);
         }
         else if (rnd.Rannyu() < 0.1)
         {
            batch[j] = full_shift(batch[kgood1], int(rnd.Rannyu() * NCITIES / 3) + 1);
         }
         else if (rnd.Rannyu() < 0.05)
         {
            batch[j] = shuffle(batch[j]); //sucks to be you lol
         }
         else if (rnd.Rannyu() < 0.4)
         {
            //crossover with middle guy
            int kgood2 = NBATCH * topfrac + int(rnd.Rannyu() * NBATCH * midfrac);
            batch[j] = cross(batch[kgood1], batch[kgood2]);
         }
         else
         {
            //crossover with TOP guy
            int kgood2 = int(rnd.Rannyu() * NBATCH * topfrac);
            batch[j] = cross(batch[kgood1], batch[kgood2]);
         }
      }
   };

   //ACTION

   ofstream lenhis("len-his.dat", ios::trunc);
   ofstream avlenhis("av-len-his.dat", ios::trunc);

   for (int i = 0; i < NSTEPS; i++)
   {
      live2();
      order();
      if (i % OUTPUTSKIP == 0)
      { // do some checks and write a path, so that the program can be interrupted before the end of the cycle and still give a result
         // NOTE that if you CTRL+C the program you can still interrupt it in the middle of
         // the writing operation and get an incomplete path, if you're unlucky enough.

         if (!is_good(batch[0]))
         {
            cout << "something bad happened: the best path does not pass through each city once" << endl;
            abort();
         }

         double lanlan = len(batch[0]);
         cout << i << "-    " << lanlan << endl;
         lenhis << i << " " << lanlan << endl;

         double avlanlan = 0;
         for (int k = 0; k < int(NBATCH / 2); k++)
         {
            avlanlan += len(batch[0]);
         }
         avlanlan = avlanlan / int(NBATCH / 2);
         avlenhis << i << " " << avlanlan << endl;

         //OUTPUT FILES
         ofstream finalpath("final-path.dat", ios::trunc);
         // ofstream finalorder("final-order.dat",ios::trunc);
         for (int i = 0; i < NCITIES; i++)
         {
            finalpath << positions[batch[0][i]][0] << " " << positions[batch[0][i]][1] << endl;
            // finalorder  << batch[0][i] << " " << batch[0][i] << endl;
         }
         finalpath.close();
         // finalorder.close();
      }
   }

   lenhis.close();
   avlenhis.close();

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
