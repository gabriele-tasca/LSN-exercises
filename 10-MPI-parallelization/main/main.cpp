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

#include "mpi.h"

//BASE PARAMETERS
#define NCITIES 100 // up here because they appear in some function definitions
#define NDIM 2

using namespace std;

int main(int argc, char *argv[])
{

   MPI::Init(argc, argv);
   int rank = MPI::COMM_WORLD.Get_rank();
   int size = MPI::COMM_WORLD.Get_size();

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open())
   {
      Primes.seekg(rank * 2); //customized
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
            for (int k = 0; k < 4; k++)
               seed[k] += 729 * rank; //customized
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

   //CITY POSITIONS
   vector<vector<double>> positions(NCITIES + 1);
   for (int jcit = 0; jcit < NCITIES; jcit++)
   {
      positions[jcit] = {0, 0};
   }

   //READ CITY CONFIG FROM FILE
   ifstream posfile("city-pos.dat");
   for (int jcit = 0; jcit < NCITIES; jcit++)
   {
      for (int kdim = 0; kdim < NDIM; kdim++)
      {
         double temp;
         posfile >> temp;
         positions[jcit][kdim] = temp;
      }
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
         vec = random_swap(vec);
      }
      return vec;
   };

   auto full_shift = [&pbc](vector<int> vec, int n) {
      auto copy = vec;
      for (int z = 0; z < NCITIES; z++)
      {
         copy[z] = vec[pbc(z + n)];
      }
      return copy;
   };

   auto invert = [](vector<int> vec) {
      for (int z = 0; z < int(NCITIES / 2); z++)
      {
         int temp = vec[z];
         vec[z] = vec[NCITIES - z - 1];
         vec[NCITIES - z - 1] = temp;
      }
      return vec;
   };

   auto shift = [&pbc](vector<int> vec, int first, int chunk, int shift) {
      auto copy = vec;
      for (int z = first; z < first + chunk; z++)
      {
         auto temp = vec[z];
         vec[z] = vec[pbc(z + shift)];
         vec[pbc(z + shift)] = temp;
      }
      return vec;
   };

   //SETUP
   int nsteps = 200000000;
   int betablock = 100000;
   int maxblock = 50 * betablock;
   int outputskip = betablock;
   int mpiskip = betablock * 20;

   double Tzero = 500;

   double bbetazero = 1. / Tzero;
   double bbetafactor = 1.5;
   double bbetaadd = 0.2;

   vector<int> path(NCITIES);
   for (int jcit = 0; jcit < NCITIES; jcit++)
   {
      path[jcit] = jcit;
   }
   auto best = path;
   double currentlen = len(path);
   double bestlen = currentlen;
   double mpibestlen = -99;
   int mpibestrk = 0;
   double bbeta = bbetazero;

   int accepted = 0;
   int attempted = 0;

   ofstream lenhis("len-his.dat", ios::trunc);

   //METROPOLIS MOVE FUNCTION
   auto metropolis = [&]() {
      auto candid = random_swap(path);

      // auto candid = path;
      // std::random_shuffle(candid.begin(), candid.end());

      double newlen = len(candid);
      if (newlen < bestlen)
      {
         bestlen = newlen;
         best = candid;
      }
      double p = min(1., exp(-bbeta * (newlen - currentlen)));
      if (rnd.Rannyu() < p)
      {
         path = candid;
         accepted++;
         currentlen = len(path);
      }
   };

   //ACTION LOOP
   for (int i = 0; i < nsteps; i++)
   {
      attempted++;
      metropolis();

      if (accepted == betablock || attempted > maxblock)
      {
         // cout << "concluded block with T=" << (1. / bbeta) << " on "<< rank <<" , with acceptance rate = " << double(accepted) / double(attempted) << endl;
         accepted = 0;
         attempted = 0;
         // if(bbeta > bbetaadd){
         //    bbeta += bbetaadd;
         // }else{
         bbeta = bbeta * bbetafactor;
         // }

         // cout << "      current length on " << rank << ": " << currentlen << "     best: " << bestlen << endl;
      }

      if (i % outputskip == 0)
      { // do some checks and write a path, so that the program can be interrupted before the end of the cycle and still give a result
         // NOTE that if you CTRL+C the program you can still interrupt it in the middle of
         // the writing operation and get an incomplete path, if you're unlucky enough

         if (rank != 0)
         {
            MPI::COMM_WORLD.Isend(&bestlen, 1, MPI_DOUBLE, 0, 1);
         }
         if (rank == 0)
         {
            double *box = new double[size];
            box[0] = bestlen;
            for (int rk = 1; rk < size; rk++)
            {
               MPI::COMM_WORLD.Recv(&box[rk], 1, MPI_DOUBLE, rk, 1);
            }

            mpibestlen = 999;
            for (int j = 0; j < size; j++)
            {
               if (box[j] < mpibestlen)
               {
                  mpibestlen = box[j];
                  mpibestrk = j;
               }
            }
            cout << i << "- overall best: " << mpibestlen << " from process " << mpibestrk << "\n";
            lenhis << i << " " << mpibestlen << " " << bbeta << endl;
         }

         //the best process writes the final path on file
         if (rank == mpibestrk)
         {
            //OUTPUT FILES
            ofstream finalpath("final-path.dat", ios::trunc);
            // ofstream finalorder("final-order.dat",ios::trunc);
            for (int kc = 0; kc < NCITIES; kc++)
            {
               finalpath << positions[best[kc]][0] << " " << positions[best[kc]][1] << endl;
               // finalorder  << path[i] << " " << path[i] << endl;
            }
            finalpath.close();
            // finalorder.close();
         }
      }
   }

   lenhis.close();

   // double honkhonk = rnd.Rannyu();
   // // cout << "hOnK " << honkhonk << " which is a double of size " << sizeof(honkhonk) << endl;
   // double *bucket = new double[size];
   // for (int k = 0; k < size; k++)
   // {
   //    bucket[k] = 27;
   // }
   // bucket[0] = 44;

   // int buffer;
   // buffer = rank;

   // // MPI_Gather(&honkhonk, 1, MPI_DOUBLE, bucket, 1, MPI_DOUBLE, 0, MPI::COMM_WORLD);

   // if (rank == 1)
   // {
   //    cout << "FIRST, " << 1 << " has " << buffer << " " << endl;
   // }

   // if (rank == 0)
   // {
   //    MPI::COMM_WORLD.Isend(&buffer, 1, MPI_INTEGER, 1, 1);
   // }
   // if (rank == 1)
   // {
   //    MPI::COMM_WORLD.Recv(&buffer, 1, MPI_INTEGER, 0, 1);
   // }

   // if (rank == 1)
   // {
   //    cout << "THEN, " << 1 << " has " << buffer << " " << endl;
   // }

   // if (rank == 1)
   // {
   //    for (int k = 0; k < size; k++)
   //    {
   //       cout<<"BUT AFTER, " << k << " gave " << bucket[k] << " " << endl;
   //    }
   // }

   rnd.SaveSeed();
   MPI::Finalize();

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
