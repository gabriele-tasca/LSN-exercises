/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h> // srand, rand: to generate random number
#include <iostream> // cin, cout: Standard Input/Output Streams Library
#include <fstream>  // Stream class to both read and write from/to files.
#include <cmath>    // rint, pow
#include "random.h"
#include "MolDyn_NVE.h"

#include <vector>

using namespace std;

int main()
{

  Input(); //Inizialization

  int nconf = 1;
  for (int istep = 1; istep <= nstep; ++istep)
  {
    Move(); //Move particles with Verlet algorithm
    if (istep % iprint == 0)
      cout << "Number of time-steps: " << istep << endl;
    if (istep % 10 == 0)
    {
      Measure(); //Properties measurement
      //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
      nconf += 1;
    }

    if (istep % (nstep / nblocks) == 0)
    {
      BlockMeasure(istep / (nstep / nblocks));
    }
  }
  ConfFinal(); //Write final configuration to restart

  return 0;
}

void Input(void)
{

  ifstream ReadInput, ReadConf, ReadConfLast;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl
       << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl
       << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //initialize rand

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

  /////end rand

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart / rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol, 1.0 / 3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblocks;
  ReadInput >> iprint;

  ReadInput >> restart;
  ReadInput >> rescale;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl
       << endl;
  ReadInput.close();

  //Prepare array for measurements
  n_props = 5; //Number of observables
  iv = 0;      //Potential energy
  ik = 1;      //Kinetic energy
  ie = 2;      //Total energy
  it = 3;      //Temperature
  press = 4;

  nbins = 100;
  bin_size = (box / 2.0) / (double)nbins;

  names = {"epot", "ekin", "etot", "temp", "pres"};
  stima_props = vector<double>(n_props, 0);
  store_props = vector<double>(n_props, 0);
  outer_sum = vector<double>(n_props, 0);
  outer_squares_sum = vector<double>(n_props, 0);

  walker_gr = vector<double>(nbins, 0);
  stima_gr = vector<double>(nbins, 0);
  store_gr = vector<double>(nbins, 0);
  outer_sum_gr = vector<double>(nbins, 0);
  outer_squares_sum_gr = vector<double>(nbins, 0);

  //Read initial configuration

  if (restart == 1) //RESTART == TRUE
  {
    cout << "RESTARTING..." << endl;
    ReadConf.open("config.final"); //read config
    ReadConfLast.open("old.final");
    cout << "Read restart configuration from file config.final and old.final" << endl
         << endl;
    for (int i = 0; i < npart; ++i)
    {
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

    for (int i = 0; i < npart; ++i) // read old config
    {
      ReadConfLast >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConfLast.close();

    if (rescale == 1) //if RESCALE == TRUE, try to rescale velocities (old x's) to match the T in input
    {
      cout << "RESCALING..." << endl;
      double v2sum = 0;
      for (int i = 0; i < npart; ++i)
      {
        vx[i] = (x[i] - xold[i]) / delta;
        vy[i] = (y[i] - yold[i]) / delta;
        vz[i] = (z[i] - zold[i]) / delta;
        v2sum += (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2));
      }
      // vmodsum = vmodsum / npart;
      // double tempnow = 3. * vmodsum;
      // double tempfs = (temp / tempnow);
      double oldt;
      double fs;
      // to make it so the program actually works, we read the last averaged temperature from ave_temp
      // rather than using the last velocities
      ifstream tempz;
      tempz.open("tlast.dat");
      if (tempz.is_open())
      {

        tempz >> oldt;
        cout << "OLDT: " << oldt << '\n'; // Display it

        double tfac = temp / oldt;
        fs = sqrt(tfac);
        tempz.close();
      }
      else
      {
        cout << "no ave_temp file from last time \n Using the old method, don't expect this to work though" << endl;
        fs = sqrt(3 * temp / v2sum) * 3; // fs = velocity scale factor
      }

      for (int i = 0; i < npart; ++i)
      {
        vx[i] = vx[i] * fs;
        vy[i] = vy[i] * fs;
        vz[i] = vz[i] * fs;

        xold[i] = x[i] - vx[i] * delta;
        yold[i] = y[i] - vy[i] * delta;
        zold[i] = z[i] - vz[i] * delta;
      }
    }
  }
  else //RESTART == FALSE
  {
    cout << "READING CONFIG.0 AND MAKING UP VELOCITIES..." << endl;
    ReadConf.open("config.0"); //read config

    for (int i = 0; i < npart; ++i) //read config
    {
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

    //generate velocities to match the T in input
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < npart; ++i)
    {
      vx[i] = rnd.Rannyu() - 0.5;
      vy[i] = rnd.Rannyu() - 0.5;
      vz[i] = rnd.Rannyu() - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim = 0; idim < 3; ++idim)
      sumv[idim] /= (double)npart;

    double sumv2 = 0.0, fs;
    for (int i = 0; i < npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2) * 3; // fs = velocity scale factor
    for (int i = 0; i < npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = x[i] - vx[i] * delta;
      yold[i] = y[i] - vy[i] * delta;
      zold[i] = z[i] - vz[i] * delta;
    }
  }

  return;
}

void Move(void)
{ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for (int i = 0; i < npart; ++i)
  { //Force acting on particle i
    fx[i] = Force(i, 0);
    fy[i] = Force(i, 1);
    fz[i] = Force(i, 2);
  }

  for (int i = 0; i < npart; ++i)
  { //Verlet integration scheme

    xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
    ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
    znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

    // there is currently no need to calculate velocities at each step
    // since the Verlet algorithm doesn't use them ...
    vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
    vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
    vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir)
{ //Compute forces as -Grad_ip V(r)
  double f = 0.0;
  double dvec[3], dr;

  for (int i = 0; i < npart; ++i)
  {
    if (i != ip)
    {
      dvec[0] = Pbc(x[ip] - x[i]); // distance ip-i in pbc
      dvec[1] = Pbc(y[ip] - y[i]);
      dvec[2] = Pbc(z[ip] - z[i]);

      dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
      dr = sqrt(dr);

      if (dr < rcut)
      {
        f += dvec[idir] * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure()
{ //Properties measurement
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  //ofstream Epot, Ekin, Etot, Temp;

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

  for (int k = 0; k < nbins; ++k)
  {
    walker_gr[k] = 0.0;
  }

  //cycle over pairs of particles
  for (int i = 0; i < npart - 1; ++i)
  {
    for (int j = i + 1; j < npart; ++j)
    {

      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx * dx + dy * dy + dz * dz;
      dr = sqrt(dr);

      if (dr < rcut)
      {
        vij = 4.0 / pow(dr, 12) - 4.0 / pow(dr, 6);
        //Potential energy
        v += vij;

        pij = (48.0 / pow(dr, 12) - 24.0 / pow(dr, 6));
        //pressure (virial)
        p += pij;
      }

      //update of the histogram of g(r)

      for (int jbin = 0; jbin < nbins; ++jbin)
      {
        if (dr > jbin * bin_size && dr < (jbin + 1) * bin_size)
        {
          walker_gr[jbin] += 2;
          break;
        }
      }
    }
  }
  ///

  //Kinetic energy
  for (int i = 0; i < npart; ++i)
    t += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

  stima_props[iv] = v / (double)npart;                                               //Potential energy
  stima_props[ik] = t / (double)npart;                                               //Kinetic energy
  stima_props[it] = (2.0 / 3.0) * t / (double)npart;                                 //Temperature
  stima_props[ie] = (t + v) / (double)npart;                                         //Total energy
  stima_props[press] = stima_props[it] * rho + (1 / (3. * vol)) * p / (double)npart; //Total energy  //here press

  // vector<ofstream> outf(n_props);// dont think anyone actually needs all the istantaneous values
  // for (int k = 0; k < n_props; k++)
  // {
  //   outf[k].open("output_" + names[k] + ".dat", ios::app);
  // }
  // for (int k = 0; k < n_props; k++)
  // {
  //   outf[k] << stima_props[k] <<  endl;
  //   outf[k].close();
  // }

  for (int k = 0; k < n_props; k++)
  {
    store_props[k] += stima_props[k] / (nstep / nblocks);
  }
  
  for (int kap = 0; kap < nbins; kap++)
  {
    store_gr[kap] += 10*walker_gr[kap] / (nstep / nblocks);
  }

  return;
}

void BlockMeasure(int iblock)
{ //block average properties
  // ofstream ave_Epot, ave_Ekin, ave_Etot, ave_Temp;
  vector<ofstream> outavef(n_props);

  for (int k = 0; k < n_props; k++)
  {
    outavef[k].open("ave_" + names[k] + ".dat", ios::app);

    outer_sum[k] += store_props[k];
    outer_squares_sum[k] += pow(store_props[k], 2);
    store_props[k] = 0;

    double average_so_far = outer_sum[k] / (iblock + 1);
    double average_of_squares_so_far = outer_squares_sum[k] / (iblock + 1);
    double uncertainty_so_far = sqrt((1 / (float(iblock))) * (average_of_squares_so_far - pow(average_so_far, 2)));

    outavef[k] << average_so_far << " " << uncertainty_so_far << endl;
    outavef[k].close();

    if (iblock == nblocks && k == it)
    {
      //final T:
      ofstream finaltf;
      finaltf.open("tlast.dat");
      finaltf << average_so_far << endl;
      cout <<"final temp: "<< average_so_far<<endl;
      finaltf.close();
    }
  }

  //same but for g(r)
  ofstream gr_history_f;
  gr_history_f.open("ave_gr.dat", ios::app);
  for (int k = 0; k < nbins; k++)
  {
    double _r = (k) * bin_size;
    double dvr = 4. * M_PI / 3. * (   pow(_r + bin_size, 3) - pow(_r, 3)  );
    store_gr[k] = store_gr[k] / (rho * npart * dvr );

    outer_sum_gr[k] += store_gr[k];
    outer_squares_sum_gr[k] += pow(store_gr[k], 2);
    store_gr[k] = 0;
    double average_so_far = outer_sum_gr[k] / (iblock + 1);
    double average_of_squares_so_far = outer_squares_sum_gr[k] / (iblock + 1);
    double uncertainty_so_far = sqrt((1 / (float(iblock))) * (average_of_squares_so_far - pow(average_so_far, 2)));
   
    gr_history_f << iblock <<" "<< _r <<" "<< average_so_far << " " << uncertainty_so_far << endl;
  }
  gr_history_f.close();



  if (iblock == nblocks)
    {
      //final g(r):
      ofstream gr_final_f;
      gr_final_f.open("final_gr.dat", ios::trunc);
      for (int k = 0; k < nbins; k++)
      {
        double _r = (k) * bin_size;

        double average_so_far = outer_sum_gr[k] / (iblock + 1);
        double average_of_squares_so_far = outer_squares_sum_gr[k] / (iblock + 1);
        double uncertainty_so_far = sqrt((1 / (float(iblock))) * (average_of_squares_so_far - pow(average_so_far, 2)));
   
        gr_final_f << iblock <<" "<< _r <<" "<< average_so_far << " " << uncertainty_so_far << endl;
      }
      gr_final_f.close();
    }

  return;
}

void ConfFinal(void)
{ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl
       << endl;
  WriteConf.open("config.final");

  for (int i = 0; i < npart; ++i)
  {
    WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
  }
  WriteConf.close();

  //second-last step conf

  ofstream OldWriteConf;

  cout << "Print final configuration to file old.final " << endl
       << endl;
  OldWriteConf.open("old.final");

  for (int i = 0; i < npart; ++i)
  {
    OldWriteConf << xold[i] / box << "   " << yold[i] / box << "   " << zold[i] / box << endl;
  }
  OldWriteConf.close();
  return;
}

void ConfXYZ(int nconf)
{ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i = 0; i < npart; ++i)
  {
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)
{ //Algorithm for periodic boundary conditions with side L=box
  return r - box * rint(r / box);
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
