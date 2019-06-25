/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <string>
#include <vector>

//parameters, observables
const int m_props=5;
int n_props;
int iv,ik,it,ie,press;
//double stima_pot, stima_kin, stima_etot, stima_temp;
std::vector<std::string> names;
std::vector<double> stima_props;
std::vector<double> store_props;
std::vector<double> outer_sum;
std::vector<double> outer_squares_sum;


double bin_size, nbins;
std::vector<double> walker_gr;
std::vector<double> stima_gr;
std::vector<double> store_gr;
std::vector<double> outer_sum_gr;
std::vector<double> outer_squares_sum_gr;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblocks, iprint, seed;
int restart, rescale; //actually a bool
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void BlockMeasure(int);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
