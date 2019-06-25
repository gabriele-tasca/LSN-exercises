#include <fstream>
using namespace std;
int main()
{
	double start = 0.1;
	double end = 0.48;
	double step = 0.02;
	ofstream f;
	f.open("T.dat",std::ios::trunc);
	for (double t=start; t<=end; t+=step)
		{
			f << t << std::endl;
		}
}
