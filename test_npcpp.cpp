#include <iostream>
#include <cmath>
#include "npcpp.h"

int main()
{
  int nu = 0;
  int c = 0;
  int z = 1;
  int a = 1;
  double E = 1;
  double u = 1e-40 * std::pow(1e13 * 5.067728853, 2);
  int T = 0;

  setreac(nu, c, z, a);
  setsigu(u);
  setsigt(T);
  double sig = signa(E);
  std::cout << sig << std::endl;

    
}
