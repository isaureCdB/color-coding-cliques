/*#include "omp.h"*/

void get_msd(
  int nc1, int nc2, int natom,
  double *c1, double *c2,
  double *msd
) 
{
  /*#pragma omp parallel for schedule(static,4)*/
  for (int i = 0; i < nc1; i++) {
    double *cc1 = c1 + i * natom * 3;
    double *cmsd = msd + i * nc2;
    for (int j = 0; j < nc2; j++) {
      double *cc2 = c2 + j * natom * 3;
      double d = 0;
      for (int k = 0; k < natom; k++) {
        double *a = cc1 + k * 3;
        double *b = cc2 + k * 3;
        for (int l = 0; l < 3; l++) {
          double dif = a[l] - b[l];
          d += dif * dif;
        }
      }
      cmsd[j] = d;
    }
  }  
}