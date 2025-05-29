#include "d.h"

int main(int argc, char **argv)
{
  if(argc<5) on_syntax_error("d.cpp",  // all distances in units of a
    "d.out b s d N [ymax=2*b dy=b/200  err=1e-12 nmax=2000]");
  bb=atof(argv[1]);
  ss=atof(argv[2]);
  dd=atof(argv[3]);
  N =atoi(argv[4]);
  numero ymax=bb/2+dd+70*dd;
//  numero ymin=-bb/2-dd;
  numero ymin=bb/2;
  numero dy=bb/200;
  if(argc>5) ymax=atof(argv[5]);
  if(argc>6) dy  =atof(argv[6]);
  if(argc>7) err =atof(argv[7]);
  if(argc>8) nmax=atoi(argv[8]);

// Initialize geometry and solve

  set_MM(ss,dd,N);

  // Loop over grid and compute potential
  for (numero y = ymin; y <= ymax; y += dy) {
      numero v = potential(0, y);
      printf("%g %g\n", y, v);
  }
  
    return 0;
}
