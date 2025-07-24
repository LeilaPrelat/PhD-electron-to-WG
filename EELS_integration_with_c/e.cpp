#include "e.h"

int main(int argc, char **argv)
{
  if(argc<12) on_syntax_error("e.cpp",  // all distances in units of w
    "e.out Wp/W d/W h/W N epsilon x0/W x1/W nx z0/W z1/W nz [err=1e-12 nmax=2000]");
  numero Wp       =atof(argv[ 1]);
  numero d        =atof(argv[ 2]);
  numero h        =atof(argv[ 3]);
  int    N        =atoi(argv[ 4]);
  numero epsilon  =atof(argv[ 5]);
  numero x0       =atof(argv[ 6]);
  numero x1       =atof(argv[ 7]);
  numero nx       =atoi(argv[ 8]);
  numero z0       =atof(argv[ 9]);
  numero z1       =atof(argv[10]);
  numero nz       =atoi(argv[11]);
  if(argc>12) err =atof(argv[12]);
  if(argc>13) nmax=atoi(argv[13]);

  set_MM(Wp,d,h,N,epsilon);

  int i,j;  numero x,z;
  for(i=0; i<nx; i++)
  for(j=0; j<nz; j++) {
    if(i==0) x=x0; else x=x0+i*(x1-x0)/(nx-1);
    if(j==0) z=z0; else z=z0+j*(z1-z0)/(nz-1);
    if(nx>1) printf("%g ",x);
    if(nz>1) printf("%g ",z);
    printf("%g\n",potential(x,z));
  }

  return 0;
}

// examples:
// a.out 0.5 0.2 0.5 200 9 0 0 1 0.501 2 300 > bo.dat  -->  plot
// a.out 0.5 0.2 0.5 200 9 0 2 200 0.501 2 300 > bo.dat ; density bo.dat; open bo.gif

//  g++ -o e.out e.cpp

