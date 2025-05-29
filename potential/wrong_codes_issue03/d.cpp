#include "d.h"

int main(int argc, char **argv)
{
  if(argc<5) on_syntax_error("d.cpp",  // all distances in units of a
    "d.out b s d N [ymax=2*b dy=b/200  err=1e-12 nmax=2000]");
  bb=atof(argv[1]);
  ss=atof(argv[2]);
  dd=atof(argv[3]);
  N =atoi(argv[4]);
  numero ymax=5*bb;
  numero dy=bb/200;
  if(argc>5) ymax=atof(argv[5]);
  if(argc>6) dy  =atof(argv[6]);
  if(argc>7) err =atof(argv[7]);
  if(argc>8) nmax=atoi(argv[8]);

// Initialize geometry and solve

  set_MM(ss,dd,N);

  // print out (y-b/2)/a and V/V0 for y in the y/a range from 1/2 to ymax
  numero y=bb/2+ds; // start slightly above the top edge of the wire
  do {
    printf("%g %g\n", y-bb/2, potential(0,y)); // is the height relative to the top edge of the wire, for easier interpretation.
    y+=dy;
  } while(y<ymax);

  return 0;
}
