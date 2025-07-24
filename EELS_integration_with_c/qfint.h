#ifndef jga_qfint         // ************************************************
#define jga_qfint 1       // ***  jga/qfint.h    Translated from FORTRAN  ***
                          // ***                                28-vi-90  ***
#include "jga.h"          // ***                                  8-x-97  ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Integration routines by quadratures                               ***
// ***                                                                    ***
// ***                                     /b                             ***
// ***      qfint (a, b, f, err, nmax) =  |     f(x)  dx   +/-  err       ***
// ***                                   /a                               ***
// ***                                                                    ***
// ***      nmax = maximum number of evaluation of f(x)                   ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  --- numero  qfint(a, b, numero  func(numero x), err, mmax);       ***
// ***  --- complex qfint(a, b, complex func(numero x), err, mmax);       ***
// ***                                                                    ***
// ***  qfint_init(nmax) initializes the quadrature, but it is invoked    ***
// ***  automatically when calling qfint(...)                             ***
// ***                                                                    ***
// **************************************************************************


numero *qfint_xv, *qfint_wv;
int qfint_mmax=0;

// --------------------------------------------------------------------------

int qfint_init(int mmax)
{
  if(mmax<=qfint_mmax) return 0;

  numero  pi2=2*pi, theta;
  int  i,j,m;

  if(qfint_mmax>0) {delete [] qfint_xv;  delete [] qfint_wv;}

  qfint_xv=new numero [mmax+2];
  qfint_wv=new numero [mmax+2];

  j=1;  m=4;  qfint_mmax=mmax;

  do {
    for(i=1; i<m; i++)  if(i%2 || m==4) {
      theta = i * (pi2 / m);
      qfint_xv [j]= (theta-sin(theta)) / pi2;
      qfint_wv [j]=  1  - cos(theta);
      j++;
    }
    m=2*m;
  } while (m<=mmax || m<=8);

  return 0;
}

// --------------------------------------------------------------------------

numero qfint(numero a, numero b, numero func(numero x), numero err, int mmax)
{
  numero  esm,sigmam;   int  i,j,m;
  numero  eta,e;   // e is the error

  if(mmax>qfint_mmax)  qfint_init(mmax);

  if(a==b)  return 0;

  m=4;  j=1;  esm=0;
  do {
    sigmam=0;
    for(i=1; i<m; i++)  if(i%2 || m==4) {
      eta=(b-a)*qfint_xv[j]+a;
      sigmam+=qfint_wv[j]*func(eta);
      j++;
    }
    sigmam=2*(b-a)/m*sigmam;
    e     =(esm-sigmam)/2;   if(e<0) e=-e;
    esm   =0.5*(esm+sigmam);
    m=2*m;
  } while ((e>err && m<=mmax)  || m<=8);

  return  esm;
}


// **************************************************************************
// *** complex functions                                                  ***
// **************************************************************************

#ifdef jga_complex

complex qfint(numero a, numero b, complex func(numero x), numero err, int mmax)
{
  complex  esm,sigmam;   int  i,j,m;
  numero   eta, e;   // e is the error

  if(mmax>qfint_mmax)  qfint_init(mmax);

  if(a==b)  return 0;

  m=4;  j=1;  esm=0;
  do {
    sigmam=0;
    for(i=1; i<m; i++)  if(i%2 || m==4) {
      eta=(b-a)*qfint_xv[j]+a;
      sigmam+=qfint_wv[j]*func(eta);
      j++;
    }
    sigmam=2*(b-a)/m*sigmam;
    e     =mod(esm-sigmam)/2;
    esm   =0.5*(esm+sigmam);
    m=2*m;
  } while ((e>err && m<=mmax)  || m<=8);

  return  esm;
}

#endif  // complex

#endif  // ******************************************************************
