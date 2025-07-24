#ifndef jga_rmatrix       // ************************************************
#define jga_rmatrix 1     // ***  jga/rmatrix.h                           ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Matrix algebra (real matrices)                                    ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Operations:                                                       ***
// ***                                                                    ***
// ***    For a nxm matrix,  0<=i<n  and  0<=j<m                          ***
// ***    A*x  x*A  A/x  -A              A(nxm)                           ***
// ***    A+B  A-B  A=B                  A(nxm) and B(nxm)                ***
// ***    x+A  A+x  A-x  x-A  A=x  x/A   A(nxn)                           ***
// ***    A*B                            A(nxm) and B(mxq)                ***
// ***    A/B                        <-  A(nxm) and B(mxm)                ***
// ***                                                                    ***
// ***    prod(A,B)                      A(nx1), B(nx1), return numero.   ***
// ***                                                                    ***
// ***    --A                        <-  A(nxn) inverse, A.sdet           ***
// ***                                                                    ***
// ***    A.a(i,j,x)                 <-  A(i,j) = x                       ***
// ***    A.a(i,j)                   <-  gives A(i,j)                     ***
// ***    A.a(i,x)                   <-  like A.a(i,0,x)                  ***
// ***    A.a(i)                     <-  like A.a(i,0)                    ***
// ***                                                                    ***
// ***    A.add(i,j,x)               <-  add x to A.a(i,j)                ***
// ***    A.add(i,x)                 <-  add x to A.a(i,0)                ***
// ***    A.alloc(n,m)               <-  memory allocation (nxm)          ***
// ***    A.alloc(n)                 <-  like A.alloc(n,1)                ***
// ***                                                                    ***
// ***    A.free()                   <-  memory de-allocation             ***
// ***                                                                    ***
// ***  Functions:                                                        ***
// ***                                                                    ***
// ***    transpose(A) = A^T         <-  A(nxm)                           ***
// ***                                                                    ***
// ***  Here A and B are 'rmatrix', x is a number of type 'numero'.       ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** declaration                                                        ***
// **************************************************************************

class rmatrix {
public:

  int n,m,nm, sel, sdet;                 // sel=0  -> general
  numero *p;                             // sel=1  -> diagonal

  rmatrix(void) {n=m=nm=sel=0; sdet=1;}
  void free(void)  {if(nm) {delete [] p;}  n=m=nm=0;}
  rmatrix(int i, int j)  {nm=0;  alloc(i,j);}
  rmatrix(int i) {nm=0;  alloc(i,1);}
  rmatrix(const rmatrix &A);
  ~rmatrix(void) {free();}
  void alloc(int i, int j);
  void alloc(int i)  {alloc(i,1);}
  void alloc_diagonal(int i) {alloc(i,1);  sel=1;}
  void a(int i, int j, numero x) {p[i*m+j]=x;}
  void a(int i, numero x) {p[i]=x;}
  void add(int i, int j, numero x) {p[i*m+j]+=x;}
  void add(int i, numero x) {p[i]+=x;}
  numero a(int i, int j);
  numero a(int i) {if(i>-1 && i<nm)  return p[i];  else  return 0;}
  rmatrix &operator=(numero x);
  rmatrix &operator=(const rmatrix &A);
  rmatrix &operator--(void);                // rmatrix inversion
  numero gaussj(const rmatrix &Mb, int n, int m, int opt);
  numero gaussj(const rmatrix &Mb, int n, int m)   // return ln(||Mb||)
    {return gaussj(Mb,n,m,0);}
//------------------------------------------------------------
  friend rmatrix transpose(rmatrix A);
  friend rmatrix operator+(numero x, rmatrix A);
  friend rmatrix operator+(rmatrix A, numero x) {return x+A;}
  friend rmatrix operator+(rmatrix A, const rmatrix &B);
  friend rmatrix operator-(rmatrix A);
  friend rmatrix operator-(numero x, rmatrix A);
  friend rmatrix operator-(rmatrix A, numero x) {return (-x)+A;}
  friend rmatrix operator-(rmatrix A, const rmatrix &B);
  friend rmatrix operator*(numero x, rmatrix A);
  friend rmatrix operator*(rmatrix A, numero x) {return x*A;}
  friend rmatrix operator*(const rmatrix &A, const rmatrix &B);
  friend rmatrix operator*(numero *r, rmatrix A);
  friend rmatrix operator*(const rmatrix &A, numero *r) {return r*A;}
  friend numero prod(rmatrix A, rmatrix B);
  friend rmatrix operator/(numero x, rmatrix A) {return x*(--A);}
  friend rmatrix operator/(rmatrix A, numero x) {return (1/x)*A;}
  friend rmatrix operator/(rmatrix A, rmatrix B) {return A*(--B);}
};


// **************************************************************************
// *** implementation                                                     ***
// **************************************************************************

rmatrix::rmatrix(const rmatrix &A)
{
  int ij;
  nm=0;  if(A.nm) alloc(A.n,A.m);  sel=A.sel;

  for(ij=0; ij<nm; ij++)  p[ij]=A.p[ij];
}

numero rmatrix::a(int i, int j)
{
  if(i>-1 && i<n) {
    if(sel==1)  if(i==j)          return p[i];
    if(sel==0)  if (j>-1 && j<m)  return p[i*m+j];
  }
  return 0;
}

void rmatrix::alloc(int i, int j)
{
  if(nm)  delete [] p;   n=i; m=j; nm=i*j;  sel=0;
  if(nm) {p=new numero [nm];  for(i=0; i<nm; i++) p[i]=0;}
}

rmatrix operator+(numero x, rmatrix A)
{
  int i,j,ij;

  if(A.sel==0) {
    if(A.n!=A.m)  on_error("rmatrix","A is not square in x+A");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)  if(i==j)  A.p[ij]+=x;
  } else

  if(A.sel==1)  for(i=0; i<A.n; i++)  A.p[i]+=x;

  return A;
}

rmatrix operator+(rmatrix A, const rmatrix &B)
{
  int i,j,ij;

  if(A.sel==0 && B.sel==0) {
    if(A.n!=B.n || A.m!=B.m) on_error("rmatrix","diff. dimensions in A+B (1)");
    for(ij=0; ij<A.nm; ij++)  A.p[ij]+=B.p[ij];
  } else

  if(A.sel==0 && B.sel==1) {
    if(A.n!=B.n || A.m!=B.n) on_error("rmatrix","diff. dimension in A+B (2)");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]+=B.p[i];
  } else

  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n) on_error("rmatrix","different dimension in A+B (3)");
    for(i=0; i<A.n; i++)  A.p[i]+=B.p[i];
  } else

  if(A.sel==1 && B.sel==0)  return B+A;

  return A;
}

rmatrix operator-(rmatrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=-A.p[ij];
  return A;
}

rmatrix operator-(numero x, rmatrix A)
{
  int i,j,ij;

  if(A.sel==0) {
    if(A.n!=A.m)  on_error("rmatrix","A is not square in x-A");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]=x-A.p[ij];  else  A.p[ij]=-A.p[ij];
  } else

  if(A.sel==1)  for(i=0; i<A.n; i++)  A.p[i]=x-A.p[i];

  return A;
}

rmatrix operator-(rmatrix A, const rmatrix &B)
{
  int i,j,ij;

  if(A.sel==0 && B.sel==0) {
    if(A.n!=B.n || A.m!=B.m) on_error("rmatrix","diff. dimensions in A-B (1)");
    for(ij=0; ij<A.nm; ij++)  A.p[ij]-=B.p[ij];
  } else

  if(A.sel==0 && B.sel==1) {
    if(A.n!=B.n || A.m!=B.n) on_error("rmatrix","diff. dimension in A-B (2)");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]-=B.p[i];
  } else

  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n) on_error("rmatrix","different dimension in A-B (3)");
    for(i=0; i<A.n; i++)  A.p[i]-=B.p[i];
  } else

  if(A.sel==1 && B.sel==0)  return -(B-A);

  return A;
}

rmatrix operator*(numero x, rmatrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=x*A.p[ij];
  return A;
}

rmatrix operator*(const rmatrix &A, const rmatrix &B)
{
  int i,j,k,ij;  numero val;  rmatrix C;
  numero *Cp,*Ap0,*Ap,*Bp,*Bp0;

  if(A.sel==0 && B.sel==0) {
    if(A.m!=B.n)  on_error("rmatrix","incompatible dimensions in A*B");
    C.alloc(A.n,B.m);
    for(i=0, Cp=C.p, Ap0=A.p; i<A.n; i++, Ap0+=A.m)
    for(j=0, Bp0=B.p; j<B.m; j++, Cp++, Bp0++) {
      for(k=0, val=0,  Ap=Ap0, Bp=Bp0;  k<A.m; k++, Ap++, Bp+=B.m)
        val+=(*Ap)*(*Bp);
      *Cp=val;
  } } else

  if(A.sel==0 && B.sel==1) {
    if(A.m!=B.n)  on_error("rmatrix","incompatible dimensions in A*B, B=diag.");
    C.alloc(A.n,A.m);
    for(i=ij=0; i<A.n; i++) for(j=0; j<A.m; j++,ij++) C.p[ij]=A.p[ij]*B.p[j];
  } else

  if(A.sel==1 && B.sel==0) {
    if(A.n!=B.n)  on_error("rmatrix","incompatible dimensions in A*B, A=diag.");
    C.alloc(B.n,B.m);
    for(i=ij=0; i<B.n; i++) for(j=0; j<B.m; j++,ij++) C.p[ij]=A.p[i]*B.p[ij];
  } else

  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n)  on_error("rmatrix","incompatible diag. dimensions in A*B");
    C.alloc_diagonal(A.n);
    for(i=0; i<A.n; i++) C.p[i]=A.p[i]*B.p[i];
  }

  return C;
}

rmatrix transpose(rmatrix A)
{
  if(A.sel==1)  return A;

  rmatrix C(A.m,A.n);
  int i,j,ij;

  for(ij=0; ij<A.nm; ij++) {
    j=ij%A.m;  i=ij/A.m;
    C.p[j*C.m+i]=A.p[ij];
  }

  return C;
}

rmatrix &rmatrix::operator--(void)        // rmatrix inversion and return ln|A|
{
  int ij;  if(sel==1) {for(ij=0; ij<nm; ij++) p[ij]=1/p[ij];  return *this;}
  rmatrix B;
  gaussj(B,n,0);

  return *this;
}

rmatrix &rmatrix::operator=(const rmatrix &A)
{
  int ij;
  if(A.nm)  alloc(A.n,A.m);  sel=A.sel;
  for(ij=0; ij<nm; ij++)  p[ij]=A.p[ij];
  return *this;
}

rmatrix &rmatrix::operator=(numero x)
{
  int i,j,ij;

  if(sel==0)
    for(i=ij=0; i<n; i++) for(j=0; j<m; j++,ij++)
      if(i==j) p[ij]=x; else  p[ij]=0;

  else if(sel==1) for(i=0; i<n; i++) p[i]=x;

  return *this;
}

rmatrix operator*(numero *r, rmatrix A)
{
  int ij;  numero *p;  p=A.p;

  for(ij=0; ij<A.nm; ij++, r++, p++)  *p=(*r)*(*p);

  return A;
}

numero prod(rmatrix A, rmatrix B)
{
  numero val=0;
  int i;

  if(A.m!=1 || B.m!=1 || A.n!=B.n) on_error("rmatrix","A*B is not a scalar");

  for(i=0; i<A.n; i++)  val+=A.p[i]*B.p[i];

  return val;
}

numero rmatrix::gaussj(const rmatrix &Mb, int nn, int mm, int opt)
{
          if(n!=m || n<nn)  on_error("rmatrix","A not square in gaussj(A)");
  if(mm)  if(m!=Mb.n || Mb.m<mm)
            on_error("rmatrix","incompatible dimensions in gaussj");

  int  i,icol,irow,j,k,l,ll, *indxc,*indxr,*ipiv,*iA,*iB;  sdet=1;
  numero  big,val;
  numero  dum,pivinv, *pirow,*picol,*pp,*ppp,*pa,*pb,*ppb, det=0;

  ipiv=new int [nn];  indxc=new int [nn];  indxr=new int [nn];
  iA=new int [nn];  iB=new int [nn];
  for(i=0; i<nn; i++)  {ipiv[i]=0;  iA[i]=i*m;  iB[i]=i*Mb.m;}

  for(i=0; i<nn; i++) {                    // main loop over columns
    big=0;
    for(j=0,pp=p; j<nn; j++,pp+=m) {       // search pivot element
      if(ipiv[j]!=1)
      for(k=0,ppp=pp; k<nn; k++,ppp++) {
        if(!ipiv[k]) {
          val=ABS(*ppp);
          if(val>=big) {
            big=val;
            irow=j;
            icol=k;
        } }  else  if(ipiv[k]>1)
                     on_error("rmatrix","1) singular rmatrix in gaussj");
    } }
    ipiv[icol]=ipiv[icol]+1;
    if(irow!=icol) {
      sdet=-sdet;
      pirow=p+iA[irow];  picol=p+iA[icol];
      for(l=0; l<nn; l++,pirow++,picol++) {
        dum=*pirow;  *pirow=*picol;  *picol=dum;
      }
      if(mm)  {
        pirow=Mb.p+iB[irow];  picol=Mb.p+iB[icol];
        for(l=0; l<mm; l++,pirow++,picol++) {
          dum=*pirow;  *pirow=*picol;  *picol=dum;
    } } }
    indxr[i]=irow;
    indxc[i]=icol;   pp=p+iA[icol]+icol;
    if(*pp==0) {
      if(opt)  return -infinity;
      else     on_error("rmatrix","2) singular rmatrix in gaussj");
    }
    det=det+log(ABS(*pp));  if(*pp<0) sdet=-sdet;
    pivinv=1/(*pp);
    *pp=1;
    for(l=0, pp=p+iA[icol]; l<nn; l++,pp++)  *pp=(*pp)*pivinv;
    if(mm)  for(l=0, pp=Mb.p+iB[icol]; l<mm; l++,pp++)  *pp=(*pp)*pivinv;
    for(ll=0,pp=p,ppb=Mb.p; ll<nn; ll++,pp+=m,ppb+=mm)  if(ll!=icol)  {
        dum=*(pp+icol);
        *(pp+icol)=0;  pa=p+iA[icol];  pb=Mb.p+iB[icol];
        for(l=0,ppp=pp; l<nn; l++,ppp++,pa++)  *ppp=*ppp-(*pa)*dum;
        if(mm)  for(l=0,ppp=ppb; l<mm; l++,ppp++,pb++)  *ppp=*ppp-(*pb)*dum;
    }
  }

  for(l=nn-1; l>=0; l--)  if(indxr[l]!=indxc[l])
  for(k=0,pp=p+indxr[l],ppp=p+indxc[l]; k<nn; k++,pp+=m,ppp+=m) {
      dum=*pp;  *pp=*ppp;  *ppp=dum;
  }

  delete [] ipiv;  delete [] indxc;  delete [] indxr;
  delete [] iA;  delete [] iB;

  return det;
}

#endif  // ******************************************************************
