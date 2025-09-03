#ifndef jga_spline        // ************************************************
#define jga_spline 1      // ***  jga/spline.h                            ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Cubic spline and linear interpolation                             ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Example of use:                                                   ***
// ***                                                                    ***
// ***     spline  a;           declare variable a (real function)        ***
// ***     splinec a;           declare variable a (complex function)     ***
// ***     a.alloc(8,0);        allocation of 8 points (0...7); spline    ***
// ***     a.alloc(8);          like a.alloc(8,0);                        ***
// ***    (a.alloc(8,1);)       allocation of 8 points (0...7); linear    ***
// ***     a.put(3,2.1,3.2);    point 3 set to (2.1,3.2); default: (0,0)  ***
// ***     a.add(5,2.1,3.2);    add 3 (2.1,3.2) to point 5                ***
// ***     ...                                                            ***
// ***     ...                                                            ***
// ***    (a.prod(v);)          multiply function values by v             ***
// ***     a.init(yp0,ypn_1);   initialize spline (first der.)            ***
// ***                          or simply a.init() for natural spline     ***
// ***                                                                    ***
// ***     a.val(4.5)           gives y(x=4.5);  0 if it is out of range  ***
// ***     a.der(4.5)           gives y'(x=4.5); 0 if it is out of range  ***
// ***     a.dder(4.5)          gives y''(x=4.5); 0 if it is out of range ***
// ***                                                                    ***
// ***     a.free()             free memory associated to a               ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***     Alternatively one can use a.init(*x,*y,n) instead of alloc     ***
// ***     put/add and init, where x[i] and y[i] are the (x,y) values     ***
// ***     for 0<=i<n (this requires n>=4). This produces spline interp.  ***
// ***     a.init(*x,*y,n,1) produces liner interpolation.                ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Comments:                                                         ***
// ***                                                                    ***
// ***  --- yp0 and ypn_1 is the first derivative at the end points;      ***
// ***      if they are larger than 0.99*infinity, the result is natural  ***
// ***      spline: y''=0; these parameters are disregarded in the case   ***
// ***      of linear interpolation                                       ***
// ***                                                                    ***
// ***  --- constant functions are stored with n=1, and the value is      ***
// ***      v=y[0]; this can be initialized by doing:                     ***
// ***                         a.alloc(1);  a.init(v,0);                  ***
// ***                                                                    ***
// ***  --- if n<=0, the function is assumed to by 0;                     ***
// ***      initialize doing:  a.n=0; or a.alloc(0);                      ***
// ***                                                                    ***
// ***  --- for cubic spline (interpol=0) at least 4 points are needed    ***
// ***                                                                    ***
// **************************************************************************


class spline {
public:
  numero *x, *y, *y2, a,b;         // (a,b) are the extremes of the interval
  int n, interpol, init_flag;
  spline(void) {n=init_flag=0;  x=y=y2=NULL;  a=b=0;}
  void alloc(int na, int type_interpol);
  void alloc(int na) {alloc(na,0);}
  void put(int na, numero xa, numero ya) {x[na]=xa;  y[na]=ya;}
  void add(int na, numero xa, numero ya) {x[na]+=xa;  y[na]+=ya;}
  void prod(numero v) {int i;  for(i=0; i<n; i++) y[i]=y[i]*v;}
  int  init(numero yp0, numero ypn_1);
  void init(void) {if(n>1) init(infinity,infinity);}
  void init(numero *xx, numero *yy, int nn);
  void init(numero *xx, numero *yy, int nn, int type_interpol);
  void normalize(numero p);
  #ifdef jga_files
    void read(char *name, int ls, int cx, int cy, int ti);
    void read(char *name, int ls, int cx, int cy) {read(name,ls,cx,cy,0);}
  #endif
  numero integ(numero aa, numero bb, int l);
  numero integ(numero aa, numero bb) {return integ(aa,bb,0);}
  numero integ(void) {return integ(a,b);}
  numero val(numero xx);
  numero der(numero xx);
  numero dder(numero xx);
  void free(void);
  spline &operator=(const spline &ss);
  // routines for internal use

  void find(numero xx, int &ii, int &jj);  // x[ii]<=xx<x[jj], |ii-jj|=1
  numero integ(int k, numero aa, numero bb, int l);
  numero dpow(numero aa, numero bb, int l);
};

// --------------------------------------------------------------------------

void spline::free(void)
{
  if(n>0)  {delete [] x;  delete [] y;}
  if(n>=4 && interpol==0)  delete [] y2;
  n=0;
}

// --------------------------------------------------------------------------

spline &spline::operator=(const spline &ss)
{
  int i;
  n=ss.n;  interpol=ss.interpol;  init_flag=ss.init_flag;  a=ss.a;  b=ss.b;
  if(n>0) {
    x=new numero [n];  y=new numero [n];
    if(n>=4 && interpol==0) y2=new numero [n];
    for(i=0; i<n; i++) {
      x[i]=ss.x[i];  y[i]=ss.y[i];
      if(y2!=NULL && init_flag)  y2[i]=ss.y2[i];
  } }

  return *this;
}

// --------------------------------------------------------------------------

void spline::alloc(int na, int type_interpol)
{
  free();  n=na;  interpol=type_interpol;
  if(n>0)  {
    x=new numero [n];  y=new numero [n];
    for(na=0; na<n; na++) {x[na]=0; y[na]=0;}
  }
  if(n>=4 && interpol==0) y2=new numero [n];
}

// --------------------------------------------------------------------------

int spline::init(numero yp0, numero ypn_1)
{
    int i,k;  numero p,qn,sig,un;  numero *u;  init_flag=1;

    if(n==0)  return 0;                  // 0 function
    if(n==1)  {y[0]=yp0;  return 0;}     // constant function

    a=x[0];  b=x[n-1];
    if(interpol)  return 0;

    u=new numero [n];

    if(yp0>0.99*infinity)  y2[0]=u[0]=0;
    else  { y2[0]=-0.5;
            u[0]=3/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0])-yp0);}

    for(i=1; i<n-1; i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2;
      y2[i]=(sig-1)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
            -(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    if(ypn_1>0.99*infinity)  qn=un=0;
    else  {qn=0.5;
           un=3.0/(x[n-1]-x[n-2])*(ypn_1-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));}

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1);

    for(k=n-2; k>=0; k--)  y2[k]=y2[k]*y2[k+1]+u[k];

    delete [] u;   return 0;
}

// --------------------------------------------------------------------------

void spline::normalize(numero p)  // the maximum absolute value of the spline
{                                 // function is set to p
  numero val=-infinite;  int i;
  for(i=0; i<n; i++) if(ABS(y[i])>val) val=ABS(y[i]);
  for(i=0; i<n; i++) y[i]=p*y[i]/val;
  init();
}

// --------------------------------------------------------------------------

#ifdef jga_files

void spline::read(char *name, int ls, int cx, int cy, int ti)
{
  int nco=number_of_columns(ls,name);  // input columns are 1, 2, ...
  int n=number_of_rows(ls,name), i,j;  numero val, x,y;
  if(cx>nco || cy>nco || cx<1 || cy<1)
    on_error("splinec::read", "column(s) out of range in file", name);
  FILE *fin;  fin=fopen(name,"r");  skip(fin,ls);  cx--;  cy--;
  alloc(n,ti);
  for(i=0; i<n; i++) {
    for(j=0; j<nco; j++) {
      val=read_numero(fin);
      if(cx==j) x=val;  if(cy==j) y=val;
    }
    put(i,x,y);
  }
  fclose(fin);
  init();
}

#endif

// --------------------------------------------------------------------------

numero spline::val(numero xx)
{
    if(n==0)  return 0;
    if(n==1)  return y[0];
    if(n<4 && interpol==0)  on_error("spline::val", "less than 4 points");
    int  klo,khi;   numero  h,bb,aa;

    find(xx, klo,khi);
    h=x[khi]-x[klo];
    if(h==0)  on_error("spline::val", "vanishing spacing");

    if(interpol)  return y[klo] + (y[khi]-y[klo])*(xx-x[klo])/h;

    aa=(x[khi]-xx)/h;
    bb=(xx-x[klo])/h;

    return aa*y[klo]+bb*y[khi]+((aa*aa*aa-aa)*y2[klo]
                               +(bb*bb*bb-bb)*y2[khi])*h*h/6;
}

// --------------------------------------------------------------------------

numero spline::der(numero xx)
{
    if(n<=1)  return 0;
    if(n<4 && interpol==0)  on_error("spline::der", "less than 4 points");
    int  klo,khi;   numero  h,bb,aa,ap,bp;

    find(xx, klo,khi);
    h=x[khi]-x[klo];
    if(h==0)  on_error("spline::der", "vanishing spacing");

    if(interpol)  return (y[khi]-y[klo])/h;

    aa=(x[khi]-xx)/h;    ap=-1/h;
    bb=(xx-x[klo])/h;    bp= 1/h;

    return ap*y[klo]+bp*y[khi]+
           ((3*aa*aa*ap-ap)*y2[klo]+(3*bb*bb*bp-bp)*y2[khi])*h*h/6;
}

// --------------------------------------------------------------------------

numero spline::dder(numero xx)
{
    if(n<=1)  return 0;
    if(n<4 && interpol==0)  on_error("spline::dder", "less than 4 points");
    int  klo,khi;   numero  h,bb,aa,ap,bp;

    find(xx, klo,khi);
    h=x[khi]-x[klo];
    if(h==0)  on_error("spline::der", "vanishing spacing");

    if(interpol)  return 0;

    aa=(x[khi]-xx)/h;    ap=-1/h;
    bb=(xx-x[klo])/h;    bp= 1/h;

    return (aa*ap*ap*y2[klo]+bb*bp*bp*y2[khi])*h*h;
}

// --------------------------------------------------------------------------

void spline::init(numero *xx, numero *yy, int nn, int type_interpol)
{
  free();  n=nn;  interpol=type_interpol;
  if(type_interpol==0)  y2=new numero [n];
  x=copy(xx,n);  y=copy(yy,n);
  init();
}

void spline::init(numero *xx, numero *yy, int nn) {init(xx,yy,nn,0);}

// --------------------------------------------------------------------------

void spline::find(numero xx, int &klo, int &khi)
{
  int s,k;

  if(x[0]>x[n-1])  s=0;  else  s=1;
  klo=0;   khi=n-1;

  while(khi-klo>1) {
    k=(khi+klo)/2;
    if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
} }

// --------------------------------------------------------------------------

numero spline::dpow(numero aa, numero bb, int l)
{
  if(l!=-1)  return (pow(bb,l+1.0)-pow(aa,l+1.0))/(l+1);
  else       return log(bb/aa);
}

// --------------------------------------------------------------------------

numero spline::integ(int k, numero aa, numero bb, int l)
{
  numero x0=x[k], x1=x[k+1];
  numero h=x1-x0;
  numero I0=dpow(aa,bb,l);
  numero I1=dpow(aa,bb,l+1);

  if(interpol)  return y[k]*I0+(y[k+1]-y[k])*(I1-x[k]*I0)/h;

  numero I2=dpow(aa,bb,l+2);
  numero I3=dpow(aa,bb,l+3);
  numero h3=h*h*h;
  numero A1=(x1*I0-I1)/h;
  numero B1=(I1-x0*I0)/h;
  numero A3=((( I0*x1-3*I1)*x1+3*I2)*x1-I3)/h3;
  numero B3=(((-I0*x0+3*I1)*x0-3*I2)*x0+I3)/h3;

  return A1*y[k]+B1*y[k+1]+((A3-A1)*y2[k]+(B3-B1)*y2[k+1])*h*h/6;
}

// --------------------------------------------------------------------------

numero spline::integ(numero aa, numero bb, int l)
{
  if(n==0) return 0;
  if(n==1) return y[0]*dpow(aa,bb,l);

  if(n<4 && interpol==0) on_error("spline::integ", "less than 4 points");

  int i, kalo,kahi,kblo,kbhi;  find(aa,kalo,kahi);  find(bb,kblo,kbhi);

  numero  val=integ(kalo,aa,x[kahi],l);
  for(i=kahi; i<kblo; i++)  val+=integ(i,x[i],x[i+1],l);
  return  val+integ(kblo,x[kblo],bb,l);
}


// **************************************************************************
// *** complex functions                                                  ***
// **************************************************************************

#ifdef jga_complex

class splinec {
public:
  numero *x, a,b;         // (a,b) are the extremes of the interval
  complex *y, *y2;
  int n, interpol;
  splinec(void) {n=0;}
  void alloc(int na, int type_interpol);
  void alloc(int na) {alloc(na,0);}
  void put(int na, numero xa, complex ya) {x[na]=xa;  y[na]=ya;}
  void add(int na, numero xa, complex ya) {x[na]+=xa;  y[na]+=ya;}
  void prod(complex v) {int i;  for(i=0; i<n; i++) y[i]=y[i]*v;}
  int init(complex yp0, complex ypn_1);
  void init(void) {if(n>1) init(infinity,infinity);}
  void init(numero *xx, complex *yy, int nn);
  #ifdef jga_files
    void read(char *name, int ls, int cx, int cr, int ci, int ti);
    void read(char *name, int ls, int cx, int cr, int ci)
      {read(name,ls,cx,cr,ci,0);}
  #endif
  complex val(numero xx);
  complex der(numero xx);
  complex dder(numero xx);
  void free(void);
};

// --------------------------------------------------------------------------

void splinec::free(void)
{
  if(n>0)  {delete [] x;  delete [] y;}
  if(n>=4 && interpol==0)  delete [] y2;
  n=0;
}

// --------------------------------------------------------------------------

void splinec::alloc(int na, int type_interpol)
{
  free();  n=na;  interpol=type_interpol;
  if(n>0)  {
    x= new numero [n];  y= new complex [n];
    for(na=0; na<n; na++) {x[na]=0; y[na]=0;}
  }
  if(n>=4 && interpol==0) y2=new complex [n];
}

// --------------------------------------------------------------------------

#ifdef jga_files

void splinec::read(char *name, int ls, int cx, int cr, int ci, int ti)
{
  int nco=number_of_columns(ls,name);  // input columns are 1, 2, ...
  int n=number_of_rows(ls,name), i,j;  numero val, x,re,im;
  if(cx>nco || cr>nco || ci>nco || cx<1 || cr<1 || ci<1)
    on_error("splinec::read", "column(s) out of range in file", name);
  FILE *fin;  fin=fopen(name,"r");  skip(fin,ls);  cx--;  cr--;  ci--;
  alloc(n,ti);
  for(i=0; i<n; i++) {
    for(j=0; j<nco; j++) {
      val=read_numero(fin);
      if(cx==j) x=val;  if(cr==j) re=val;  if(ci==j) im=val;
    }
    put(i,x,complex(re,im));
  }
  fclose(fin);
  init();
}

#endif

// --------------------------------------------------------------------------

int splinec::init(complex yp0, complex ypn_1)
{
    int  i,k;   complex  p,qn,sig,un;   complex  *u;

    if(n==0)  return 0;                  // 0 function
    if(n==1)  {y[0]=yp0;  return 0;}     // constant function

    a=x[0];  b=x[n-1];
    if(interpol)  return 0;

    u=new complex [n];

    if(real(yp0)>0.99*infinity)  y2[0]=u[0]=0;
    else  { y2[0]=-0.5;
            u[0]=3/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0])-yp0);}

    for(i=1; i<n-1; i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2;
      y2[i]=(sig-1)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
            -(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    if(real(ypn_1)>0.99*infinity)  qn=un=0;
    else  {qn=0.5;
           un=3.0/(x[n-1]-x[n-2])*(ypn_1-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));}

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1);

    for(k=n-2; k>=0; k--)  y2[k]=y2[k]*y2[k+1]+u[k];

    delete [] u;   return 0;
}

// --------------------------------------------------------------------------

complex splinec::val(numero xx)
{
    if(n==0)  return 0;
    if(n==1)  return y[0];
    if(n<4 && interpol==0)  on_error("splinec::val", "less than 4 points");
    int  klo,khi,k,s;   numero  h,bb,aa;

    if(x[0]>x[n-1])  s=0;  else  s=1;
    klo=0;   khi=n-1;

    while(khi-klo>1) {
      k=(khi+klo)/2;
      if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
    }

    h=x[khi]-x[klo];
    if(h==0)  on_error("splinec::val", "vanishing spacing");

    if(interpol)  return y[klo] + (y[khi]-y[klo])*(xx-x[klo])/h;

    aa=(x[khi]-xx)/h;
    bb=(xx-x[klo])/h;

    return aa*y[klo]+bb*y[khi]+((aa*aa*aa-aa)*y2[klo]
                               +(bb*bb*bb-bb)*y2[khi])*h*h/6;
}

// --------------------------------------------------------------------------

complex splinec::der(numero xx)
{
    if(n<=1)  return 0;
    if(n<4 && interpol==0)  on_error("splinec::der", "less than 4 points");
    int  klo,khi,k,s;   numero  h,bb,aa,ap,bp;

    if(x[0]>x[n-1])  s=0;  else  s=1;
    klo=0;   khi=n-1;

    while(khi-klo>1) {
      k=(khi+klo)/2;
      if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
    }

    h=x[khi]-x[klo];
    if(h==0)  on_error("splinec::der", "vanishing spacing");

    if(interpol)  return (y[khi]-y[klo])/h;

    aa=(x[khi]-xx)/h;   ap=-1/h;
    bb=(xx-x[klo])/h;   bp= 1/h;

    return ap*y[klo]+bp*y[khi]+
           ((3*aa*aa*ap-ap)*y2[klo]+(3*bb*bb*bp-bp)*y2[khi])*h*h/6;
}

// --------------------------------------------------------------------------

complex splinec::dder(numero xx)
{
    if(n<=1)  return 0;
    if(n<4 && interpol==0)  on_error("splinec::dder", "less than 4 points");
    int  klo,khi,k,s;   numero  h,bb,aa,ap,bp;

    if(x[0]>x[n-1])  s=0;  else  s=1;
    klo=0;   khi=n-1;

    while(khi-klo>1) {
      k=(khi+klo)/2;
      if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
    }

    h=x[khi]-x[klo];
    if(h==0)  on_error("splinec::der", "vanishing spacing");

    if(interpol)  return 0;

    aa=(x[khi]-xx)/h;   ap=-1/h;
    bb=(xx-x[klo])/h;   bp= 1/h;

    return (aa*ap*ap*y2[klo]+bb*bp*bp*y2[khi])*h*h;
}

// --------------------------------------------------------------------------

void splinec::init(numero *xx, complex *yy, int nn)
{
  free();  n=nn;  interpol=0;
  x=copy(xx,n);  y=copy(yy,n);  y2=new complex [n];
  init();
}

#endif  // complex

#endif  // ******************************************************************
