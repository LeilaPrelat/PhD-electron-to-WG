// Charge distribution on (and potential around) a biased
// perfect-conductor wire of rectangular cross section with rounding edges
// placed near a planar perfect-conductor surface.
// We assume that the wire is at a potential V0 relative to the plane.
// The wire center is at (0,0), and the rectangle sides are along x and y
// with full lengths a and b, respectively (a along x, b along y).
// The planar surface is parallel to the x-z plane, at a distance d from
// the closest x-parallel surface of the wire.
// All lengths (including b) are given in units of the x side-length a.
// All potentials are in units of V0. The plane is set at 0 potential.
// --- (c) F. Javier Garcia de Abajo, ICFO, May 27, 2025 ---

#include "qfint.h"
#include "rmatrix.h"

numero ss;  // corner rounding radius (external param. in units of a)
int    N=0; // number of parametrization points (external parameter)
numero err=1e-12;  int nmax=2000;    // integration parameters for neighbors
numero L,ds,*xs,*ys;  // perimeter, element length, discretization points
numero bb;            // side length along y divided by side length along x
numero dd;            // wire-plane distance in units of a
rmatrix MM,nn;        // secular matrix, induced charge

// --- surface contour coordinates as a function of length s along perimeter,
void Rs(numero s, numero &xs, numero &ys)    // the rounding radius is ss,
{                                            // center at (0,0), sides || x,y,
  if(s>L/2) {Rs(L-s,xs,ys); ys=-ys;} else {  // s=0 and L at (xs,ys)=(-0.5,0);
    if(s>L/4) {Rs(L/2-s,xs,ys); xs=-xs;} else {  // counter-clockwise param. starting from middle of bottom side.
      numero s1=0.5*bb-ss;
      if(s<=s1) {xs=-0.5;  ys=-s;} else {
        s=s-s1;  s1=pi*ss/2;
        if(s<=s1) {xs=-0.5+ss*(1-cos(s/ss));  ys=-0.5*bb+ss*(1-sin(s/ss));}
        else {s=s-s1; xs=-0.5+ss+s;  ys=-0.5*bb;}
} } } }

// --- determine whether point (x,y) is inside (1) or outside (0) the wire
// Uses symmetry (mirrors all points to the first quadrant).
int isinside(numero x, numero y)
{
  if(x<0) x=-x;  if(y<0) y=-y;
  if(x>0.5) return 0;
  if(y>0.5*bb) return 0;
  if(x<=0.5-ss) return 1;
  if(y<=0.5*bb-ss) return 1;
  if(sqr(x-(0.5-ss))+sqr(y-(0.5*bb-ss))>ss*ss) return 0;
  return 1;
}

// --- calculation of the secular matrix MM
numero Mssp(numero s, numero sp)
{
  numero x1,y1,x2,y2;  Rs(s,x1,y1);  Rs(sp,x2,y2);
  return log((sqr(x1-x2)+sqr(y1-y2))/(sqr(x1-x2)+sqr(y1+y2+bb+2*dd)));
}
numero Mij_s;  numero Mijint(numero sp) {return Mssp(Mij_s,sp);}
numero Mij(int i, int j)
{
  Mij_s=(i+0.5)*ds;  numero sj=(j+0.5)*ds;
  if(ABS(i-j)%N>20) return Mssp(Mij_s,sj)*ds; // middle point if distance is large
  if(i!=j) return qfint(sj-ds/2,sj+ds/2,Mijint,err,nmax); // for close neighbors : accurate integration
  return qfint(sj-ds/2,sj,Mijint,err,nmax)
        +qfint(sj,sj+ds/2,Mijint,err,nmax);
}
void set_MM(numero ss_, numero dd_, int N_)
{
  if(N>0) {delete [] xs;  delete [] ys;}
  N=N_;  if(N%4) N=(N/4)*4;                   // make N multiple of 4
  ss=ss_;  L=2*(1+bb)-8*ss+2*pi*ss;  ds=L/N;  // L=perimeter; ds=element length
  dd=dd_;
  xs=new numero[N]; ys=new numero[N];
  int i,j;
  for(i=0; i<N; i++) Rs((i+0.5)*ds,xs[i],ys[i]);  // parametrization coodinates
  MM.alloc(N,N);
  for(i=0; i<N; i++)
  for(j=0; j<N; j++) MM.a(i,j,Mij(i,j));
  --MM;
  rmatrix V0(N,1);                            // we will calculate V(x,y)/V0
  for(i=0; i<N; i++) V0.a(i,0,1.0);
  nn=MM*V0;
}

// --- total potential (in units of V0)

numero Msp_xx,Msp_yy;
numero V_reference = -1;  // used to normalize the potential

numero Msp(numero sp)
{
  numero x=Msp_xx,y=Msp_yy, xsp,ysp;  Rs(sp,xsp,ysp);
  return log((sqr(x-xsp)+sqr(y-ysp))/(sqr(x-xsp)+sqr(y+ysp+bb+2*dd)));
}
numero potential(numero x, numero y)
{  
  if(y<-bb/2-dd) return 0;
  if(isinside(x,y)) return 1;
  numero val=0, si;  int i;
  Msp_xx=x;  Msp_yy=y;
  for(i=0; i<N; i++) {
    si=(i+0.5)*ds;
    val-=nn.a(i,0)*qfint(si-ds/2,si+ds/2,Msp,err,nmax);
  }
  
  
   // Normalize once using reference point on the conductor
   // You need to compute the potential relative to a known point on the conductor where the potential should be exactly 1 (e.g., the center of the bottom   edge: x=0,y=−bb/2). 
   // This avoids inconsistency due to a mis-scaled Green's function or numerical artifacts.
   // We fix this by dividing the whole potential by the value it gives at a known point on the conductor. This way, everything is scaled correctly.
  if (V_reference < 0) {
    numero xref = 0.0;
    numero yref = -bb / 2 + 1e-6;  // Just inside the conductor
    numero val_ref = 0;
    Msp_xx = xref;
    Msp_yy = yref;
    for (int i = 0; i < N; i++) {
      si = (i + 0.5) * ds;
      val_ref -= nn.a(i, 0) * qfint(si - ds/2, si + ds/2, Msp, err, nmax);
    }
    V_reference = val_ref;
  }

  return val / V_reference;
}


// converged when comparing "a.out 2 0.1 1 N > bo.dat" with N=200 and 400

// --- --- ---
