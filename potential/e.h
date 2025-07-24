// Potential distribution near an in-plane gated structure consisting of
//    1) a central rectangular wire (width W, height h)
//    2) two identical rectangular wires to the sides of the central one
//       (widths W', same height h, surface-to-surface distance to central ribbon d)
// The bottom side of the wires are all lying on the surface (z=0) of
// a dielectric substrate (permittivity epsilon). The wires are perfect conductors
// placed on the z>0 side (vacuum), aligned along y, and with x=z=0 at the center
// of the lower surface of the central one.
// We use mirror symmetry with respect to the x=0 plane (2025-07-09 notes).
// The central wire is biased to a potential V0.
// The side wires are biased to a potential V0+U.
// All distances are in units of W.
// All potentials are in units of U.
// --- (c) NPT group, ICFO, July 9, 2025 ---
// Domain (x,y): as long as the test point (x,z) is outside all conductors. Avoid evaluating near/in conductors, which may cause artifacts due to divergence.

#include "qfint.h"
#include "rmatrix.h"

int    N=0;           // number of parametrization points
numero err=1e-12;  int nmax=2000;    // integration parameters in each interval
numero W=1,Wp,h,d;    // physical dimensions
numero epsilon;       // permittivity
numero *xs,*zs,*hs;   // discretization points and spacings
int    *tys;          // inteval label depending on which wire and ver/hor orientation
rmatrix MM,SS;        // secular matrix, induced charge

// --- calculation of the secular matrix MM
numero Mssp(numero x, numero z, numero xp, numero zp)
{
  return -log((sqr(x-xp)+sqr(z-zp))*(sqr(x+xp)+sqr(z-zp)))
         -log((sqr(x-xp)+sqr(z+zp))*(sqr(x+xp)+sqr(z+zp)))
          *(1-epsilon)/(1+epsilon);
}
numero Mij_x,Mij_z,Mij_j,Mij_vx,Mij_vz;
numero Mijint(numero ll)
{
  numero x=Mij_x, z=Mij_z;  int j=Mij_j;
  return Mssp(x,z,xs[j]+ll*Mij_vx,zs[j]+ll*Mij_vz);
}
numero Mij(numero x, numero z, int j)
{
  Mij_x=x;  Mij_z=z;  Mij_j=j;
  if(tys[j]==1 || tys[j]==2) {Mij_vx=0; Mij_vz=1;}
  else                       {Mij_vx=1; Mij_vz=0;}
  return qfint(-hs[j]/2,0,Mijint,err,nmax)
        +qfint(0,hs[j]/2, Mijint,err,nmax);
}
numero Mij(int i, int j) {return Mij(xs[i],zs[i],j);}

void set_MM(numero Wp_, numero d_, numero h_, int N_, numero epsilon_)
{
  if(N>0) {delete [] xs;  delete [] zs;  delete [] hs;  delete [] tys;}
  Wp=Wp_*W;  d=d_*W;  h=h_*W;  N=N_;  epsilon=epsilon_;
  numero L=W+3*h+2*Wp;
  int NW=(W/2)/L*N+1, Nh=h/L*N+1, Np=Wp/L*N+1, Nt=2*(NW+Np)+3*Nh, N=Nt; // #points in W/2,h,Wp
  numero hW=(W/2)/NW, hh=h/Nh, hp=Wp/Np;
  xs=new numero[N];  zs=new numero[N];  hs=new numero[N];  tys=new int[N];
  int i,j,tN=0;  numero tl=0;

  // --- discretization points and lengths of intervals
  for(j=0; j<NW; j++) {xs[j+tN]=(j+0.5)*hW;  zs[j+tN]=0;           hs[j+tN]=hW;  tys[j+tN]=0;}  tN+=NW;
  for(j=0; j<NW; j++) {xs[j+tN]=(j+0.5)*hW;  zs[j+tN]=h;           hs[j+tN]=hW;  tys[j+tN]=0;}  tN+=NW;
  for(j=0; j<Nh; j++) {xs[j+tN]=W/2;         zs[j+tN]=(j+0.5)*hh;  hs[j+tN]=hh;  tys[j+tN]=1;}  tN+=Nh;
  for(j=0; j<Nh; j++) {xs[j+tN]=W/2+d;       zs[j+tN]=(j+0.5)*hh;  hs[j+tN]=hh;  tys[j+tN]=2;}  tN+=Nh;
  for(j=0; j<Nh; j++) {xs[j+tN]=W/2+d+Wp;    zs[j+tN]=(j+0.5)*hh;  hs[j+tN]=hh;  tys[j+tN]=2;}  tN+=Nh;
  for(j=0; j<Np; j++) {xs[j+tN]=W/2+d+(j+0.5)*hp;  zs[j+tN]=0;     hs[j+tN]=hp;  tys[j+tN]=3;}  tN+=Np;
  for(j=0; j<Np; j++) {xs[j+tN]=W/2+d+(j+0.5)*hp;  zs[j+tN]=h;     hs[j+tN]=hp;  tys[j+tN]=3;}  tN+=Np;

  MM.alloc(N+1,N+1);    // eigenmatrix
  for(i=0; i<N; i++)
  for(j=0; j<N; j++) MM.a(i,j,Mij(i,j));
  for(j=0; j<N; j++) {MM.a(N,j,1.0); MM.a(j,N,-1.0);}
  --MM;
  rmatrix VV(N+1,1);    // potentials
  //for(j=0; j<N; j++) if(tys[j]>=2) VV.a(j,0,1.0); // ONLY BOUNDARY CONDITIONS FOR THE TWO SIDE PEC
  // Not a constant value inside the middle PEC
  
  // Explicit boundary condition enforcement
   for(j = 0; j < N; j++) {
     if (tys[j] == 0 || tys[j] == 1) {
         VV.a(j, 0, -1.0);  // Set central wire to -V0/U
     } else if (tys[j] >= 2) {
         VV.a(j, 0, 1.0); // Side wires to V0/U
    }
 }
 
  // Explicit boundary condition enforcement
//for (j = 0; j < N; j++) {
//    if (tys[j] == 0 || tys[j] == 1) {
 //       for (int k = 0; k < N + 1; k++) MM.a(j, k, 0.0);  // Clear row
  //      MM.a(j, j, -1.0);  // Set diagonal
//        VV.a(j, 0, -1.0);   // Set corresponding value in VV
 //   } else if (tys[j] >= 2) {
  //      for (int k = 0; k < N + 1; k++) MM.a(j, k, 0.0);
   //     MM.a(j, j, 1.0);
    //    VV.a(j, 0, 1.0);
   // }
//}

  SS=MM*VV;
//  printf("V0/U=%g\n",SS.a(N+1));
//numero test=0;
//for(j=0; j<N; j++) test+=SS.a(j,0);
//  printf("sum=%g\n",test);
}

// --- total potential (in units of U)
numero potential(numero x, numero z) // for (x,z) outside the wires and z>0
{
  numero val=0;  int j;
  for(j=0; j<N; j++) val+=Mij(x,z,j)*SS.a(j,0);
  return val;
}
