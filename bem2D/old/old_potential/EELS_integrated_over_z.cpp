#include "qfint.h"
#include "rmatrix.h"

// Declare the function prototype (you’d replace this with an actual header normally)
extern double P_integrated(double zmin, double theta_rad, double V0, double energy_eV, int a, double dd, double bb, int ze, int xe, int Eelectron_keV, int N);

int main() {
    // Example input parameters:
    double zmin = 0.0;             // Starting z value in normalized units (z/a)
    double theta_rad = 0.0;        // Scattering angle in radians
    double V0 = 100.0;             // Potential in eV
    double energy_eV = 1.5;        // Photon energy in eV
    int a = 20;                    // Width in nm
    double dd = 0.2;               // WG-gate distance normalized to a
    double bb = 2.0;               // Aspect ratio h/a
    int ze = 5;                    // Impact parameter
    int xe = 40;                   // Electron x-position in %W
    int Eelectron_keV = 200;       // Electron energy in keV
    int N = 400;                   // Points used in Sonic/EELS calculation

    double result = P_integrated(zmin, theta_rad, V0, energy_eV, a, dd, bb, ze, xe, Eelectron_keV, N);

    printf("Integrated probability P(omega): %g\n", result);

    return 0;
}

