#include "complex.h"
#include "jga.h"
#include "spline.h"

int Nz_program; // depends on what is inside dy.cpp 
int Nz=400; // points for z inside EELS . this needs to be checked 
int Nint = 10 * Nz; // Integration grid

// Global arrays
numero listz_norm_program;  // data from c++ program
numero listV_norm_program;  // data from c++ program 
numero listz_norm_dat; // EELS data, lenght Nz
numero listEELS_dat;   // EELS data, lenght Nz
numero listEnergy_dat  // EELS data, lenght Nz
numero listq_dat       // EELS data, lenght Nz
numero listx_dat       // EELS data, lenght Nz
numero listz_norm_dat  // EELS data, lenght Nz

// Always load dy.out and fill the data arrays
void load_dy(double bb, double ss, double dd, int N) {
    FILE *potential = fopen("./dy.out", "r");
    if (!potential) {
        fprintf(stderr, "Failed to open dy.out\n");
        exit(1);
    }

    numero listz_norm, listV_norm;
    int i = 0;
    while (fscanf(potential, "%lf %lf\n", &listz_norm, &listV_norm) == 2) {
        listz_norm_program[i] = listz_norm;
        listV_norm_program[i] = listV_norm;
        i++;
    }
    fclose(potential);

    Nz_program = i;
}


//  ------------------------------------------------- Load EELS data -------------------------------------------------   //  
// energy_eV: photon energy in eV // a: width in nm // h: height in nm // ze: impact parameter of the electron is ze+h/2 in nm //
// xe: position of the electron in %W. possible values are [0,40,80] //  Eelectron_keV: energy of the electron in keV // N is points used in sonic
void load_EELS_from_file(double energy_eV, int w, int h, int xe, int ze, int Eelectron_keV, int N){
    FILE *inputFile;
    char filename[256]; // Define file names for important the data (strings of lenght up to 1000)
    char labelxe[20]; 

    // Define more variables: empty arrays to store data from the files
    numero *listEnergy, *listq, *listx, *listz, *listEELS;
    
   // Corrected labelxe
    if (xe == 40) {
        sprintf(labelxe, "_xe40W");
    } else if (xe == 80) {
        sprintf(labelxe, "_xe80W");
    } else {
        sprintf(labelxe, "_xe0");
    }
    
    // Open the file
    sprintf(filename, "datfiles_EELS_along_z/EELS_along_z_Si_N%i_a%inm_h%inm_ze%inm_Ee%ikeV%s_energy%.8f.dat",N,w,h,ze,Eelectron_keV,labelxe,energy_eV);    
    printf("%s\n",filename);
    inputFile=fopen(filename,"r"); 

    // Read the file line by line and store data of the EELS
    int i=0;
    while ( fscanf(inputFile, "%lf %lf %lf %lf %lf\n", &listEnergy, &listq, &listx, &listz, &listEELS)==5) {   
    listEnergy_dat[i]=listEnergy; listq_dat[i]=listq; listx_dat[i]=listx; listz_norm_dat[i]=listz/a; listEELS_dat[i]=listEELS; // save the z normalized to a
    i++;
}
    fclose(inputFile);

}
 

//  -------------  Calculate P(omega) integral over z combining EELS from sonic and potential from ./dy.out ------------- //
// zmin: minimal value of the integration over z in nm // theta_rad: theta in radians // V0: voltage in eV // 
// energy_eV: photon energy in eV // a: width in nm //  dd: distance of WG-gate normalized to a // bb: aspect ratio h/a 
// ze: impact parameter of the electron is ze+h/2 in nm // xe: position of the electron in %W. possible values are [0,40,80]
//  Eelectron_keV: energy of the electron in keV // N is points used in sonic
double P_integrated(double zmin, double energy_eV, double theta_rad, double V0, int w, double bb, int xe, int ze,
    int Eelectron_keV, int N) {
    
    // Constants
    double me_c2_eV = 510998.95069;
    double Ee_eV = Eelectron_keV * 1e3;
    double beta = sqrt(1.0 - pow(1.0 + Ee_eV / me_c2_eV, -2.0));
    double gamma_e = 1.0 / sqrt(1 - beta * beta);
    double aux = beta * sin(theta_rad);
    
    double ss=0.1*a;// fixed value
    
    // load listz_norm_program and listV_norm_program
    load_dy_if_needed(bb, ss, dd, N); 

    // Interpolate potential from ./dy.out
    spline s1;
    s1.alloc(Nint, 1);               // allocate space for n points, degree 1 = linear
    for (int i = 0; i < Nz_program; i++)
      s1.put(i, listz_norm_program[i], listV_norm_program[i]);
    s1.init();                   // prepares internal structures
    
    
    // Load EELS data from sonic
    double h;
    h=bb*w;
    load_EELS_from_file(energy_eV, w, h, ze, xe, Eelectron_keV, N); 
    
    // Interpolate EELS
    spline s2;
    s2.alloc(Nint, 1);               // allocate space for n points, degree 1 = linear
    for (int i = 0; i < Nz; i++)
      s2.put(i, listz_norm_dat[i], listEELS_dat[i]);
    s2.init();                   // prepares internal structures
    
    
    double zmax = listz_norm_dat[Nz - 1];
    double dz = (zmax - zmin) / (Nint - 1);

    double integral = 0.0;
    numero list_integrand[Nint]; 

    for (int i = 0; i < Nint; i++) {
        double zval = zmin + i * dz;
        double Vz_norm = s1.val(zval); 
        double EELSz = s2.val(zval);

        double arg = ABS(aux * aux + 2.0 * Vz_norm * V0 / (me_c2_eV * gamma_e));
        double denom = sqrt(arg);
        double integrand = 2.0 * beta * EELSz / denom;

        integral += integrand * dz;
        list_integrand[i]=integrand;
    }

    return integral * w; // because dz is in z/w
}

