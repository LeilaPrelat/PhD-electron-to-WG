// 1-./run_bmin_vs_theta_V0.sh to run python code, bmin(V0,theta), for different bmin, xe, d/W, h/W + potential data
// 2-./split_positive_eels.sh to separate the EELS files in energies
// 3-g++ -o EELS_integrated_over_z.out EELS_integrated_over_z.cpp
// 4-./EELS_integrated_over_z.out

#include "complex.h"
#include "jga.h"
#include "spline.h"

// Constants
#define Nz_EELS 500         // maximum points for z inside EELS . this needs to be checked 
#define Nz_program 500         // maximum points for z (nz) used in c code . this needs to be checked 
#define MAX_SOLUTIONS 1000       // maximum points for solutions V0,theta for a fixed bmin 
#define MAX_ENERGY_POINTS 16000   // maximum points in energy for the eels (steps in energy of 0.005 eV)

int Ncode=400;    // points for convergency of c code . fixed
//int Nint=10*Nz_EELS;   // Integration grid
int epsilon=9;    // fix the permittivity of the substrate for all codes 
double x0_norm_w=0;  // position of the electron: 0,0.2,0.4 
int Nz_EELS_real;      // lenght of the EELS to know the zmax (last value of list_z_dat) 
int Nz_program_real;   // lenght of the c+ program. we need the exact number interpolation.          
//int Nz_code_potential_real; // lenght of the code that stores z/w V(z/w)/V0
 
// Global arrays 
extern double list_V0_sols[MAX_SOLUTIONS];
extern double list_theta_mrad_sols[MAX_SOLUTIONS];
extern double list_energy_BEM[MAX_ENERGY_POINTS];

/////////////////// from run_potential_norm_V0.sh /////////////////////////
extern double listz_norm_w_program[Nz_program];   // data from c++ program. needs to be loaded
extern double listV_norm_V0_program[Nz_program];  // data from c++ program. needs to be loaded
///////////////////////////////////////////////////////////////////////////

extern double listq_dat[Nz_EELS];  // EELS>0 data, lenght Nz
extern double listx_dat[Nz_EELS];  // EELS>0 data, lenght Nz
extern double listEnergy_dat[Nz_EELS];  // EELS>0 data, lenght Nz
extern double listz_norm_dat[Nz_EELS];  // EELS>0 data, lenght Nz
extern double listz_dat[Nz_EELS];       // EELS>0 data, lenght Nz
extern double listEELS_dat[Nz_EELS];   // EELS>0 data, lenght Nz

// interpolation of the data
double linear_interp(double *x, double *y, int n, double x_query) {
    if (x_query <= x[0]) return y[0];
    if (x_query >= x[n - 1]) return y[n - 1];

    for (int i = 0; i < n - 1; i++) {
        if (x_query >= x[i] && x_query <= x[i + 1]) {
            double t = (x_query - x[i]) / (x[i + 1] - x[i]);
            return y[i] + t * (y[i + 1] - y[i]);
        }
    }
    return 0.0; // Should never happen
}

// load potential from run_potential_norm_V0.sh
void load_potential_from_file(double wp_norm_w, double d_norm_w, double h_norm_w, double xe_norm_w){
    FILE *inputFile_potential;
    char filename_potential[120];     // Define file names for important the data (strings of lenght up to 1000)
    char labelxe_name[40];
    
    sprintf(filename_potential, "potential_data/potential_wp%.1f_d%.1f_h%.1f_N%i_xe%.1f.txt",wp_norm_w,d_norm_w,h_norm_w,Ncode,xe_norm_w);    
    inputFile_potential=fopen(filename_potential,"r");
    if (inputFile_potential == NULL) { // CHANGED: Added file open check
        printf("Error opening %s\n", filename_potential);
        return;
    }else{
//printf("opening %s\n", filename_potential);
}
     
    char header[200];
    fgets(header, sizeof(header), inputFile_potential); // Skip header

    // Read the file line by line and store data of the EELS
    int i = 0;  
    while (i < Nz_program && fscanf(inputFile_potential, "%lf %lf", &listz_norm_w_program[i], &listV_norm_V0_program[i]) == 2) {
        i++;
    }
    fclose(inputFile_potential);
     
     Nz_program_real=i;
}
    

//  ------------------------------------------------- Load EELS data -------------------------------------------------   //  
// energy_eV: photon energy in eV // a: width in nm // h: height in nm // ze: impact parameter of the electron is ze+h/2 in nm //
// xe: position of the electron in %W. possible values are [0,40,80] //  Eelectron_keV: energy of the electron in keV // N is points used in sonic
// EELS later will depend on L (lenght of the virtual geometry in nm)
void load_EELS_from_file(double energy_eV, int w, int h, double x0_norm_w, int Eelectron_keV, int N){
    FILE *inputFile;
    char filename[500]; // Define file names for important the data (strings of lenght up to 1000)
 //   char labelxe[20]; 

    int L = 500; // nm . convergency parameter for the virtual geometry
    double xe=x0_norm_w*w; // nm
//    sprintf(labelxe, "_xe%.2fW", xe);
    
    // Round energy to 7 decimal digits for filename match
    double energy_rounded = round(energy_eV * 1e7) / 1e7;
    
    // Open the file
    //sprintf(filename, "datfiles_EELS_along_z/new_files/EELS_along_z_Si_N%i_a%inm_h%inm_ze%inm_Ee%ikeV%s_energy_%.7f.txt",N,w,h,ze,Eelectron_keV,labelxe,energy_rounded);    
    sprintf(filename, "datfiles_EELS/new_files/EELS_along_z_N%i_W%inm_h%inm_L%inm_Ee%ikeV_xe%.2f_energy_%.4f.dat",N,w,h,L,Eelectron_keV,x0_norm_w,energy_rounded);    
    
    inputFile=fopen(filename,"r"); 
 if (inputFile == NULL) { // CHANGED: Added file open check
        printf("Error opening %s\n", filename);
        return;
}else{
//printf("opening %s\n", filename);
}

    // Read the file line by line and store data of the EELS
    int i = 0;
    double eV_temp, q_temp, x_temp, z_temp, EELS_temp; // x_temp should be equal to the position of the electron xe
    while (fscanf(inputFile, "%lf %lf %lf %lf %lf", &eV_temp, &q_temp, &x_temp, &z_temp, &EELS_temp) == 5) {
    if (i >= Nz_EELS) {
        printf("[ERROR] Too many data points in %s (max allowed = %d)\n", filename, Nz_EELS);
        break;
    }
     
    listz_norm_dat[i] = z_temp / w;
    listz_dat[i] = z_temp;
    listEELS_dat[i] = EELS_temp;
    i++;
    
}
    Nz_EELS_real=i;

    if (x_temp != xe) { 
     printf("[ERROR] position of electron in x = %.0f is not matching the one in bem files ( = %.0f)\n", xe, x_temp);
    }

    fclose(inputFile);

}

//  -------------  Calculate P(omega) integral over z combining EELS from sonic and potential from ./dy.out ------------- //
// zmin: minimal value of the integration over z in nm // theta_rad: theta in radians // V0: voltage in eV // 
// energy_eV: photon energy in eV // a: width in nm //  dd: distance of WG-gate normalized to a // bb: aspect ratio h/a 
// ze: impact parameter of the electron is ze+h/2 in nm // xe: position of the electron in %W. possible values are [0,40,80]
//  Eelectron_keV: energy of the electron in keV // N is points used in sonic
double P_integrated(double zmin_norm_w, double energy_eV, double w, double wp_norm_w, double h_norm_w, double d_norm_w, double x0_norm_w, int Eelectron_keV, int N, double theta_mrad, double V0, int save_P_vs_z) {
    
    double theta_rad=theta_mrad*1e-3;
    // Constants
    double me_c2_eV = 510998.95069;
    double Ee_eV = Eelectron_keV * 1e3;
    double beta = sqrt(1.0 - pow(1.0 + Ee_eV / me_c2_eV, -2.0));
    double gamma_e = 1.0 / sqrt(1 - beta * beta);
    double aux = beta * sin(theta_rad); // sin is in degree or radians in javier's codes ? 
   //  printf("sin(theta) = %.10f",  sin(theta_rad));
   double integral = 0.0;
   
    // Load data of potential along z
    load_potential_from_file(wp_norm_w,d_norm_w,h_norm_w,x0_norm_w);
 
    // Load EELS data from sonic for the array of integration over z : we obtain zmax_norm_w
    double h=h_norm_w*w;
    load_EELS_from_file(energy_eV, w, h, x0_norm_w, Eelectron_keV, N); 
    
    int Nint = Nz_EELS_real*10; // increase numbers to integrate as a sum using linear interpolation of the data
    double zmax_norm_w = listz_norm_dat[Nz_EELS_real - 1];
    double dz = (zmax_norm_w - zmin_norm_w) / (Nint - 1);
    if (dz<0){
    printf("[ERROR] some bug in the code because dz<0\n");
    }
    
    
    double list_integrand[Nint]; 
// integration over z
    for (int i = 0; i < Nint; i++) {
        double zval = zmin_norm_w + i * dz;
        
        if (zval>zmax_norm_w){
        printf("[ERROR] some bug in the code because z/w inside the integration is bigger than zmax/w\n");
        }
        
        // interpolation to evalute in zval
      //  double Vz_norm_V0 =  linear_interp(listz_norm_w_program, listV_norm_V0_program, Nz_program_real, zmin_norm_w); // THIS WITH PHI(bmin)
        double Vz_norm_V0 =  linear_interp(listz_norm_w_program, listV_norm_V0_program, Nz_program_real, zval);        // OR THIS WITH PHI(z) ????
        double EELSz = linear_interp(listz_norm_dat, listEELS_dat, Nz_EELS_real, zval);
        
        
      //  printf("z/w=%.8f, EELS(z/w)=%.8f", zval, EELSz);

        double arg = ABS(aux * aux + 2.0 * Vz_norm_V0 * V0 / (me_c2_eV * gamma_e));  
        double denom = sqrt(arg);
        
        double integrand = 2.0 * beta * EELSz / denom;

//printf("eV=%.4feV, z/w=%.7f, zmax/w = %.8f , zmin/w = %.8f\n", energy_eV, zval, zmax_norm_w, zmin_norm_w);

	//printf("arg=%.8f, denominator=%.8f, zval = %.8f, EELS=%.2e, Vz_norm_V0=%.8f, integrand=%.8f, dz * w = %.8f\n", arg, denom, zval, EELSz, Vz_norm_V0, integrand, dz * w);
	
        integral += integrand * dz * w; // because dz is in z/w
        list_integrand[i]=integrand;
    }
    
    if(save_P_vs_z == 1){
   
   char label_txt[200];
   sprintf(label_txt, "_Ee%ikeV_wp%.2f_h%.2f_d%.2f_xe%.2f", Eelectron_keV, wp_norm_w, h_norm_w, d_norm_w, x0_norm_w);
    
   char info_bmin[310];
   sprintf(info_bmin, "_V0_%.2feV_theta%.2fmrad_zmin%.2f", V0, theta_mrad, zmin_norm_w);

   char all_info_label[830];
   sprintf(all_info_label, "_energy%.4feV_w%.0fnm%s%s", energy_eV, w, label_txt, info_bmin);
   
   char filename_Pintegrand[880];
   sprintf(filename_Pintegrand, "EELS_integrated/P_integrand_vs_z%s.dat", all_info_label);
   
    // Save z and integrand values to file
    FILE *outfile;
    
    outfile = fopen(filename_Pintegrand, "w");
    
    
    if (outfile == NULL) {
        printf("Error creating output file integrand_vs_z.txt\n");
    } else {
        fprintf(outfile, "# z/w\tIntegrand for V0 = %.4f eV, theta = %.4f mrad\n", V0,theta_mrad);
        for (int i = 0; i < Nint; i++) {
            double zval = zmin_norm_w + i * dz;
            fprintf(outfile, "%.8f\t%.8e\n", zval, list_integrand[i]);  // zval is normalized by w, so z = zval * w
        }
        fclose(outfile);
        printf("[INFO] Saved integrand vs z\n");
    }

    
    }
    
return integral; 
}



