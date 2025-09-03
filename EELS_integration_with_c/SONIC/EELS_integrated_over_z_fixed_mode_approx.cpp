#include "complex.h"
#include "jga.h"
#include "EELS_integrated_over_z_fixed_mode_approx.h"
 
double listz_norm_w_program[Nz_program];   // data from c++ program. needs to be loaded
double listV_norm_V0_program[Nz_program];  // data from c++ program. needs to be loaded
///////////////////////////////////////////////////////////////////////////

double listq_dat[Nz_EELS];         // EELS>0 data, lenght Nz
double listx_dat[Nz_EELS];        // EELS>0 data, lenght Nz
double listEnergy_dat[Nz_EELS];  // EELS>0 data, lenght Nz
double listz_norm_dat[Nz_EELS];  // EELS>0 data, lenght Nz
double listz_dat[Nz_EELS];       // EELS>0 data, lenght Nz
double listEELS_dat[Nz_EELS];   // EELS>0 data, lenght Nz
 
// importing EELS integrated over z
double P_integrated_approx(double zmin, double energy_eV, int w, double wp_norm_w, double h_norm_w, double x0_norm_w, int Eelectron_keV, int N, double theta_mrad, double V0);


// Global arrays init
double list_V0_sols[MAX_SOLUTIONS]; // solutions for bmin. needs to be loaded 
double list_theta_mrad_sols[MAX_SOLUTIONS]; // solutions for bmin. needs to be loaded  
double list_energy_BEM[MAX_ENERGY_POINTS]; // list of energy to iterate 

int Nsolutions;      // Will be set when loading solutions. starts from 1 because of the headers
int num_energy_points=0; // number of energies from bem2D with EELS>0

// load solutions V0, theta_mrad from run_potential_norm_V0.sh : maybe the solutions are not needed 
void load_solutions_V0_theta_from_file(double wp_norm_w, double d_norm_w, double h_norm_w, double x0_norm_w, int Ee_electron_keV, double bmin_norm_w) {
    FILE *inputFile_V0_theta_mrad_solutions;
    char filename_V0_theta[256]; // Safer length
    char label_txt[200];

    // Define temporary variables for reading
    double list_V0_load_data, list_theta_mrad_load_data;

    sprintf(label_txt, "_for_bmin%.2f_Ee%ikeV_wp%.2f_h%.2f_d%.2f_xe%.1f", bmin_norm_w, Ee_electron_keV, wp_norm_w, h_norm_w, d_norm_w, x0_norm_w);
    sprintf(filename_V0_theta, "potential_data/V0_theta_mrad_sols%s.txt", label_txt);
    inputFile_V0_theta_mrad_solutions = fopen(filename_V0_theta, "r");
    if (inputFile_V0_theta_mrad_solutions == NULL) {
        printf("Error: Cannot open file %s\n", filename_V0_theta);
        return;
    } else{
//printf("opening %s\n", filename_V0_theta);
}
    
    char header[200];
    fgets(header, sizeof(header), inputFile_V0_theta_mrad_solutions); // Skip header

    int i=0; 
    while (fscanf(inputFile_V0_theta_mrad_solutions, "%lf %lf", &list_V0_load_data, &list_theta_mrad_load_data) == 2) {
        list_V0_sols[i] = list_V0_load_data;
        list_theta_mrad_sols[i] = list_theta_mrad_load_data;
        i++;
    }
    Nsolutions = i;
    fclose(inputFile_V0_theta_mrad_solutions);
}

void load_energy_from_EELS_files(int w, double h, double x0_norm_w, int Eelectron_keV, int N, int mode){

    FILE *inputFile_energy;
    char filename_energy[256]; 
    int L = 1000; // nm . convergency parameter for the virtual geometry
  //  double xe=x0_norm_w*w; // nm
    
    
 //   char labelxe[20];
    
 //   if (xe == 0.0) {
 //   sprintf(labelxe, "_xe0");
 //   } else {
 //   sprintf(labelxe, "_xe%.0fW", xe);
 //   }
    
    // Open the energy file to loop over all energies
    sprintf(filename_energy, "datfiles_EELS_along_z/new_files/energy_list_mode%i_N%i_W%.0dnm_h%.0fnm_L%inm_Ee%ikeV_xe%.2f.txt",mode,N,w,h,L,Eelectron_keV,x0_norm_w);    
       
    inputFile_energy=fopen(filename_energy,"r"); 
 if (inputFile_energy == NULL) { // CHANGED: Added file open check
        printf("Error opening %s\n", filename_energy);
        return;
} else{
//printf("opening %s\n", filename_energy);
}

    int i=0;
    double list_energy_load_data;
    while (fscanf(inputFile_energy, "%lf\n", &list_energy_load_data) == 1) {
    if (i >= MAX_ENERGY_POINTS) {
        printf("[ERROR] Too many energy points in %s (max = %d)\n", filename_energy, MAX_ENERGY_POINTS);
        break;
    }
    list_energy_BEM[i] = list_energy_load_data;
    i++;
}

    num_energy_points=i; 
    fclose(inputFile_energy);
 }

int main() {
    double wp_norm_w = 0.5;  // example value
    double d_norm_w = 0.2;   // example value
    double h_norm_w = 0.4;   // example value
    double x0_norm_w = 0.0;  // example value
    int Eelectron_keV = 200; // example value
    int N = 500;             // example value
    double w = 500;             // example value
    double bmin_norm_w = 50.0/w; // example value
    double z_min_val_norm_w = bmin_norm_w + h_norm_w;
    double h=h_norm_w*w ;
    int save_integrand=0; // save z/w vs P(z/w)
    int mode=1; // mode one is the highest peak. more resolution in energy for each mode. see "data_modes.txt" in bem2d in local computer 
    
    char label_txt[200];
    sprintf(label_txt, "_mode%i_Ee%ikeV_wp%.2f_h%.2f_d%.2f_xe%.2f", mode, Eelectron_keV, wp_norm_w, h_norm_w, d_norm_w, x0_norm_w);

    // Load V0 and theta         
    load_solutions_V0_theta_from_file(wp_norm_w, d_norm_w, h_norm_w, x0_norm_w, Eelectron_keV, bmin_norm_w);
    if (Nsolutions == 0) {
    printf("ERROR: No V0/theta solutions loaded.\n");
    return 1;
}
 
    load_energy_from_EELS_files(w, h, x0_norm_w, Eelectron_keV, N, mode);
    if (num_energy_points == 0) {
    printf("ERROR: No energy points loaded.\n");
    return 1;
}
      
    // Loop over all solutions of bmin(V0,theta)
    for (int j = 0; j < Nsolutions; j++) {
        double V0 = list_V0_sols[j];
        double theta_mrad = list_theta_mrad_sols[j];

        char info_bmin[310];
        sprintf(info_bmin, "_V0_%.5feV_theta%.5fmrad_bmin%.2f", V0, theta_mrad, bmin_norm_w);

        char all_info_label[830];
        sprintf(all_info_label, "_w%.0fnm%s%s", w, label_txt, info_bmin);
        
        double list_P_integrated_over_z[num_energy_points];
        double list_energies_P_integrated_over_z[num_energy_points];
        
                                                        
        for (int k = 0; k < num_energy_points; k++) {
            double energy = list_energy_BEM[k];        
            double result = P_integrated_approx(z_min_val_norm_w, energy, w, wp_norm_w, h_norm_w, d_norm_w, x0_norm_w, Eelectron_keV, N, theta_mrad, V0, save_integrand);
                      
            list_P_integrated_over_z[k]=result;
            list_energies_P_integrated_over_z[k]=energy;
            }
        
        // Save the results (energy, Pintegrated over z) for each solution
        char filename_Pintegrated[880];
        sprintf(filename_Pintegrated, "EELS_integrated/P_integrated_over_z_approx%s.dat", all_info_label);
        FILE* output_file = fopen(filename_Pintegrated, "w");
        
        fprintf(output_file, "# Energy(eV)\tP_integrated approx for V0 = %.5f eV, theta = %.5f mrad\n", V0,theta_mrad);
        for (int i = 0; i < num_energy_points; i++) {
        fprintf(output_file, "%.10f %.10f\n", list_energies_P_integrated_over_z[i], list_P_integrated_over_z[i]);
        }

        fclose(output_file); 
        printf("Data saved for j = %d, V0 = %.5f eV, theta_mrad = %.5f\n", j, V0, theta_mrad);
  
    }
        
    return 0;
}

// 1-./run_bmin_vs_theta_V0.sh to run python code, bmin(V0,theta), for different bmin, xe, d/W, h/W + potential data
// 2-./split_positive_eels.sh to separate the EELS files in energies
// 3-g++ -o EELS_integrated_over_z_fixed_mode_approx.out EELS_integrated_over_z_fixed_mode_approx.cpp
// 4-./EELS_integrated_over_z_fixed_mode_approx.out

