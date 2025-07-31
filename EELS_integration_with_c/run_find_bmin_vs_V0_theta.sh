#!/bin/bash

# Unified Script: Generate potential data and run Python processing

# --- Settings ---
w=400 ## nm
wp_norm_w=0.5
d_norm_w=0.2
h_norm_w=0.5
Np=400
x0_norm_w=0 # Options: 0.0, 0.2, 0.4, 0.6, respect to the width of the middle waveguide
x1_norm_w=$x0_norm_w
z0_norm_w=$(echo "$h_norm_w * 1.001" | bc) # z0 close to h (z=0 at the interface between waveguide and substrate)
z1_norm_w=$(echo "$z0_norm_w + 400 / $w" | bc -l) # z1 = z0 + 400 / w 
nz=200
nx=1
epsilon=9 ## DC of saphire for the voltage of the waveguides
Ee_electron_keV=200
bmin_norm_w=$(echo "scale=8; 50/$w" | bc)

# Set xe label
#labelxe="_xe$(printf "%.0f" "$(echo "$x0_norm_w * 100" | bc -l)")W"

# Input and output paths
base_dir="potential_data"
mkdir -p "$base_dir"

# Build output filename
filename_potential="potential_wp${wp_norm_w}_d${d_norm_w}_h${h_norm_w}_N${Np}_xe${x0_norm_w}.txt"

# Step 1: Get normalization factor (absolute value of first V)
normalization_factor=$(./e.out "$wp_norm_w" "$d_norm_w" "$h_norm_w" "$Np" "$epsilon" "0" "0" "1" "$z0_norm_w" "$z0_norm_w" "2" | \
awk 'NR==1 { print ($2 < 0) ? -$2 : $2 }')

echo "Normalization factor: $normalization_factor"

# Step 2: Run e.out to create normalized potential data
./e.out "$wp_norm_w" "$d_norm_w" "$h_norm_w" "$Np" "$epsilon" "$x0_norm_w" "$x1_norm_w" "$nx" "$z0_norm_w" "$z1_norm_w" "$nz" | \
awk -v factor="$normalization_factor" 'BEGIN { print "# z/w V/V0" } { printf "%s %.10f\n", $1, $2 / factor }' > "${base_dir}/${filename_potential}"

echo "Potential data saved to ${base_dir}/${filename_potential}"

# Step 3: Check Python script exists
PYTHON_SCRIPT="find_bmin_vs_V0_theta.py"
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: $PYTHON_SCRIPT not found in the current directory."
    exit 1
fi

# Step 4: Run Python script with all required parameters
echo "Running Python script with input parameters:"
echo "w=$w, wp_norm_w=$wp_norm_w, d_norm_w=$d_norm_w, h_norm_w=$h_norm_w, Np=$Np, epsilon=$epsilon, x0_norm_w=$x0_norm_w, Ee_electron_keV=$Ee_electron_keV, bmin_norm_w=$bmin_norm_w"

python3 "$PYTHON_SCRIPT" "$w" "$wp_norm_w" "$d_norm_w" "$h_norm_w" "$Np" "$epsilon" "$x0_norm_w" "$Ee_electron_keV" "$bmin_norm_w"

echo "All steps completed."

# End of unified script
#chmod +x run_find_bmin_vs_V0_theta.sh
#./run_find_bmin_vs_V0_theta.sh


