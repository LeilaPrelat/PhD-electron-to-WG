#!/bin/bash

# Unified Script: Generate potential data and run Python processing

# --- Settings ---
w=500 ## nm
wp_norm_w=0.5
d_norm_w=0.2
h_norm_w=0.4
Np=400
x0_norm_w=0.0 # Options: 0.0, 0.2, 0.4, 0.6, respect to the width of the middle waveguide
x1_norm_w=$x0_norm_w
epsilon=9 ## DC of saphire for the voltage of the waveguides
b_norm_w=$(echo "scale=8; 50/$w" | bc)
delta=0.1 ## length of the z values for the fitting is 2*delta

# Set xe label
#labelxe="_xe$(printf "%.0f" "$(echo "$x0_norm_w * 100" | bc -l)")W"

# Input and output paths
base_dir="potential_data"
mkdir -p "$base_dir"
 
# Step 1: Check Python script exists
PYTHON_SCRIPT="fit_potential_as_linear.py"
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: $PYTHON_SCRIPT not found in the current directory."
    exit 1
fi

# Step 2: Run Python script with all required parameters
echo "Running Python script with input parameters:"
echo "w=$w, wp_norm_w=$wp_norm_w, d_norm_w=$d_norm_w, h_norm_w=$h_norm_w, Np=$Np, epsilon=$epsilon, x0_norm_w=$x0_norm_w, b_norm_w=$b_norm_w, delta=$delta"

python3 "$PYTHON_SCRIPT" "$w" "$wp_norm_w" "$d_norm_w" "$h_norm_w" "$Np" "$epsilon" "$x0_norm_w" "$b_norm_w" "$delta"

echo "All steps completed."

# End of unified script
#chmod +x run_potential_approx.sh
#./run_potential_approx.sh


