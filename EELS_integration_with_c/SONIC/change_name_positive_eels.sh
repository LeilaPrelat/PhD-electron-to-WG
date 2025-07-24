#!/bin/bash

# --- Settings ---
N=500
w=400       # nm
h=200       # nm
L=500       # nm. convergency parameter for the virtual geometry
Ee=200      # keV
xe_norm_w=0.00  # normalized to w: options are 0.00,0.40,0.80

list_energy=($(seq 0.2 0.005 3)) ## eV

# Input and output paths
base_dir="datfiles_EELS"
output_dir="${base_dir}/new_files"
mkdir -p "$output_dir"
 
# Construct input filename
filename="EELS_along_z_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_energy${energy}.dat"
input_file="${base_dir}/${filename}"

# Output file for energy list
energy_list_file="$output_dir/energy_list_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_xe${xe_norm_w}.txt"
> "$energy_list_file"  # Clear any old content

# --- Main loop ---
for energy in "${list_energy[@]}"; do
    # Convert energy to meV for input file naming
    energy_meV=$(printf "%04d" "$(echo "$energy * 1000" | bc -l | awk '{printf "%d", $1}')")
    
    # Construct input filename
    input_filename="EELS_along_z_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_energy${energy_meV}.dat"
    input_file="${base_dir}/${input_filename}"

    # Check if input file exists
    if [[ -f "$input_file" ]]; then
        # Format energy to 3 decimal places for output filename
        formatted_energy=$(printf "%.3f" "$energy")
        
        # Construct output filename
        output_filename="EELS_along_z_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_xe${xe_norm_w}_energy${formatted_energy}.dat"
        output_file="${output_dir}/${output_filename}"

        # Filter rows where the 5th column is positive and write to output
        awk '$5 > 0' "$input_file" > "$output_file"

        # Append the energy to the energy list file
        # If filtered file is not empty, add energy to list
        if [[ -s "$output_file" ]]; then
            echo "$formatted_energy" >> "$energy_list_file"
        else
            # Optionally remove empty filtered file
            rm -f "$output_file"
        fi
    else
        echo "Warning: File not found: $input_file"
    fi
done

echo "Processing complete. Filtered files saved to $output_dir"
echo "Energy list (only positive EELS): $energy_list_file"

#chmod +x change_name_positive_eels.sh

#./change_name_positive_eels.sh


