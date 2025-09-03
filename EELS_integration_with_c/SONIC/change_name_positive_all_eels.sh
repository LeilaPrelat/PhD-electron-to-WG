#!/bin/bash
# Here we only save the energy if ALL the EELS is positive for different z

# --- Settings ---
N=500
w=500       # nm
h=200       # nm
L=1000      # nm. convergency parameter for the virtual geometry
Ee=200      # keV
xe_norm_w=0.00  # normalized to w: options are 0.00,0.40,0.80
mode=3 ## mode 1 is the highest. more resolution in energy for each mode. see "data_modes.txt" in bem2d in local computer 

# Energy ranges depending on mode
if [[ $mode -eq 1 ]]; then
    list_energy=($(seq 0.7 0.001 0.84))
elif [[ $mode -eq 2 ]]; then
    list_energy=($(seq 0.28 0.001 0.44))
elif [[ $mode -eq 3 ]]; then
    list_energy=($(seq 1.45 0.001 1.61))
fi

# Input and output paths
base_dir="datfiles_EELS_along_z"
output_dir="${base_dir}/new_files"
mkdir -p "$output_dir"
 
# Output file for energy list
energy_list_file="$output_dir/energy_list_mode${mode}_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_xe${xe_norm_w}.txt"
> "$energy_list_file"  # Clear any old content

# --- Main loop ---
for energy in "${list_energy[@]}"; do
    # Convert energy to meV for input file naming
    energy_meV=$(printf "%04d" "$(echo "$energy * 1000" | bc -l | awk '{printf "%d", $1}')")
    
    # Construct input filename
    input_filename="EELS_comsol2_vs_z_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_energy${energy_meV}.dat"
    input_file="${base_dir}/${input_filename}"

    # Check if input file exists
    if [[ -f "$input_file" ]]; then
        # Check if ALL values in 5th column are positive
        if awk 'NF && $5 <= 0 {exit 1}' "$input_file"; then
            formatted_energy=$(printf "%.4f" "$energy")
            output_filename="EELS_along_z_N${N}_W${w}nm_h${h}nm_L${L}nm_Ee${Ee}keV_xe${xe_norm_w}_energy${formatted_energy}.dat"
            output_file="${output_dir}/${output_filename}"

            # Copy file as-is (all values pass)
            cp "$input_file" "$output_file"

            # Append to energy list
            echo "$formatted_energy" >> "$energy_list_file"
        fi
    else
        echo "Warning: File not found: $input_file"
    fi
done

echo "Processing complete. Files with all-positive 5th column saved to $output_dir"
echo "Energy list: $energy_list_file"

#chmod +x change_name_positive_all_eels.sh

#./change_name_positive_all_eels.sh


