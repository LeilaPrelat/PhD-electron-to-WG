#!/bin/bash

# Define fixed parameters
N=700
W=500
h=300
L=500
Ee=200

list_energy=($(seq 0.2 0.005 3))

# Output file
output="datfiles_EELS/EELS_comsol2_spectrum_N${N}_W${W}nm_h${h}nm_L${L}nm_Ee${Ee}keV.dat"
> "$output"  # Clear or create the file

for energy in "${list_energy[@]}"; do
    # Format energy for filename: e.g. 0.200 -> 0200
    energy_formatted=$(printf "%04d" "$(echo "$energy * 1000" | bc | cut -d'.' -f1)")

    filename="datfiles_EELS/EELS_comsol2_N${N}_W${W}nm_h${h}nm_L${L}nm_Ee${Ee}keV_energy${energy_formatted}.dat"

    if [[ -f "$filename" ]]; then
       eels=$(awk 'NR==1 {print $5}' "$filename")
    if [[ -n "$eels" ]]; then
        printf "%.8f\t%s\n" "$energy" "$eels" >> "$output"
    else
        echo "File exists but EELS value is empty: $filename" >&2
    fi
else
    echo "Missing file: $filename" >&2
fi

done

#chmod +x unify_dat_files_in_energy.sh

#./unify_dat_files_in_energy.sh
