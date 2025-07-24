#!/bin/bash

# [EXTRA. NO NEED]
# --- Settings ---
N=500
w=500       # nm
h=300       # nm
ze=50       # nm
Ee=200      # keV
xe_percent=0  # Options: 0, 40, or 80

# Set xe label
if [ "$xe_percent" -eq 40 ]; then
    labelxe="_xe40W"
elif [ "$xe_percent" -eq 80 ]; then
    labelxe="_xe80W"
else
    labelxe="_xe0"
fi

# Input and output paths
base_dir="datfiles_EELS_along_z"
output_dir="${base_dir}/new_files"
mkdir -p "$output_dir"
 
# Construct input filename
filename="EELS_along_z_Si_N${N}_a${w}nm_h${h}nm_ze${ze}nm_Ee${Ee}keV${labelxe}.dat"
input_file="${base_dir}/${filename}"

# Output file for energy list
energy_list_file="$output_dir/energy_list_N${N}_a${w}nm_h${h}nm_ze${ze}nm_Ee${Ee}keV${labelxe}.txt"
> "$energy_list_file"  # Clear any old content

#Check if energy list exists and is not empty
if [[ ! -s "$energy_list_file" ]]; then
    echo "❌ Energy list file not found or is empty: $energy_list_file"
    exit 1
fi

while read -r energy; do
    [[ -z "$energy" ]] && continue  # skip empty lines

    file="${output_dir}/${filename%.dat}_energy_${energy}.txt"

    if [[ -f "$file" ]]; then
        echo "Checking $file"
        if awk '{ if ($5 <= 0) { print "❌ Non-positive EELS in: " FILENAME; exit 1 } }' "$file"; then
            echo "✅ All EELS values positive in: $file"
        fi
    else
        echo "⚠️ File missing for energy: $energy"
    fi
done < "$energy_list_file"


#chmod +x check_EELS_positive.sh

#./check_EELS_positive.sh


