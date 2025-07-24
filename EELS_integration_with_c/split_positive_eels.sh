#!/bin/bash

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

# Check input file
if [ ! -f "$input_file" ]; then
    echo "File not found: $input_file"
    exit 1
fi

# Temporary file for energy-EELS chunks
temp_prefix=$(mktemp -u)
awk -v outdir="$output_dir" -v base="${filename%.dat}" -v temp="$temp_prefix" '
{
    energy = sprintf("%.7f", $1)
    data[energy] = data[energy] $0 "\n"
}
END {
    for (e in data) {
        file = temp "_energy_" e ".tmp"
        printf "%s", data[e] > file
    }
}
' "$input_file"

# Process each temporary file
for temp_file in "${temp_prefix}"_energy_*.tmp; do
    energy=$(basename "$temp_file" | sed -E "s/.*_energy_(.*)\.tmp/\1/")
    
    # Check that *all* EELS values (5th column) are strictly positive
    if awk '{if ($5 <= 0) exit 1}' "$temp_file"; then
        final_file="${output_dir}/${filename%.dat}_energy_${energy}.txt"
        mv "$temp_file" "$final_file"
        echo "$energy" >> "$energy_list_file"
    else
        rm "$temp_file"
    fi
done

echo "Filtered EELS files saved to: $output_dir"
echo "Energy list (only positive EELS): $energy_list_file"

#chmod +x split_eels.sh

#./split_eels.sh


