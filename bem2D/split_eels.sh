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

# Temporary file for collecting energy values
temp_energies=$(mktemp)

# Split based on unique energy values (1st column, rounded to 8 decimals)
awk -v outdir="$output_dir" -v base="${filename%.dat}" -v temp="$temp_energies" '
{
    energy = sprintf("%.7f", $1)
    print energy >> temp
    data[energy] = data[energy] $0 "\n"
}
END {
    for (e in data) {
        file = outdir "/" base "_energy_" e ".txt"
        printf "%s", data[e] > file
    }
}
' "$input_file"

# Sort and remove duplicates to get unique energy list
sort -u "$temp_energies" > "$energy_list_file"
rm "$temp_energies"

echo "EELS file split into per-energy files in: $output_dir"
echo "Energy list saved to: $energy_list_file"

#chmod +x split_eels.sh

#./split_eels.sh


