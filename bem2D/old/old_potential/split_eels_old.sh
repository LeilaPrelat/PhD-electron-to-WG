#!/bin/bash

# --- Settings ---
N=500
w=500       # nm
h=200       # nm
ze=50       # nm
Ee=100      # keV
xe_percent=40  # Options: 0, 40, or 80

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
energy_list_file="$output_dir/energy_list_N${N}_a${w}nm_h${h}nm_ze${ze}nm_Ee${Ee}keV${labelxe}.txt"
mkdir -p "$output_dir"
> "$energy_list_file"  # Clear previous list

# Precision for energy values (first column)
precision=7
temp_energy_file=$(mktemp)


file="$base_dir/EELS_along_z_Si_N${N}_a${w}nm_h${h}nm_ze50nm_Ee${Ee}keV${labelxe}.dat"
[ -f "$file" ] || continue

base_filename=$(basename "$file" .dat)  # Keep full name, including zeXXXnm

echo "Processing $file → base name: $base_filename"

awk -v precision="$precision" \
-v output_dir="$output_dir" \
-v base="$base_filename" \
-v temp_file="$temp_energy_file" '
{
if (NF == 0) next;
key = sprintf("%.*f", precision, $1);
groups[key] = groups[key] $0 "\n";
print key >> temp_file;
}
END {
for (k in groups) {
    fname = output_dir "/" base "_energy_" k ".txt";
    print groups[k] > fname;
    close(fname);
}
}' "$file"


# Save sorted unique energies
sort -u "$temp_energy_file" > "$energy_list_file"
rm "$temp_energy_file"

echo "List of unique energies saved to: $energy_list_file"

#chmod +x split_eels2.sh

#./split_eels2.sh

### run the potential and save it in data to use it for the EELS integration (see EELS_integrated_over_z.h) 

## filename_potential="potential_wp${wp_norm_w}_d${d_norm_w}_h${h_norm_w}_N${Np}_xe${x0}nm_z0${z0}nm_z1${z0}nm.txt"
###  ./e.out wp_norm_w d_norm_w h_norm_w N epsilon x0 z0 z1 x1 nx nz > potential_data.dat



