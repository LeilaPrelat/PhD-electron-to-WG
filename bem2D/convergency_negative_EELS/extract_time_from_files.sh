#!/bin/bash

# list of (L, N) pairs
list_LN=(
  "700 500"
  "1000 800"
  "1000 1000"
  "1500 800"
  "1500 1000"
  "1000 1500"
  "1500 1500"
  "2000 1000"
)

# output file
output="results.txt"
> "$output"   # clear file

# add header
echo -e "L(nm)\tN\tTime(min)" >> "$output"

for pair in "${list_LN[@]}"; do
    L=$(echo $pair | awk '{print $1}')
    N=$(echo $pair | awk '{print $2}')

    pattern="EELS_comsol2_N${N}_W500nm_h200nm_L${L}nm_Ee200keV_energy035.sh.o*"
    files=($pattern)

    if [[ -f "${files[0]}" ]]; then
        # grab minutes (number before 'm')
        value=$(grep "^real" "${files[0]}" | awk '{print $2}' | sed 's/m.*//')
        echo -e "$L\t$N\t$value" >> "$output"
    else
        echo -e "$L\t$N\tFILE_NOT_FOUND" >> "$output"
    fi
done

