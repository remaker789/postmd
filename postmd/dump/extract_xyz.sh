#!/bin/bash

# Extract specific elements from XYZ file while preserving XYZ format
# For Linux systems
# Usage: ./extract_xyz.sh -i input.xyz -e "Na K" [-o output.xyz]

# Initialize variables
input_file=""
output_file=""
elements=""

# Parse command line arguments
while getopts "i:o:e:" opt; do
    case $opt in
        i) input_file="$OPTARG";;
        o) output_file="$OPTARG";;
        e) elements="$OPTARG";;
        \?) echo "Invalid option: -$OPTARG" >&2
            exit 1;;
        :) echo "Option -$OPTARG requires an argument." >&2
            exit 1;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_file" ] || [ -z "$elements" ]; then
    echo "Usage: $0 -i input.xyz -e \"element1 element2 ...\" [-o output.xyz]"
    echo "Options:"
    echo "  -i: Input XYZ file"
    echo "  -e: Space-separated list of elements to extract"
    echo "  -o: Output XYZ file (optional)"
    exit 1
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Generate default output filename if not provided
if [ -z "$output_file" ]; then
    base_name="${input_file%.*}"
    elements_str=$(echo "$elements" | tr ' ' '_')
    output_file="${base_name}_${elements_str}.xyz"
fi

# Create AWK script for XYZ format extraction
cat > extract_xyz.awk << EOF
BEGIN {
    atoms_count = 0
    frame_count = 0
    in_frame = 0
    buffer = ""
}

# Read number of atoms (first line of each frame)
/^[0-9]+$/ {
    if (frame_count > 0) {
        # Print previous frame if we found matching atoms
        if (atoms_count > 0) {
            print atoms_count
            print buffer
        }
    }
    total_atoms = \$1
    in_frame = 1
    atoms_count = 0
    buffer = ""
    frame_count++
    getline  # Skip comment line
    buffer = buffer \$0 "\\n"
    next
}

# Process atom lines
{
    if (in_frame) {
        split(\$0, fields)
        element = fields[1]
        for (e in elements_array) {
            if (element == e) {
                atoms_count++
                buffer = buffer \$0 "\\n"
                break
            }
        }
    }
}

END {
    # Print last frame if we found matching atoms
    if (atoms_count > 0) {
        print atoms_count
        print buffer
    }
}
EOF

# Convert space-separated elements to AWK array initialization
elements_init=$(echo "$elements" | awk '{
    printf "elements_array["
    for(i=1; i<=NF; i++) {
        if(i>1) printf ","
        printf "\"%s\"",$i
    }
    printf "]=\"\""
}')

# Run AWK script with elements array
awk -v "$elements_init" -f extract_xyz.awk "$input_file" > "$output_file"

# Clean up temporary file
rm extract_xyz.awk

echo "Extracted XYZ file has been saved to $output_file"