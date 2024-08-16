#!/bin/bash

# Define the base folder
base_folder="/gpfs/home/rs1521"

# Define folders A and B
folders=("T_kodakarensis" "P_furiosus")

# Iterate through each folder A and B
for folder in "${folders[@]}"; do
    # Read each line (ID) in the top_hit_ids.txt file
    while IFS= read -r id; do
        # Extract the part before the last "__"
        bait_id="${id%__*}"

        # Define the path to the source folder C and destination folder D
        predictions_folder="${base_folder}/${folder}/${bait_id}_predictions/"
        predictions_folder_alt="${base_folder}/${folder}/predictions_${bait_id}_cf155/"
        top_structures_folder="${base_folder}/${folder}/top_structures/"
        
        for file in "${predictions_folder}${id}"_*.pdb; do
            # If the file exists, copy it to folder D
            if [[ -e "$file" ]]; then
                cp "$file" "$top_structures_folder"
                echo "Copied $file to $top_structures_folder"
            fi
        done

        for file in "${predictions_folder_alt}${id}"_*.pdb; do
            # If the file exists, copy it to folder D
            if [[ -e "$file" ]]; then
                cp "$file" "$top_structures_folder"
                echo "Copied $file to $top_structures_folder"
            fi
        done

    done < "${base_folder}/${folder}/top_hit_ids.txt"
done
