#!/usr/bin/env bash

# Extract raw datas
tar xzvf repository.tar.gz  -C ./datas
# Go in folder containning datas
cd datas/repository
# Loop on first level of directories
for first_dir in "$(ls -d */ | cut -f1 -d'/')"; do
    for directory in $first_dir; do
        cd "$directory"
        # Loop on subdrirectories
        for second_dir in "$(ls -d */ | cut -f1 -d'/')"; do
            for sub_directory in $second_dir; do
                cd "$sub_directory"
                # Loop on files
                for file in *.txt ; do
                    # Get basename
                    basename=${file%.txt}
                    # Replace tab by commas
                    sed 's/\t/,/g' "$file" > temp_"$basename".csv
                    # Get last column phenotype column
                    cat temp_"$basename".csv | awk -F',' '{print $NF}' > phenotype_"$basename".csv
                    # Delete last column in genotype data
                    cat temp_"$basename".csv | awk -F',' '{OFS=","};NF{NF-=1};1' <temp_"$basename".csv >genotype_"$basename".csv
                done
                # Remove raw datas
                rm *.txt
                rm temp*
                cd ..
            done
        done
        cd ..
    done
done
