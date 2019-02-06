#!/usr/bin/env bash
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
#
# This script must be run next to the repository.tar.gz archive that collects gamete data.
# It will then decompress the archive and reformat data so that they can be used directly by our programs.
mkdir datas
# Extract raw datas
tar xzvf repository.tar.gz  -C ./datas
# Go in folder containning datas
cd datas/repository
echo -e '\033[33mParse files to be suitable with input of methods \033[0m'
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
