#!/usr/bin/env bash

#extraction
# tar xzvf repository.tar.gz  -C ./datas
cd datas/repository
for first_dir in "$(ls -d */ | cut -f1 -d'/')"; do
    for directory in $first_dir; do
        cd "$directory"
        for second_dir in "$(ls -d */ | cut -f1 -d'/')"; do
            for sub_directory in $second_dir; do
                cd "$sub_directory"
                for file in "$(ls)"; do
                    # Replace tab by commas
                    sed 's/\t/,/g' "$file" > genotype_"$file".csv
                    # Get last column phenotype column
                    cat genotype_"$file".csv | awk -F, '{print $NF}' > phenotype_"$file".csv
                done
                cd ..
            done
        done
        cd ..
    done
done
