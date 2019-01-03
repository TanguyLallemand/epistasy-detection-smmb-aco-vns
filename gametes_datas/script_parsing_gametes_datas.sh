#!/usr/bin/env bash

#extraction
tar xzvf repository.tar.gz  -C ./datas
cd datas/repository
for first_dir in "$(ls -d */ | cut -f1 -d'/')"; do
    for directory in $first_dir; do
        cd "$directory"
        for second_dir in "$(ls -d */ | cut -f1 -d'/')"; do
            for sub_directory in $second_dir; do
                cd "$sub_directory"
                for file in *.txt ; do
                    # Get basename
                    basename=${file%.txt}
                    # Replace tab by commas
                    sed 's/\t/,/g' "$file" > temp_"$basename".csv
                    # Get last column phenotype column
                    cat temp_"$basename".csv | awk -F',' '{print $NF}' > phenotype_"$basename".csv
                    cat temp_"$basename".csv | awk -F',' '{OFS=","};NF{NF-=1};1' <temp_"$basename".csv >genotype_"$basename".csv
                done
                rm *.txt
                cd ..
            done
        done
        cd ..
    done
done
