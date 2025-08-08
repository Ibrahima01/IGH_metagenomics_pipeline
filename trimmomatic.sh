#!/bin/bash

# Define the directory containing all the FASTQ files
base_dir="Samples"

# Define the output directory for cleaned files
output_dir="trimmed_reads"
mkdir -p ${output_dir}

# Define the path to the Kraken database
kraken_db="/home/ibrahima/BasespaceCLI/EID_LASV_Zoonotic_27012023/Data/Intensities/BaseCalls/kraken_db_virus"

# Define the path to the Kraken executable
kraken_executable="/home/ibrahima/anaconda3/bin/kraken2"

# Path to the Trimmomatic adapters file (adjust this path as necessary)
adapters="/home/ibrahima/anaconda3/share/trimmomatic/adapters/TruSeq3-PE-2.fa"

# Loop through all the FASTQ files in the directory
for R1 in ${base_dir}/*_R1*.fastq.gz; do
    if [ -f "${R1}" ]; then
        # Define the corresponding R2 file
        R2=$(echo ${R1} | sed 's/_R1/_R2/')

        # Extract the sample name before "_R1" or "_R2" including everything before those patterns
        sample_name=$(basename ${R1} | sed 's/_R1.*//')

        # Check if both R1 and R2 files exist
        if [[ -f "${R1}" && -f "${R2}" ]]; then
            echo "Processing sample: ${sample_name}"

            # Define the output file names for paired and unpaired reads
            output_paired_R1="${output_dir}/${sample_name}_paired_R1.fastq.gz"
            output_unpaired_R1="${output_dir}/${sample_name}_unpaired_R1.fastq.gz"
            output_paired_R2="${output_dir}/${sample_name}_paired_R2.fastq.gz"
            output_unpaired_R2="${output_dir}/${sample_name}_unpaired_R2.fastq.gz"

            # Run Trimmomatic
            trimmomatic PE -phred33 \
                ${R1} ${R2} \
                ${output_paired_R1} ${output_unpaired_R1} \
                ${output_paired_R2} ${output_unpaired_R2} \
                ILLUMINACLIP:${adapters}:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

            # Define the output file names for Kraken
            kraken_output="${output_dir}/${sample_name}_kraken_output.txt"
            kraken_report="${output_dir}/${sample_name}_kraken_report.txt"

            # Run Kraken
            ${kraken_executable} --db ${kraken_db} --paired ${output_paired_R1} ${output_paired_R2} --output ${kraken_output} --report ${kraken_report}
        else
            echo "Paired files not found for sample: ${sample_name}"
        fi
    fi
done

echo "Trimmomatic processing completed."
