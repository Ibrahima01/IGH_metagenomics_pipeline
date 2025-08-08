#!/bin/bash

# Dossier contenant les fichiers FASTQ
fastq_dir="trimmed_reads"

# Dossier de sortie pour les assemblages SPAdes
output_dir="SPAdes_Assemblies"
mkdir -p "${output_dir}"

# Chemin vers l'exécutable de SPAdes
spades_exec="/home/ibrahima/Metagenomic/SPAdes-4.2.0-Linux/bin/spades.py"

# Boucle sur tous les fichiers paired-end R1
for R1_file in ${fastq_dir}/*_paired_R1.fastq.gz; do
    # Extraire le nom du sample (avant le premier "_paired")
    sample_name=$(basename "${R1_file}" | sed 's/_paired_R1.*//')

    # Définir le fichier R2 correspondant
    R2_file="${fastq_dir}/${sample_name}_paired_R2.fastq.gz"

    # Dossier de sortie pour ce sample
    sample_output_dir="${output_dir}/${sample_name}"

    # Sauter si l'assemblage existe déjà
    if [[ -f "${sample_output_dir}/contigs.fasta" ]]; then
        echo "✅ Skipping ${sample_name} (assembly already exists)"
        continue
    fi

    # Vérifier que les deux fichiers existent
    if [[ -f "${R1_file}" && -f "${R2_file}" ]]; then
        echo "🔄 Processing sample: ${sample_name}"

        mkdir -p "${sample_output_dir}"

        # Lancer SPAdes
        python3 "${spades_exec}" \
            -1 "${R1_file}" \
            -2 "${R2_file}" \
            -o "${sample_output_dir}" \
            --threads 16 --memory 128

        echo "✅ Assembly completed for sample: ${sample_name}"

    else
        echo "❌ Paired files not found for sample: ${sample_name}"
    fi
done

echo "🎉 SPAdes pipeline completed for all samples."

