#!/bin/bash

# Répertoire contenant les assemblages SPAdes
assembly_dir="SPAdes_Assemblies"

# Répertoire de sortie pour les résultats Kraken2
kraken_output_dir="Kraken_Contigs_005_confidence"
mkdir -p "${kraken_output_dir}"

# Base de données Kraken2
kraken_db="/home/ibrahima/kraken_db_standard"

# Chemin vers Kraken2
kraken_executable="/home/ibrahima/anaconda3/bin/kraken2"

# Boucle sur tous les fichiers contigs.fasta
for contigs_file in ${assembly_dir}/*/contigs.fasta; do
    # Extraire le nom du sample depuis le chemin
    sample_name=$(basename $(dirname "${contigs_file}"))

    # Définir les fichiers de sortie
    kraken_output="${kraken_output_dir}/${sample_name}_kraken_contigs_output.txt"
    kraken_report="${kraken_output_dir}/${sample_name}_kraken_contigs_report.txt"

    # Sauter si le résultat existe déjà
    if [[ -f "${kraken_report}" ]]; then
        echo "✅ Skipping ${sample_name} (already classified)"
        continue
    fi

    echo "🔬 Classifying contigs for sample: ${sample_name}"

    # Exécuter Kraken2
    ${kraken_executable} --db ${kraken_db} \
                         --output ${kraken_output} \
                         --report ${kraken_report} \
                         --use-names \
                         --threads 50 \
			 --confidence 0.05\
                         "${contigs_file}"

    echo "✅ Classification completed for sample: ${sample_name}"
done

echo "🎯 Kraken2 classification done for all contigs."
