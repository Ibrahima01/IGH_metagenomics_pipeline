#!/bin/bash

# === PATHS ===
fastq_dir="trimmed_reads"
reference_dir="genomes"
output_base_dir="Assemblies"
csv_path="virus_sample_10_contigs_counts.csv"

# === PRÉPARATION ===
mkdir -p "${output_base_dir}"

# === EXTRAIRE LA MAPPING VIRUS → SAMPLES ===
declare -A virus_to_samples

# === EXTRAIRE LA MAPPING VIRUS → SAMPLES ===
declare -A virus_to_samples

# Lecture correcte dans le shell parent
while IFS=',' read -r virus taxid sample count; do
    [[ "$virus" == "Virus" ]] && continue  # ignorer l'en-tête
    virus_key=$(echo "$virus" | tr ' ' '_')
    virus_to_samples["$virus_key"]+="${sample} "
done < "${csv_path}"

# === BOUCLE SUR LES GENOMES DE RÉFÉRENCE ===
for reference_genome in "${reference_dir}"/*.fasta; do
    ref_name=$(basename "${reference_genome}" .fasta)

    # Vérifier si des samples sont associés à ce virus
    if [[ -z "${virus_to_samples[$ref_name]}" ]]; then
        echo "⏭️  Aucun sample à traiter pour le virus ${ref_name}. Skip..."
        continue
    fi

    echo "🦠 Traitement du virus : ${ref_name}"

    output_dir="${output_base_dir}/${ref_name}"
    mkdir -p "${output_dir}"

    # Indexation BWA si nécessaire
    if [[ ! -f "${reference_genome}.bwt" ]]; then
        echo "🔧 Indexation BWA de ${reference_genome}..."
        bwa index "${reference_genome}"
    fi

    # Boucler sur les samples associés à ce virus
    for sample_name in ${virus_to_samples[$ref_name]}; do
        R1_file="${fastq_dir}/${sample_name}_paired_R1.fastq.gz"
        R2_file="${fastq_dir}/${sample_name}_paired_R2.fastq.gz"
        consensus_file="${output_dir}/${sample_name}.fa"

        # Skip si déjà fait
        if [[ -f "${consensus_file}" ]]; then
            echo "✅ ${sample_name} déjà assemblé pour ${ref_name}. Skip..."
            continue
        fi

        if [[ -f "${R1_file}" && -f "${R2_file}" ]]; then
            echo "🔄 Assemblage de ${sample_name} sur ${ref_name}..."

            output_bam="${output_dir}/${sample_name}.bam"

            bwa mem -t 50 "${reference_genome}" "${R1_file}" "${R2_file}" | \
                samtools view -u -@ 3 - | \
                samtools sort -@ 3 -o "${output_bam}"

            echo "🧬 Mapping terminé pour ${sample_name}"

            samtools mpileup -A -d 0 -Q 0 -B "${output_bam}" | \
                ivar consensus -p "${output_dir}/${sample_name}" -t 0.5 -m 10

            echo "✅ Consensus généré pour ${sample_name} → ${ref_name}"
        else
            echo "❌ Fichiers manquants pour ${sample_name}"
        fi
    done
done

echo "🎉 Assemblage terminé pour tous les samples sélectionnés."
