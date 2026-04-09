#!/bin/bash
#SBATCH -J drep-biologic
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t24:00:00
#SBATCH -q inferno
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-13%5
#SBATCH -o /storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/log/040926_drep-biologic_out.16
#SBATCH -e /storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/log/040926_drep-biologic_err.16

source /storage/home/hcoda1/0/skapasi3/miniconda3/etc/profile.d/conda.sh
cd /storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/1_rawdata/cleaned_MAG_files/
conda activate drep_env

PARENT_DIR="/storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/1_rawdata/cleaned_MAG_files"
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_list_updated.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_ID}"

# biological dRep - 99%

BIO_DREP_DIR="${PARENT_DIR}/biological_drep"
mkdir -p "${BIO_DREP_DIR}/pooled_99" \
         "${BIO_DREP_DIR}/combined_genome_info"
# Initialize pooled genome info file
echo "genome,completeness,contamination" \
    > "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_99.csv"

ALL_BINS=""

while IFS= read -r SAMPLE_ID; do
    SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_ID}"
    # Collect non-norm dereplicated bins
    NONNORM_DEREP="${SAMPLE_DIR}/1_drep_output/dereplicated_genomes"
    if [ -d "$NONNORM_DEREP" ]; then
        ALL_BINS="$ALL_BINS $(find "$NONNORM_DEREP" -name "*.fasta")"
        tail -n +2 "${SAMPLE_DIR}/genome_info.csv" \
            >> "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_99.csv"
    else
        echo "WARNING: Missing $NONNORM_DEREP" >&2
    fi

    # Collect norm dereplicated bins
    NORM_DEREP="${SAMPLE_DIR}/1_norm_drep_output/dereplicated_genomes"
    if [ -d "$NORM_DEREP" ]; then
        ALL_BINS="$ALL_BINS $(find "$NORM_DEREP" -name "*.fasta")"
        tail -n +2 "${SAMPLE_DIR}/norm_genome_info.csv" \
            >> "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_99.csv"
    else
        echo "WARNING: Missing $NORM_DEREP" >&2
    fi
done < sample_list_updated.txt
if [ -z "$ALL_BINS" ]; then
    echo "ERROR: No dereplicated bins found" && exit 1
fi

echo "Running biological dRep (99% ANI) on pooled bins..."
dRep dereplicate "${BIO_DREP_DIR}/pooled_99/" -g ${ALL_BINS} --genomeInfo "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_99.csv" -comp 50 -con 10 -pa 0.9 -sa 0.99 -nc 0.30 -p 32
echo "Done: $SAMPLE_ID biological dRep at 99% ANI"

# biological dRep - 95%

mkdir -p "${BIO_DREP_DIR}/pooled_95"
echo "genome,completeness,contamination" \
    > "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_95.csv"

ALL_BINS=""

while IFS= read -r SAMPLE_ID; do
    SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_ID}"
    # Collect non-norm dereplicated bins
    NONNORM_DEREP="${SAMPLE_DIR}/1_drep_output/dereplicated_genomes"
    if [ -d "$NONNORM_DEREP" ]; then
        ALL_BINS="$ALL_BINS $(find "$NONNORM_DEREP" -name "*.fasta")"
        tail -n +2 "${SAMPLE_DIR}/genome_info.csv" \
            >> "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_95.csv"
    else
        echo "WARNING: Missing $NONNORM_DEREP" >&2
    fi
    # Collect norm dereplicated bins
    NORM_DEREP="${SAMPLE_DIR}/1_norm_drep_output/dereplicated_genomes"
    if [ -d "$NORM_DEREP" ]; then
        ALL_BINS="$ALL_BINS $(find "$NORM_DEREP" -name "*.fasta")"
        tail -n +2 "${SAMPLE_DIR}/norm_genome_info.csv" \
            >> "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_95.csv"
    else
        echo "WARNING: Missing $NORM_DEREP" >&2
    fi
done < sample_list_updated.txt

if [ -z "$ALL_BINS" ]; then
    echo "ERROR: No dereplicated bins found" && exit 1
fi

echo "Running biological dRep (95% ANI) on pooled bins..."
dRep dereplicate "${BIO_DREP_DIR}/pooled_95/" -g ${ALL_BINS} --genomeInfo "${BIO_DREP_DIR}/combined_genome_info/pooled_genome_info_95.csv" -comp 50 -con 10 -pa 0.9 -sa 0.95 -nc 0.30 -p 32
echo "Done: $SAMPLE_ID biological dRep at 95% ANI"
