#!/bin/bash
#SBATCH -J glyp-gene-search
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 8
#SBATCH --mem=300G
#SBATCH -t24:00:00
#SBATCH -q inferno
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-14
#SBATCH -o /storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/log/0522_gene-search_out.16
#SBATCH -e /storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/log/0522_gene-search_err.16

source /storage/home/hcoda1/0/skapasi3/r-ktk3-0/miniconda3/etc/profile.d/conda.sh
cd /storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/2_analyses/glyp-gene-search
conda activate binning_env

SAMPLE_LIST="/storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/1_rawdata/cleaned_MAG_files/sample_list_updated.txt"
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
HOMEDIR="/storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/1_rawdata/cleaned_MAG_files"

cat MAGs/*.fasta > all_mags_combined.fasta
makeblastdb -in all_mags_combined.fasta -dbtype nucl -out mags_db

# Search each gene file against MAGs
blastn -query "/storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/2_analyses/glyp-gene-search/refs/yhhS.fna" -db mags_db -out yhhS_vs_mags.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -perc_identity 90
blastn -query "/storage/home/hcoda1/0/skapasi3/r-ktk3-0/HIAMI/2_analyses/glyp-gene-search/refs/epsps.fna" -db mags_db -out epsps_vs_mags.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -perc_identity 90

mkdir -p blast_results

for mag in MAGs/*.fasta; do
    mag_name=$(basename "$mag" .fasta)
    # Build per-MAG database
    makeblastdb -in "$mag" -dbtype nucl -out "blast_results/${mag_name}_db"
    # Search both gene files
    for gene in refs/yhhS.fna refs/epsps.fna; do
        gene_name=$(basename "$gene" .fna)
        blastn -query "$gene" -db "blast_results/${mag_name}_db" -out "blast_results/${gene_name}_vs_${mag_name}.tsv" -outfmt 6 -evalue 1e-5 -perc_identity 90
    done
done

cat refs/yhhS.fna refs/epsps.fna > genes_combined.fasta
bowtie2-build genes_combined.fasta genes_index

R1="${HOMEDIR}/${SAMPLE_NAME}/R1_lcombined.fq.gz"
R2="${HOMEDIR}/${SAMPLE_NAME}/R2_lcombined.fq.gz"

bowtie2 -x genes_index -1 "$R1" -2 "$R2" -S "${SAMPLE_NAME}_vs_genes.sam" --no-unal -p 8 2> "${SAMPLE_NAME}_mapping_stats.txt"
# Convert to sorted BAM
samtools sort "${SAMPLE_NAME}_vs_genes.sam" -o "${SAMPLE_NAME}_vs_genes.bam"
samtools index "${SAMPLE_NAME}_vs_genes.bam"
        
# Get coverage stats per gene
samtools coverage "${SAMPLE_NAME}_vs_genes.bam" > "${SAMPLE_NAME}_coverage.tsv"
echo "Done: $SAMPLE_NAME"
