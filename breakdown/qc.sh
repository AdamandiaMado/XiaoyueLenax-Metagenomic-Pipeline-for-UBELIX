#!/usr/bin/env bash
#SBATCH --cpus-per-task=60
#SBATCH --mem=80GB
#SBATCH --time=3-0:00:00
#SBATCH --job-name=Human3_kraken2_test
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_Human3_kraken2_test_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_Human3_kraken2_test_%j.e
#SBATCH --partition=epyc2
#SBATCH --array=0-6

# ===========================================================================

#           Module setup - Loeading all required packages

#============================================================================
    #module load Anaconda3
    
    #   Load all vital-IT packages preinstalled
        # metawrap is the placeholder name for your metawrap envirnment! 
        conda activate metawrap
    module load vital-it/7
    module load UHTS/Quality_control/fastqc/0.11.9
    module load UHTS/Analysis/trimmomatic/0.36

# ====================================================================================

#                       User Input: Sets working directory

# ====================================================================================
WORKDIR=""   # Master directory where you want everything to be in there, placeholder
RAW_DATA_DIR=""   # Raw data directory where the raw fastq files are , placeholder
        #WARNING: make sure your data is in _R1.fastq format!

# Replace these with your sample names so it can be analysed automatically in a loop.        
sample_names=("Human1" "Human2" "Human3") # Recommend to manually type in all your sample names to loop through later..
REF_GENOME="/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/UHGG_reps.fasta" 



# Step 1 - FastQC--------------------------------------------------------------------------------------------------------------

#       The raw reads must be named with suffix *_R_1.fastq and *_R_2.fastq !!
for sample_id in "${sample_names[@]}"; do
    metawrap read_qc -1 $RAW_DATA_DIR/${sample_id}_R_1.fastq -2 -1 $RAW_DATA_DIR/${sample_id}_R_2.fastq -t 24 -o ${OUT_DIR}
    #Options:

    #    -1 STR          forward fastq reads
    #    -2 STR          reverse fastq reads
    #    -o STR          output directory
    #    -t INT          number of threads (default=1)
    #    -x STR          prefix of host index in bmtagger database folder (default=hg38)

    #    --skip-bmtagger         dont remove human sequences with bmtagger
    #    --skip-trimming         dont trim sequences with trimgalore
    #    --skip-pre-qc-report    dont make FastQC report of input sequences
    #    --skip-post-qc-report   dont make FastQC report of final sequences

done
        # Feedback Check:
        if [ $? -eq 0 ]; then
        echo "Step 1: QC completed successfully, output files at $OUT_DIR/QC"
        else
        echo "Step 1: QC failed. We recommend you restart again."
        exit 1
        fi


