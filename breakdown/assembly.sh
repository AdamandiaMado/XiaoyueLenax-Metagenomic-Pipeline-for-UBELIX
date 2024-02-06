#!/usr/bin/env bash
#SBATCH --cpus-per-task=60
#SBATCH --mem=80GB
#SBATCH --time=3-0:00:00
#SBATCH --job-name=assembly
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_assembly_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_assembly_%j.e
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
WORKDIR=""   # Master directory where you want everything to be in there
RAW_DATA_DIR=""   # Raw data directory where the raw fastq files are 
        #WARNING: make sure your data is in _R1.fastq format!

# Replace these with your sample names so it can be analysed automatically in a loop.        
sample_names=("Human1" "Human2" "Human3") # Recommend to manually type in all your sample names to loop through later..
REF_GENOME="/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/UHGG_reps.fasta" 



# Step 2 - Assembly: MetaSPAdes& megahit -------------------------------------------------------------------

    # Warning: This is the most computationally intense stage.
    # Make sure to allocate enough storage!
    for sample_id in "${sample_names[@]}"; do

        metawrap assembly -1 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq -2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq -o ${ASSEMBLY_DIR}/megahit2/${sample_id} -t 60 -m 80
        metawrap assembly --metaspades -1 ${RAW_DATA_DIR}/${sample_id}_R_1.fastq -2 ${RAW_DATA_DIR}//${sample_id}_R_2.fastq -o ${ASSEMBLY_DIR}/metaspades2/${sample_id} -t 60 -m 60 

        # Usage: metaWRAP assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir
        # Options:

        # -1 STR          forward fastq reads
        # -2 STR          reverse fastq reads
        # -o STR          output directory
        # -m INT          memory in GB (default=24)
        # -t INT          number of threads (defualt=1)
        # -l INT          minimum length of assembled contigs (default=1000)

        # --megahit       assemble with megahit (default)
        # --metaspades    assemble with metaspades instead of megahit (better results but slower and higher memory requirement)    
    
    done
    
    # Wait for all parallel tasks to finish
    wait
    
        if [ $? -eq 0 ]; then
        echo "Step 2: assembly completed. Check output at $OUT_DIR/ASSEMBLY"
        else
        echo "Step 2: assembly has failed, please check the input files again or try changing allocated resources."
        exit 1
        fi

# Step 2.5  MetaQUAST -----------------------------------------------------------------------------------------------------------
    
    # Run MetaQUAST to assess quality of the assemblies, independent from metaWRAP
    # With reference
        for sample_id in "${sample_names[@]}"; do
            metaquast.py ${ASSEMBLY_DIR}megahit/$sample_id/final_assembly.fasta --output-dir ${QC_DIR}/QUAST/refs/$sample_id -R ${REF_GENOME} -t 50
            # WIP: copy the metaquast.py to here
        done
        if [ $? -eq 0 ]; then
            echo "Step 2: metaQUAST completed. Check output at ${QC_DIR}/QUAST/refs/$sample_id"
        else
            echo "Step 2: failed, check error report"
            exit 1
        fi

    