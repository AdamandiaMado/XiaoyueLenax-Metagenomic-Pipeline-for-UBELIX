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
