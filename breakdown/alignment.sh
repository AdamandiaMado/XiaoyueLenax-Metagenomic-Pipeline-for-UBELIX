#!/usr/bin/env bash
#SBATCH --cpus-per-task=60
#SBATCH --mem=80GB
#SBATCH --time=3-0:00:00
#SBATCH --job-name=alignment
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_alignment_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_alignment_%j.e
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
    module load Blast/blast/2.2.26
    module load Blast/ncbi-blast/2.9.0+
    module load UHTS/Aligner/bowtie2/2.3.4.1
    module load UHTS/Assembler/SPAdes/3.15.4
    module load UHTS/Assembler/megahit/1.1.4
    module load UHTS/Quality_control/quast/4.6.0
    module load UHTS/Analysis/samtools/1.10
    module load UHTS/Analysis/metabat/2.12.1

# ====================================================================================

#                       User Input: Sets working directory

# ====================================================================================
    WORKDIR=""   # Master directory where you want everything to be in there
    RAW_DATA_DIR=""   # Raw data directory where the raw fastq files are 
            #WARNING: make sure your data is in _R1.fastq format!

    # Replace these with your sample names so it can be analysed automatically in a loop.        
    sample_names=("Human1" "Human2" "Human3") # placeholder, replace later
    REF_GENOME="/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/UHGG_reps.fasta" #placeholder, replace later

    OUT_DIR=$WORKDIR/OUTPUT
    QC_DIR=$OUT_DIR/QC
    ASSEMBLY_DIR=$OUT_DIR/ASSEMBLY
    BINNING_DIR=$OUT_DIR/BINNING
    BB_DIR=$OUT_DIR/BLOBOLOGY
    KRAKEN2_DIR=$OUT_DIR/KRAKEN2
    Blast_db=/storage/scratch/users/rj23k073/programs/BLAST/Database
    BLAST_DIR=$OUT_DIR/ANNOTATION    #Here we uses Russ' database to perform blast later
    echo "Input directory = $RAW_DATA_DIR, output = $OUT_DIR"

# Step 4: MEGA blast -----------------------------------------------------------------------------------------
#       This is a built-in metawrap module, classify_bins.  
#       Classifies all contigs in a set of bins by aligning them to the NCBI database with MEGABLAST, 
#       pruning the resulting hits, and assigning final taxonomy with taxator-k. 
#       The consensus taxonomy of each bin is called by contructing aweighted consensus tree, 
#       and traversing the tree by the most likely path.

    CLASSIFY_BINS_DIR=${SEL_BIN_DIR}/Classified_bins
    for sample_id in "${sample_names[@]}"; do
        # Ensure output directories exist
        mkdir -p "${BLAST_DIR}/${sample_id}"
        metawrap classify_bins -b ${SEL_BIN_DIR} -o ${CLASSIFY_BINS_DIR} -t 50
        #   Usage: metaWRAP classify_bins [options] -b bin_folder -o output_dir
        #   Options:

        #    -b STR          folder with the bins to be classified (in fasta format)
        #    -o STR          output directory
        #    -t INT          number of threads
        # Performance check for each BLAST command
            if [ $? -eq 0 ]; then
                echo "Step 4: classify_bin for ${sample_id} completed successfully, output files at ${BLAST_DIR}/${sample_id}/${Output_name}"
            else
                echo "Step 4: classify_bin for ${sample_id} has failed. Check the error report for more details."
                # Consider whether to exit or not based on your preference for handling errors
            fi
    done

        #   Performance check
            if [ $? -eq 0 ]; then
            echo "Step 4: classify_bins module completed successfully, output files at ${CLASSIFY_BINS_DIR}"
            else
            echo "Step 4: classify_bins module has failed !! Check the error report for more details."
            exit 1
            fi

# Step 4.5: Regular blast ------------------------------------------------------------------------------------
#       This is a metawrap independent section, only utilizes blastn, instead of megablast.
#       Choose either one of these to perform. 

        # Input a bin of your choice here.placeholder.
        # I have not figured out a way to perform blastn on all bins without them being shut down
        # due to limited computation cost.
        Selected_best_bin="" 
        Output_name=""


#       I cannot make them run in a loop. A tiny bin already takes ages.
#       So, make sure to select one bins per sample, and we can run them for each sample
#       For example, you can choose the longest bin. 
#       This is not a great way to deal with them, but we need a method to bypass the limits...

        for sample_id in "${sample_names[@]}"; do
            # Ensure output directories exist
            mkdir -p "${BLAST_DIR}/${sample_id}"
    
            ${blastn_script} -num_threads 50\
            -db /storage/scratch/users/rj23k073/programs/BLAST/Database/nt\  #placeholder, input blastdb
            -outfmt 6\
            -query ${Selected_best_bin} > ${BLAST_DIR}/${sample_id}/${Output_name}

                # Performance check for each BLAST command
                    if [ $? -eq 0 ]; then
                        echo "Step 4.5: blastn for ${sample_id} completed successfully, output files at ${BLAST_DIR}/${sample_id}/${Output_name}"
                    else
                        echo "Step 4.5: blastn for ${sample_id} has failed. Check the error report for more details."
                        # Consider whether to exit or not based on your preference for handling errors
                    fi
        done 

                    #   Performance check
                    if [ $? -eq 0 ]; then
                    echo "Step 4。5: blastn module completed successfully, output files at ${CLASSIFY_BINS_DIR}"
                    else
                    echo "Step 4.5: blastn module has failed !! Check the error report for more details."
                    exit 1
                    fi
