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
#Load all vital-IT packages preinstalled
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
    sample_names=("Human1" "Human2" "Human3") # Recommend to manually type in all your sample names to loop through later..
    REF_GENOME="/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/UHGG_reps.fasta" 




# Step 3 : Binning ------------------------------------------------------------------------------------------------------------
    
    # Will automatically perform binning on both megahit and metaspades file. 
    # If you only used 1 assembler, hash out the one you're not using.

    MEGAHIT_DIR=${ASSEMBLY_DIR}/MEGAHIT
    METASPADES_DIR=${ASSEMBLY_DIR}/MetaSPAdes
    CONCOCT_DIR=${OUT_DIR}/CONCOCT
    
    # Pipeline Option 1 - Russ's Concoct, independent of metawrap
    #  For MEGAHIT
        perl /storage/scratch/users/rj23k073/programs/Maxbin2/MaxBin-2.2.7/run_MaxBin.pl
        cut_up_fasta.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta -c 10000 -o 0 --merge_last -b 6_2_deer_10K.bed > 6_2_deer_10K.fa
        concoct_coverage_table.py 6_2_deer_10K.bed /storage/scratch/users/rj23k073/04_DEER/07_BAM/6_2_filter.sorted.bam > 6_2_coverage_table.tsv
        concoct --composition_file 6_2_deer_10K.fa --coverage_file 6_2_coverage_table.tsv -b concoct_output/
        merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        mkdir concoct_output/fasta_bins
        extract_fasta_bins.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
    # For MetaSPAdes
        perl /storage/scratch/users/rj23k073/programs/Maxbin2/MaxBin-2.2.7/run_MaxBin.pl
        cut_up_fasta.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta -c 10000 -o 0 --merge_last -b 6_2_deer_10K.bed > 6_2_deer_10K.fa
        concoct_coverage_table.py 6_2_deer_10K.bed /storage/scratch/users/rj23k073/04_DEER/07_BAM/6_2_filter.sorted.bam > 6_2_coverage_table.tsv
        concoct --composition_file 6_2_deer_10K.fa --coverage_file 6_2_coverage_table.tsv -b concoct_output/
        merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        mkdir concoct_output/fasta_bins
        extract_fasta_bins.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

    for sample_id in "${sample_name[@]}"; do

        # ALL FOR METAHIT RESULTS
        # You may not want metabat 1 now there's metabat2, but this is an option
        metawrap binning -o ${BINNING_DIR}/metabat2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}/final_assembly.fasta --metabat2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        # metawrap binning -o ${BINNING_DIR}/metabat1/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}final_assembly.fasta --metabat1 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/maxbin2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}final_assembly.fasta --maxbin2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/concoct/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}final_assembly.fasta --concoct ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        
        # ALL FOR METASPADES RESULTS
        metawrap binning -o ${BINNING_DIR}/metabat2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --metabat2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        # metawrap binning -o ${BINNING_DIR}/metabat1/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --metabat1 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/maxbin2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --maxbin2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/concoct/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --concoct ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        #  Parameters used:
        #       -a STR          metagenomic assembly file
        #        -o STR          output directory
        #        -t INT          number of threads (default=1)
        #        -m INT          amount of RAM available (default=4)
        #        -l INT          minimum contig length to bin (default=1000bp). Note: metaBAT will default to 1500bp minimum

        #        --metabat2      bin contigs with metaBAT2
        #        --metabat1      bin contigs with the original metaBAT
        #        --maxbin2       bin contigs with MaxBin2
        #        --concoct       bin contigs with CONCOCT

        #        --universal     use universal marker genes instead of bacterial markers in MaxBin2 (improves Archaea binning)
        #        --run-checkm    immediately run CheckM on the bin results (requires 40GB+ of memory)
        #        --single-end    non-paired reads mode (provide *.fastq files)
        #        --interleaved   the input read files contain interleaved paired-end reads
    done
        # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3: Binning completed successfully, output files at ${BINNING_DIR}"
        else
        echo "Step 3: Binning has failed !! Check the error report for more details."
        exit 1
        fi



# Step 3.5 checkm and GC plot ---------------------------------------------------------------------------------
#                      To add or remove options, use this guide below of checkm.

#                               ...::: CheckM v1.2.2 :::...

#           Lineage-specific marker set:
#               tree         -> Place bins in the reference genome tree
#               tree_qa      -> Assess phylogenetic markers found in each bin
#               lineage_set  -> Infer lineage-specific marker sets for each bin

#           Taxonomic-specific marker set:
#               taxon_list   -> List available taxonomic-specific marker sets
#               taxon_set    -> Generate taxonomic-specific marker set

#           Apply marker set to genome bins:
#               analyze      -> Identify marker genes in bins
#               qa           -> Assess bins for contamination and completeness

#           Common workflows (combines above commands):
#               lineage_wf   -> Runs tree, lineage_set, analyze, qa
#               taxonomy_wf  -> Runs taxon_set, analyze, qa

#           Reference distribution plots:
#               gc_plot      -> Create GC histogram and delta-GC plot
#               coding_plot  -> Create coding density (CD) histogram and delta-CD plot
#               tetra_plot   -> Create tetranucleotide distance (TD) histogram and delta-TD plot
#               dist_plot    -> Create image with GC, CD, and TD distribution plots together

#           General plots:
#               nx_plot      -> Create Nx-plots
#               len_hist     -> Sequence length histogram
#               marker_plot  -> Plot position of marker genes on sequences
#               gc_bias_plot -> Plot bin coverage as a function of GC

#           Bin exploration and modification:
#               unique       -> Ensure no sequences are assigned to multiple bins
#               merge        -> Identify bins with complementary sets of marker genes
#                outliers     -> [Experimental] Identify outlier in bins relative to reference distributions
#               modify       -> [Experimental] Modify sequences in a bin

#           Utility functions:
#               unbinned     -> Identify unbinned sequences
#               coverage     -> Calculate coverage of sequences
#               tetra        -> Calculate tetranucleotide signature of sequences
#               profile      -> Calculate percentage of reads mapped to each bin
#               ssu_finder   -> Identify SSU (16S/18S) rRNAs in sequences


for sample_id in "${samples[@]}"; do
	
    checkm lineage_wf -t 50 -x fa ${BIN_DIR}/${sample_id}/metabat2_bins ${OUT_DIR}/${sample_id}/checkm
    #runs tree, lineage_set, analyze, qa
    # usage: checkm lineage_wf [-h] [-r] [--ali] [--nt] [-g] [-u UNIQUE] [-m MULTI] [--force_domain] [--no_refinement] 
    #                     [--individual_markers] [--skip_adj_correction] [--skip_pseudogene_correction] [--aai_strain AAI_STRAIN] 
    #                     [-a ALIGNMENT_FILE] [--ignore_thresholds]
    #                     [-e E_VALUE] [-l LENGTH] [-f FILE] [--tab_table] [-x EXTENSION] [-t THREADS] [--pplacer_threads PPLACER_THREADS] [-q]
    #                     [--tmpdir TMPDIR]
    #                     bin_input output_dir

    checkm gc_plot --dpi 1000 -x fa ${BIN_DIR}/${sample_id}/metabat2_bins ${OUT_DIR}/${sample_id}/checkm 1
    # gc_plot: create GC histogram and delta-GC plot
    # usage: checkm gc_plot [-h] [--image_type {eps,pdf,png,ps,svg}] [--dpi DPI] [--font_size FONT_SIZE] 
    #                   [-x EXTENSION] [--width WIDTH] [--height HEIGHT]
    #                  [-w GC_WINDOW_SIZE] [-b GC_BIN_WIDTH] [-q]
    #                  bin_input output_dir dist_value [dist_value ...]
done

    # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3.5: Checkm completed successfully, output files at ${OUT_DIR}/${sample_id}/checkm"
        else
        echo "Step 3.5: Checkm has failed !! Check the error report for more details."
        exit 1
        fi

# Step 3.6 QUANT_bins----------------------------------------------------------------------------------------------------
#       This is a build in metawrap module.
#       This module takes in a set of bins and any number of paired-end read sets from metagenomic samples, and estimates the abundance of each bin in each 
#       sample with salmon. It then uses Seaborn to make clustered heatmaps genome abundances.
#
#       Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
#       For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.


# Set here the bins that you would like to use, by default is the metabat2

    SEL_BIN_DIR=${BINNING_DIR}/metabat2
    # WIP: give an option to skil salmon?
#       For megahit
for sample_id in "${sample_name[@]}"; do
    #metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins -o ${OUT_DIR}/${sample_id} -a ${BIN_DIR}/${sample_id}/binned_assembly/assembly.fa ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins \
    -o ${OUT_DIR}/${sample_id} \
    -a ${MEGAHIT_DIR}/${sample_id}/final_assembly.fasta ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    #   Usage: metaWRAP quant_bins [options] -b bins_folder -o output_dir -a assembly.fa readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]
    #   Options:
    #        -b STR          folder containing draft genomes (bins) in fasta format
    #        -o STR          output directory
    #        -a STR          fasta file with entire metagenomic assembly (strongly recommended!)
    #        -t INT          number of threads
done

#        For metaSPADES
for sample_id in "${sample_name[@]}"; do
    #metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins -o ${OUT_DIR}/${sample_id} -a ${BIN_DIR}/${sample_id}/binned_assembly/assembly.fa ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins \
    -o ${OUT_DIR}/${sample_id} \
    -a ${METASPADES_DIR}/${sample_id}/final_assembly.fasta ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    #   Usage: metaWRAP quant_bins [options] -b bins_folder -o output_dir -a assembly.fa readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]
    #   Options:
    #        -b STR          folder containing draft genomes (bins) in fasta format
    #        -o STR          output directory
    #        -a STR          fasta file with entire metagenomic assembly (strongly recommended!)
    #        -t INT          number of threads
done

    # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3.6: Quant_bins completed successfully, output files at ${OUT_DIR}/${sample_id}/checkm"
        else
        echo "Step 3.6: Quant_bins has failed !! Check the error report for more details."
        exit 1
        fi


# Step 3.7 Bin refinement -------------------------------------------------------
#
#   This is a metawrap built in module that run on the outputs of binning.sh pipeline 
#   to analyze the metagenomic bins and arrive at the best possible putative genomes.
#   There are several options to give additional binning results for comparison. 

    for sample_id in "${sample_name[@]}"; do    
        metawrap bin_refinement -o ${SEL_BIN_DIR}/${sample_id}/refined_bins -A ${SEL_BIN_DIR}/${sample_id}/metabat2_bins -1 ${RAW_DATA_DIR}/${sample_id}/*_1.fastq -2 ${RAW_DATA_DIR}/${sample_id}/*_2.fastq -t 50

        #  Usage: metaWRAP bin_refinement [options] -o output_dir -A bin_folderA [-B bin_folderB -C bin_folderC]
        #       Note: the contig names in different bin folders must be consistant (must come from the same assembly).

        #   Options:

        #    -o STR          output directory
        #    -t INT          number of threads (default=1)
        #    -m INT          memory available (default=40)
        #    -c INT          minimum % completion of bins [should be >50%] (default=70)
        #    -x INT          maximum % contamination of bins that is acceptable (default=10)

        #    -A STR          folder with metagenomic bins (files must have .fa or .fasta extension)
        #    -B STR          another folder with metagenomic bins
        #    -C STR          another folder with metagenomic bins

        #    --skip-refinement       dont use binning_refiner to come up with refined bins based on combinations of binner outputs
        #    --skip-checkm           dont run CheckM to assess bins
        #    --skip-consolidation    choose the best version of each bin from all bin refinement iteration
        #    --keep-ambiguous        for contigs that end up in more than one bin, keep them in all bins (default: keeps them only in the best bin)
        #    --remove-ambiguous      for contigs that end up in more than one bin, remove them in all bins (default: keeps them only in the best bin)
        #    --quick                 adds --reduced_tree option to checkm, reducing runtime, especially with low memory

    done

    # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3.7: bin_refinement module completed successfully, output files at ${SEL_BIN_DIR}/${sample_id}/refined_bins"
        else
        echo "Step 3.7: bin_refinement module has failed !! Check the error report for more details."
        exit 1
        fi
