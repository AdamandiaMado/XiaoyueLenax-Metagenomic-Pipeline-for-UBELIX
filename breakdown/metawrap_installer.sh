#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=10GB
#SBATCH --time=10:00:00
#SBATCH --job-name=meta_installer
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_installer_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_installer_%j.e
#SBATCH --partition=epyc2


# Installing metawrap, do this in your environment of choice! 
conda install -y -c ursky metawrap-mg
conda install -y blas=2.5=mkl
mamba install -y metawrap

NCBI_NT_DIR="Your database here - placeholder, Requires manual configuration"
NCBI_TAX_DIR="Your database here - placeholder, Requires manual configuration"

wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz" -P ${NCBI_NT_DIR}
    for a in nt.*.tar.gz
        do tar xzf $a
    done

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P ${NCBI_TAX_DIR}
cd ${NCBI_TAX_DIR}
    tar -xvf taxdump.tar.gz


cd ${INSTALL_DIR}
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz -P /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/BMTAGGER_INDEX
    gunzip *fa.gz
    cat *fa > hg38.fa
    rm chr*.fa

    # If using human genome: install the human database.
    #bmtool -d hg38.fa -o hg38.bitmask
    #srprism mkindex -i hg38.fa -o hg38.srprism -M 100000