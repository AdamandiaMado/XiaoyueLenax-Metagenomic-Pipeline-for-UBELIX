#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --job-name=folder_tree_building
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_folder_tree_building_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_folder_tree_building_%j.e
#SBATCH --partition=epyc2


# -------------------------------------------------------------------------------------------

#          Setting up pipeline structure + install metawrap packages

#                                                           Reference runtime: 00:19:22 

# -------------------------------------------------------------------------------------------
# Data Strutcure -  makes directories if they do not already exist.
# Replace the work directory here with yours, placeholder.
WORKDIR=/storage/scratch/users/xd22m086/98_temp_dir/testing_wd


cd $WORKDIR
    
    mkdir -p Metawrap_Pipeline
    mkdir -p Scripts
    mkdir -p OUTPUT
    mkdir -p db

        # Here, where this script and the other scripts should be
        scp /storage/scratch/users/xd22m086/04_metawrap_testground/1_scripts_meta/meta_all.sh .
        scp /storage/scratch/users/xd22m086/04_metawrap_testground/1_scripts_meta/ncbi.sh .
        

    cd Metawrap_Pipeline
        mkdir -p OUTPUT; mkdir -p metaWRAP;
    cd ../


    cd metaWRAP

        # Install the metawrap pipeline
        git clone https://github.com/bxlab/metaWRAP.git
        echo 'export PATH="${WORKDIR}/metaWRAP/bin:$PATH"' >> ~/.bashrc


        # Here, echo and check whether metawrap is correctly configured and you can see it in the PATH.
        echo $PATH


            # In case you need to set up the environment on your own -- check the metawrap documentation
            # In case the documentation is too confusing, use this:
            #installations
            your_env_name="Insert_here"
            
            conda install -y mamba
            mamba create -y -n ${your_env_name} python=2.7
            conda activate your_env_name

            # configure channels: may not be necessary
            conda config --add channels defaults
            conda config --add channels conda-forge
            conda config --add channels bioconda
            conda config --add channels ursky

            # Install metawrap dependencies - may not be necessary
            mamba install --only-deps -c ursky metawrap-mg

            # Install metawrap
            mamba install --metawrap

            # Important step here after installation: Manually locate WIP (00 file)
    cd ../
    cd db
        # Caution: If you already have these databases installed, soft link it for faster performance, 
        #          since these databases are very large.
        mkdir -p kraken2
            cd kraken2
            # WIP - It is undecided whether to install the kraken or not, since with other packages storage is already full.
            # kraken-build --standard --threads 24 --db MY_KRAKEN_DATABASE
            # kraken-build --db MY_KRAKEN_DATABASE --clean
            #           Update your db link to the metawrap config file.
            # KRAKEN_DB=/path/to/my/database/MY_KRAKEN_DATABASE
            cd ../
        mkdir -p NCBI_nt
            cd NCBI_nt
                # Softlink the database from the repository, avoiding double download, while the account is still active.
                ln -s /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt/* .
            cd ../

        mkdir -p NCBI_tax
            cd NCBI_tax

                # softlink from repository
                ln -s /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_tax/* .
            cd ../

        mkdir -p BMTAGGER_INDEX
        mkdir -p PROGRAMS
            cd PROGRAMS

            # Smaller program scripts are copied dicrectly for convinence. 
                scp /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt/ncbi-blast-2.14.0+/bin/blastn .
                blastn_script='/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt/ncbi-blast-2.14.0+/bin/blastn'
            cd ../
    
    cd OUTPUT
        mkdir -p QC
            cd QC
                mkdir -p QUAST
            cd ../
        
        mkdir -p ASSEMBLY
            cd ASSEMBLY
                mkdir -p MetaSPAdes
                mkdir -p MEGAHIT
            cd ../
       
        mkdir -p BINNING
        mkdir -p BLOBOLOGY
        mkdir -p KRAKEN2
            cd BINNING
                mkdir -p Metabat2
                mkdir -p Maxbin2
                mkdir -p CONCOCT
            cd ../

        mkdir -p ANNOTATION
            mkdir -p blastn
            mkdir -p megablast
            cd ../
        
# Sub folders set up complete
echo "Finished Setting up folder structures."