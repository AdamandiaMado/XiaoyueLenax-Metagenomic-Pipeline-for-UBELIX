# Metagenomic-Pipeline-for-UBELIX
This is a *brief* readme file explaining the main things about this pipeline, that you should know before running the script.
A more detailed explaination about packages, their purposes and results obtained is in the report.

<!-- GETTING STARTED -->
## Getting Started
This pipeline at the current stage consists of four main processes in metagenomic analysis: **QC**, **assembly**, **binning** and **alignment**. 
<p align="center"><img src="https://github.com/XiaoyueLenax/Metagenomic-Pipeline-for-UBELIX/assets/122524295/1196dc53-f2cc-40c8-881e-441854bbd3e5" alt="alt text" width="600"/></p>

A goal is to automatize this process so we get a result from raw inputs to _mature_ outputs. This is a bit challenging due to UBELIX resources on computationally intense steps such as assembly and alignments. In addition, the human gut microbiome data used in this project is very large (tens of GBs). Thus the pipeline offers two choices:
1. For people with a small data set / low number of samples: Use the all_in_one.sh to get results from beginning to end! (requires beta testing, urgently)
2. For people who are not so confident with the cluster requirements, check the individual module in the breakdown folder. There you can select specific modules you would like to test.


## Caution
1. You need to provide a directory for the package and databases, reference files, sample IDs, raw read names... etc. These are left blank in the script for you to fill in. You can quickly find them by using ctrl + F and searching for the string 'placeholder'.
2. With the test data, only breakdown versions of scripts were run because of computation power and avoiding repetition.
3. May folders will be constructed with this pipeline, but not all have their purpose at the current stage. Some modules from metaWRAP are not yet fully functional. But the folder is there if you would like to use it further. The required packages (except the databases of Kraken) are already installed and only require configuration.
4. Also make sure to change the parameters of slurm jobs based on your requirements. 
5. If any questions, please contact xiaoyue.deng@students.unibe.ch. 


## Installing metaWRAP
In the beginning, this pipeline will make a tree of folders for future use. In the breakdown this is done by the _folder_constructor.sh_ that can be executed separately. Feel free to change the format if you do not like this layout. Then, script 'metawrap_installer.sh' will install metaWRAP in the working directory given by the user.  It also installs metaWRAp if you fill in the name of the environment in the script - However, as advised by the metaWRAP author, the best way is still to manually configure metaWRAP. 
**how_to_install_metawrap.txt** contains a brief description of how this can be done.

## Running the script
The options and parameters available for each module are documented in the script as well as the metaWRAP interface. For detailed instructions, please refer to the script.
I hope you can find this pipeline useful for metagenomic analysis on UBELIX.
