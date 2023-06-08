#!/bin/sh
#INITIALISE FOLDERS
my_folder=/nfs/research/marioni/alsu
out_folder=${my_folder}/clust_out/lung
err_folder=${my_folder}/clust_err/lung

#SELECT SCRIPT
#If you change this, you MUST update the wrapper's grep
script_name=dadea

#CHOOSE PARAMETERS
#RAM in megabytes
memory=600000
r_command="rusage[mem=${memory}]"
#num_processors
nproc=24


smg=/hps/software/users/marioni/alsu/singularity/alsu_miloDE.simg
script=/nfs/research/marioni/alsu/hubmap_metaRef/miloDE_analysis/analysis/da_de/sims_splatter/run_rmd.R

bsub -q bigmem -e ${err_folder}/${script_name} \
-o ${out_folder}/${script_name} \
-M $memory -R $r_command -n $nproc -J ${script_name} \
"singularity exec $smg Rscript $script"