#!/bin/sh
#INITIALISE FOLDERS
my_folder=/nfs/research/marioni/alsu
out_folder=${my_folder}/clust_out/lung
err_folder=${my_folder}/clust_err/lung

#SELECT SCRIPT
#If you change this, you MUST update the wrapper's grep
script_name=add_mnn

#CHOOSE PARAMETERS
#RAM in megabytes
memory=300000
r_command="rusage[mem=${memory}]"
#num_processors
nproc=6


smg=/hps/software/users/marioni/alsu/singularity/alsu_miloDE.simg
script=/nfs/research/marioni/alsu/hubmap_metaRef/miloDE_analysis/generate_data/dev/embedding_analysis/add_mnn/run_rmd.R

bsub -q bigmem -e ${err_folder}/${script_name} \
-o ${out_folder}/${script_name} \
-M $memory -R $r_command -n $nproc -J ${script_name} \
"singularity exec $smg Rscript $script"
