#INITIALISE FOLDERS
my_folder=/nfs/research/marioni/alsu
out_folder=${my_folder}/clust_out/lung
err_folder=${my_folder}/clust_err/lung

#CHOOSE PARAMETERS
#RAM in megabytes
memory=90000


for hvg_file in "all" "wt"
do
    script_name=scvi_sup
    bsub -e ${err_folder}/${script_name} \
    -o ${out_folder}/${script_name} \
    -q gpu -gpu "num=1:gmem=10000" -M 90000 -J ${script_name} \
    "source /hps/software/users/marioni/andrian/miniconda3/bin/activate && conda activate /nfs/research/marioni/andrian/conda_environment/test_env && python3 /nfs/research/marioni/alsu/hubmap_metaRef/miloDE_analysis/generate_data/dev/embedding_analysis/add_scvi/add_scarches_counts.py $hvg_file"
done
