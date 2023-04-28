#INITIALISE FOLDERS
my_folder=/nfs/research/marioni/alsu
out_folder=${my_folder}/clust_out/lung
err_folder=${my_folder}/clust_err/lung

#CHOOSE PARAMETERS
#RAM in megabytes
memory=100000

script_name=add_scarch_mouse_chim
bsub -o az.out -q gpu -gpu "num=1:gmem=10000" -M 100000 -R rusage[mem=100000] "python "python /nfs/research/marioni/alsu/hubmap_metaRef/am_hubmapMetaRef/generate_data/datasets/mouse_embryo_chimera/add_scarches/add_scarches.py"

