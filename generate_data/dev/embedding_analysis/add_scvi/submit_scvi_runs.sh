#INITIALISE FOLDERS
my_folder=/nfs/research/marioni/alsu
out_folder=${my_folder}/clust_out/lung
err_folder=${my_folder}/clust_err/lung

#CHOOSE PARAMETERS
#RAM in megabytes
memory=90000

#out_folder=/nfs/research1/marioni/alsu/geneBasis/data/processed_time/scmer/spleen

#cd $out_folder

n_layers=(2 3)
counts_assay=("counts_all_genefull" "counts_all_genefull_corrected" "counts_unspliced" "counts_unspliced_corrected")
ref_type=("cells nuclei")

for n_layers in 2 3
do
	for counts_assay in "counts_all_genefull" "counts_all_genefull_corrected" "counts_unspliced" "counts_unspliced_corrected"
	do
	    script_name=arc_thy
            bsub -e ${err_folder}/${script_name} \
            -o ${out_folder}/${script_name} \
            -q gpu -gpu "num=1:gmem=10000" -M 90000 -J ${script_name} \
            "source /hps/software/users/marioni/andrian/miniconda3/bin/activate && conda activate /nfs/research/marioni/andrian/conda_environment/test_env && python3 /nfs/research/marioni/alsu/snc/snc_analysis/gendata/thymus_analysis/integration/scvi/add_scvi_counts.py $n_layers $counts_assay"
	done
done
