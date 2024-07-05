#!/bin/bash

# cd into dance directory
data_dir="./data"
out_dir="./output/RNA2ADT/Lymph_3slices"
cvs=('cv1' 'cv2' 'cv3')

# Loop through each folder using their names as keys
for cv in "${cvs[@]}"; do
    # Fetch the path from the associative array
    tgtpath="${data_dir}/openproblems_bmmc_cite_phase2_rna" 
    curpath="${data_dir}/Lymph_${cv}"

    # Temporary name
    # echo $curpath
    # echo $tgtpath

    # Rename folder to temporary name
    mv "${curpath}" "${tgtpath}"

    # output folder
    cur_out_dir="${out_dir}_${cv}"
    mkdir -p "${cur_out_dir}"
    echo $cur_out_dir

    CUDA_VISIBLE_DEVICES=0 python ./examples/multi_modality/predict_modality/babel.py --subtask openproblems_bmmc_cite_phase2_rna --device cuda --o "${cur_out_dir}"

    # Rename back to original name
    mv "${tgtpath}" "${curpath}"
done
