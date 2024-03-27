source /cluster/pixstor/chengji-lab/bml_casp16/tools/gate/mambaforge/bin/activate

conda activate gate

export PYTHONPATH=/cluster/pixstor/chengji-lab/bml_casp16/tools/gate

cd /cluster/pixstor/chengji-lab/bml_casp16/tools/gate

mode=$1

if [ "$1" = "monomer" ]; then
    python inference_monomer.py --fasta_path $2 --input_model_dir $3 --output_dir $4 --contact_map_file $5 --dist_map_file $6
else
    python inference_multimer.py --fasta_path $2 --input_model_dir $3 --output_dir $4 --pkldir $5 --use_af_feature=True
fi
