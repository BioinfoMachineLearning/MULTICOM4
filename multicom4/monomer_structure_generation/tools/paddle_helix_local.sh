source /bmlfast/bml_casp16/anaconda3/bin/activate base

conda activate paddle

cd /bmlfast/bml_casp16/tools/PaddleHelix/apps/protein_folding/helixfold-single

export LD_LIBRARY_PATH=/bmlfast/bml_casp16/anaconda3/envs/paddle/lib/:$LD_LIBRARY_PATH

python helixfold_single_inference.py --init_model=./helixfold-single.pdparams --fasta_file $1 --output_dir $2