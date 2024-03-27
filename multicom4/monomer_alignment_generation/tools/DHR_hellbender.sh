source /bmlfast/bml_casp16/anaconda3/bin/activate base

conda activate alphafold_DHR

cd /cluster/pixstor/chengji-lab/bml_casp16/tools/Dense-Homolog-Retrieval

python retrieve.py --input_path $1 --database_path $2 --output_a3m_path $3
