source /bmlfast/bml_casp16/anaconda3/bin/activate base

conda activate esmfold

python /bmlfast/bml_casp16/MULTICOM4/multicom4/monomer_structure_generation/tools/esm_msa_alg1.py --ina3m $1 --outfile $2