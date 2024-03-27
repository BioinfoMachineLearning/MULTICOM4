source /cluster/pixstor/chengji-lab/bml_casp16/anaconda3/bin/activate base

conda activate esmfold

#python /cluster/pixstor/chengji-lab/bml_casp16/MULTICOM4/multicom4/monomer_structure_generation/tools/esmfold.py --fasta $1 --outpdb $2 --num_recycle $3 --chunk_size $4

python /cluster/pixstor/chengji-lab/bml_casp16/MULTICOM4/multicom4/monomer_structure_generation/tools/esmfold.py --fasta $1 --outpdb $2 --num_recycle $3