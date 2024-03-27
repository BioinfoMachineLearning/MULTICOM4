source /cluster/pixstor/chengji-lab/bml_casp16/anaconda3/bin/activate base

conda activate megafold

cd /cluster/pixstor/chengji-lab/bml_casp16/tools/mindscience/MindSPONGE/applications/MEGAProtein/

export LD_LIBRARY_PATH=/cluster/pixstor/chengji-lab/bml_casp16/anaconda3/envs/megafold/lib/:$LD_LIBRARY_PATH

python main.py --model_config ./config/model.yaml --run_platform GPU --input_path $1 --checkpoint_path ./MEGA_Fold_1.ckpt --output_dir $2 --data_config $3

cp $2/seq/* $2