# MULTICOM4
This is the installation note for the tools


- DeepFold (https://github.com/newtonjoo/deepfold), surprisingly can use our conda environment to run

- MEGAFold (https://gitee.com/mindspore/mindscience/tree/master/MindSPONGE/applications/MEGAProtein)

     ```
     # Installation

     conda create -n magafold python=3.7

     conda install -c anaconda cudnn 

     conda install nvidia::cuda-nvcc

     conda install mindspore -c mindspore -c conda-forge 

     python -c "import mindspore;mindspore.set_context(device_target='GPU');mindspore.run_check()" 

     cd /bmlfast/bml_casp16/tools/mindscience/MindSPONGE/output/ 

     pip install mindsponge_gpu-1.0.0rc2-py3-none-any.whl 

     cd  /bmlfast/bml_casp16/tools/mindscience/MindSPONGE/ 

     pip install -r requirements.txt 

     mamba install -y -c bioconda hhsuite==3.3.0 kalign2

     ```

- ESMFold (https://github.com/facebookresearch/esm)

     ```
     # Installation
     cd /bmlfast/bml_casp16/tools/esm

     conda env create -f environment.yml

     conda activate esmfold

     pip install "fair-esm[esmfold]"

     pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

     ```

- Paddle-Helix (https://github.com/PaddlePaddle/PaddleHelix/tree/dev)

     ```
     # Installation
     conda create -n paddle python=3.7
     conda activate paddle
     cd /bmlfast/bml_casp16/tools/PaddleHelix/apps/protein_folding/helixfold-single
     pip install paddlepaddle_gpu-2.4.2.post117-cp37-cp37m-linux_x86_64.whl
     conda install ml-collections dm-tree biopython scipy

     python helixfold_single_inference.py --init_model=./helixfold-single.pdparams --fasta_file=data/7O9F_B.fasta --output_dir="./output"

     # If it returns error of cannot find /usr/local/cuda/lib64/libcudnn.so
     conda install -c anaconda cudnn
     export LD_LIBRARY_PATH=YOUR_ENV_DIR/lib/:$LD_LIBRARY_PATH

     ```     

- Foldseek (need to rebuild the template database) 
     ```
     conda install -c conda-forge -c bioconda foldseek
     ```

- DeepMSA2 (https://zhanggroup.org/DeepMSA/download/)

     ```
     # Installation
     conda create -n DeepMSA2 python=3.8
     conda activate DeepMSA2
     conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch

     conda create -n deepmsa2af python=3.8
     conda activate deepmsa2af
     conda update -n base conda
     conda install -y -c conda-forge openmm==7.5.1 cudnn==8.2.1.32 cudatoolkit==11.3.1 cudatoolkit-dev==11.3.1 pdbfixer==1.7
     conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04
     pip install --upgrade jax==0.2.14 jaxlib==0.1.69+cuda111 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
     pip install absl-py==0.13.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.4 dm-tree==0.1.6 immutabledict==2.0.0  ml-collections==0.1.0 numpy==1.19.5 scipy==1.7.0 tensorflow==2.5.0 pandas==1.3.4 tensorflow-cpu==2.5.0
     cd /bmlfast/bml_casp16/anaconda3/envs/deepmsa2af/lib/python3.8/site-packages/
     patch -p0 < /bmlfast/bml_casp16/tools/DeepMSA2/bin/alphafold/docker/openmm.patch
     ```

- CombFold (https://github.com/dina-lab3D/CombFold)

     ```
     cd /bmlfast/bml_casp16/tools/
     git clone https://github.com/dina-lab3D/CombFold.git
     cd CombFold

     # download boost library
     wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
     tar -zxvf boost_1_84_0.tar.gz
     cd /bmlfast/bml_casp16/tools/CombFold/boost_1_84_0

     ./bootstrap.sh  --prefix=/bmlfast/bml_casp16/tools/CombFold/boost_1_84_0

     ./b2

     ./b2 install

     # install software in CombFold
     cd /bmlfast/bml_casp16/tools/CombFold/CombinatorialAssembler
     
     # change the path of BOOST_INCLUDE and BOOST_LIB in Makefile:
     # vi Makefile
     # BOOST_INCLUDE = /bmlfast/bml_casp16/tools/CombFold/boost_1_84_0/include
     # BOOST_LIB = /bmlfast/bml_casp16/tools/CombFold/boost_1_84_0/lib/

     make

     # test running CombFold
     ./AF2trans.out

     # if return error: ./AF2trans.out: error while loading shared libraries: libboost_program_options.so.1.84.0: cannot open shared object file: No such file or directory, run:
     # export LD_LIBRARY_PATH=/bmlfast/bml_casp16/tools/CombFold/boost_1_84_0/lib
     
     ```
