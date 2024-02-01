# MULTICOM4
This is the development repository for the MULTICOM protein structure prediction system for the CASP16 competition. It includes all the programs used in the CASP15 experiment and new programs developed for the CASP16 experiment. 

Before CASP16:

1. Add new predictors into MULTICOM4

     a. Pawan's predictor (integrated into the system)

     b. DeepFold (https://github.com/newtonjoo/deepfold), surprisingly can use our conda environment to run

     c. MEGAFold (https://gitee.com/mindspore/mindscience/tree/master/MindSPONGE/applications/MEGAProtein)

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

     d. ESMFold (https://github.com/facebookresearch/esm)

     ```
     # Installation
     cd /bmlfast/bml_casp16/tools/esm

     conda env create -f environment.yml

     conda activate esmfold

     pip install "fair-esm[esmfold]"

     pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

     ```

     e. Paddle-Helix (https://github.com/PaddlePaddle/PaddleHelix/tree/dev)

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
     
     ~~f. DMFold (https://zhanggroup.org/DMFold/)~~

2. Update tools 

     a. Foldseek (need to rebuild the template database) 
     ```
     conda install -c conda-forge -c bioconda foldseek
     ```

     b. DeepMSA2 (https://zhanggroup.org/DeepMSA/download/)

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

     c. AlphaFold software (if available)
   
3. Update sequence/template databases

     a. AlphaFold2 databases (Updated on Hellbender on December 27)
   
     *   [BFD](https://bfd.mmseqs.com/),
     *   [MGnify v2023_02]([https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/),
     *   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
     *   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
     *   [PDB seqres](https://www.rcsb.org/)
     *   [UniRef30 v2023_02](https://gwdu111.gwdg.de/~compbiol/uniclust/2023_02/),
     *   [UniProt](https://www.uniprot.org/uniprot/),
     *   [UniRef90](https://www.uniprot.org/help/uniref).

     b. ~~DeepMSA2 (tianqi)~~
          
     ~~*   [TaraDB](https://zenodo.org/records/3380712)~~

     ~~*   [MetaSourceDB v2023_02](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/)~~
   
     c. Inhouse template databases (pdb_sort90, pdb_complex)

     d. AlphaFoldDB template database (https://github.com/google-deepmind/alphafold/tree/main/afdb) ~ 23TB?

     e. ColabFold Database (2021_08) (https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh)
     
4. Update interaction databases (lower priority)

     a. [String database v12.0](https://string-db.org/cgi/download?sessionId=bZgGNwyipdWy)
 
     ```
     # Update

     cd /bmlfast/bml_casp16/databases

     wget https://stringdb-downloads.org/download/protein.links.v12.0.txt.gz

     gunzip -f protein.links.v12.0.txt.gz

     wget https://stringdb-downloads.org/download/protein.aliases.v12.0.txt.gz

     gunzip -f protein.aliases.v12.0.txt.gz

     python uniprot2string.py --string_links_file protein.links.v12.0.txt --string_aliases_file protein.aliases.v12.0.txt --string2uniprot_map string2uniprot.map

     ```

     b. Uniprot to pdb id mapping

5. Possible tools

6. Current predictors

     a. Monomer (31 predictors)

     | Predictor  | Note |
     | -------------| -----| 
     |default | default AlphaFold2 |
     |default_seq_temp | Inhouse template database|
     |original | search MSA seperately (remove?) |
     |ori_seq_temp | Inhouse template database|
     |colabfold | ColabFold Alignment | 
     |colab_seq_temp | Inhouse template database|
     |img | IMG alignment (keep?) |
     |img_seq_temp | Inhouse template database|

     Added:
     
     | Predictor  | Note |
     | -------------| -----| 
     |Paddle | Paddle-Helix| 
     |ESMFold | ESMFold|
     |DeepFold | DeepFold |
     |MEGAFold | MEGAFold |
     |DMFold | DMFold|
     |16 kinds of MSA in DeepMSA2 | DeepMSA2| 
     |5 variants in AFsample | AFsample| 
     | Foldseek refinement as a predictor | use the top-ranked models generated by default alphafold2 |

     b. Multimer

     | Predictor  | Note |
     | -------------| -----| 
     |default_multimer | default AlphaFold-Multimer |
     |default_variants | default_variants|
     |spec_variants | spec_variants |
     |str_variants | str_variants (hetero-multimer only) |
     |pdb_variants| pdb_variants | 
     |unidist_variants | unidist_variants (hetero-multimer only) |
     |Iterative methods | Structural-alignment based methods |

     To be added:

     | Predictor  | Note |
     | -------------| -----| 
     |ESMFold | ESMFold|
     |DeepFold | DeepFold |
     |~~MEGAFold~~ | ~~MEGAFold~~ |
     |~~DMFold~~ | ~~DMFold~~|
     |? kinds of MSA in DeepMSA2 | DeepMSA2| 
     |5 variants in AFsample | AFsample| 

6. Quality assessment methods

     a. Gate

     b. GCPNet-EMA


7. To be solved:

     a. How to use multiple gpus to run a target, especially run it on the hellbender?

          i. Generate all the required MSAs on our server. For example, for monomer targets, N1, N2 needs to be run first, and for multimer targets, N1, N2, N3 (run default AlphaFold-Multimer for each subunits), N4, N5 (may have some issues for the structure-based template search pipeline) needs to be ran first.

          ii. For each gpu, run the whole system (need to check the docker-based predictors) for the target or individual predictors?

          iii. The directories need to be organized (e.g., N6~N10) to run the prediction for the same target with different gpus on the same server/node. 

          iv. Need to write some scripts to collect the models and run QA on them.

     Monomer:
     ```
     $OUTPUT_DIR/            
          N1_monomer_alignments_generation/          
          N1_monomer_alignments_generation_img/      
          N2_monomer_template_search/                
          N3_monomer_structure_generation/           
          N4_monomer_structure_evaluation/         
          N5_monomer_structure_refinement_avg/       
          N5_monomer_structure_refinement_avg_final/ 
     ```

     Multimer:
     ```
     $OUTPUT_DIR/
          N1_monomer_alignments_generation/
               - Subunit A
               - Subunit B
               - ...
          N1_monomer_alignments_generation_img/
               - Subunit A
               - Subunit B
               - ...
          N2_monomer_template_search/
               - Subunit A
               - Subunit B
               - ...
          N3_monomer_structure_generation/
               - Subunit A
               - Subunit B
               - ...
          N4_monomer_alignments_concatenation/ 
          N5_monomer_templates_search/            
          N6_multimer_structure_generation/         
          N7_monomer_structure_evaluation            
               - Subunit A
               - Subunit B
               - ...
          N8_multimer_structure_evaluation 
          N9_multimer_structure_refinement           
          N9_multimer_structure_refinement_final     
     ```
