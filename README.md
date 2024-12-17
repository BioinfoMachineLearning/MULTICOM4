# **MULTICOM4**

MULTICOM4 protein structure prediction system achieved remarkable success in the 16th world-wide Critical Assessment of Techniques for Protein Structure Prediction (#CASP16), ranking 1st in protein complex structure prediction without stoichiometry information (Phase 0), 2nd in protein tertiary structure prediction, 2nd in estimating the global fold accuracy of protein complex structures, 3rd in protein complex structure prediction with stoichiometry information (Phase 1), and 5th in protein-ligand structure and binding affinity prediction.

## **Overall workflow for the MULTICOM Protein structure prediction system**

![MULTICOM4](imgs/MULTICOM4.png)

# **Download MULTICOM4 package**

```
git clone --recursive https://github.com/BioinfoMachineLearning/MULTICOM4
```

# **Installation (non Docker version)**

## **Install AlphaFold/AlphaFold-Multimer and other required third-party packages (modified from [alphafold_non_docker](https://github.com/kalininalab/alphafold_non_docker))**

### **Install miniconda**

``` bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
```

### **Create a new conda environment and update**

``` bash
conda create --name multicom4 python==3.8
conda update -n base conda
```

### **Activate conda environment**

``` bash
conda activate multicom4
```

### **Install dependencies**

- Change `cudatoolkit==11.2.2` version if it is not supported in your system

``` bash
conda install -y -c conda-forge openmm==7.5.1 cudatoolkit==11.2.2 pdbfixer
conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2
```

- Change `jaxlib==0.3.25+cuda11.cudnn805` version if this is not supported in your system

``` bash
pip install absl-py==1.0.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.9 dm-tree==0.1.6 immutabledict==2.0.0 jax==0.3.25 ml-collections==0.1.0 numpy==1.21.6 pandas==1.3.4 protobuf==3.20.1 scipy==1.7.0 tensorflow-cpu==2.9.0

pip install --upgrade --no-cache-dir jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### **Download chemical properties to the common folder**

``` bash

# Replace $MULTICOM4_INSTALL_DIR with your MULTICOM4 installation directory

wget -q -P $MULTICOM4_INSTALL_DIR/tools/alphafold-v2.3.2/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
```

### **Apply OpenMM patch**

``` bash
# Replace $MULTICOM4_INSTALL_DIR with your MULTICOM4 installation directory

cd ~/anaconda3/envs/multicom4/lib/python3.8/site-packages/ && patch -p0 < $MULTICOM4_INSTALL_DIR/tools/alphafold-v2.3.2/docker/openmm.patch

# or

cd ~/miniconda3/envs/multicom4/lib/python3.8/site-packages/ && patch -p0 < $MULTICOM4_INSTALL_DIR/tools/alphafold-v2.3.2/docker/openmm.patch
```

### **Install other required third-party packages**

```
conda install tqdm
conda install -c conda-forge -c bioconda mmseqs2=14.7e284 -y
```

### **Download Genetic databases in AlphaFold2/AlphaFold-Multimer**

```
# Replace $MULTICOM4_INSTALL_DIR with your MULTICOM4 installation directory

bash $MULTICOM4_INSTALL_DIR/tools/alphafold-v2.3.2/scripts/download_all_data.sh <YOUR_ALPHAFOLD_DB_DIR>
```

### **Install the MULTICOM4 addon system and its databases**

```
# Note: here the parameters should be the absolute paths
python download_database_and_tools.py --multicom4db_dir <YOUR_MULTICOM4_DB_DIR>

# Configure the MULTICOM4 system
# Replace $MULTICOM4_INSTALL_DIR with your MULTICOM4 installation directory
# Replace $YOUR_ALPHAFOLD_DB_DIR with your downloaded AlphaFold databases directory

python configure.py --envdir ~/miniconda3/envs/multicom4 --multicom4db_dir <YOUR_MULTICOM4_DB_DIR> --afdb_dir <YOUR_ALPHAFOLD_DB_DIR>

# e.g, 
# python download_database_and_tools.py \
# --multicom4db_dir /home/multicom4/tools/multicom4_db

# python configure.py \
# --multicom4db_dir /home/multicom4/tools/multicom4_db \
# --afdb_dir /home/multicom4/tools/alphafold_databases/
```
The configure.py python script will 
* Copy the alphafold_addon scripts
* Create the configuration file (bin/db_option) for running the system

# **Genetic databases used by MULTICOM4**

Assume the following databases have been installed as a part of the AlphaFold2/AlphaFold-Multimer installation
*   [BFD](https://bfd.mmseqs.com/),
*   [MGnify](https://www.ebi.ac.uk/metagenomics/),
*   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
*   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
*   [PDB seqres](https://www.rcsb.org/)
*   [UniRef30](https://uniclust.mmseqs.com/),
*   [UniProt](https://www.uniprot.org/uniprot/),
*   [UniRef90](https://www.uniprot.org/help/uniref).

Additional databases will be installed for the MULTICOM system by setup.py:
*   [AlphaFoldDB](https://alphafold.ebi.ac.uk/): ~53G
*   [Metaclust](https://metaclust.mmseqs.org/current_release/): ~114G
*   [STRING](https://string-db.org/cgi/download?sessionId=bgV6D67b9gi2): ~129G
*   [pdb_complex](https://www.biorxiv.org/content/10.1101/2023.05.16.541055v1): ~38G
*   [pdb_sort90](https://www.biorxiv.org/content/10.1101/2023.05.01.538929v1): ~48G
*   [Uniclust30](https://uniclust.mmseqs.com/): ~87G

# **Important parameter values in db_option**

```
# AlphaFold2 parameters
monomer_num_ensemble = 1
monomer_num_recycle = 3
num_monomer_predictions_per_model = 1
monomer_model_preset = monomer

# AlphaFold-Multimer parameters
multimer_num_ensemble = 1
multimer_num_recycle = 3
num_multimer_predictions_per_model = 5
multimer_model_preset = multimer

# Common parameters
alphafold_benchmark = True
use_gpu_relax = True
models_to_relax = ALL # ALL, BEST, NONE
max_template_date = 2024-06-01
```
Please refer to [AlphaFold2](https://github.com/deepmind/alphafold) to understand the meaning of the parameters. The parameter values stored in bin/db_option file are applied to all the AlphaFold2/AlphaFold-Multimer variants in the MULTICOM4 system to generate predictions. 

For Docker version installation, you can change the default parameter values in [docker/db_option](docker/db_option).

For non Docker version of the installation, the default bin/db_option file is created automatically by configure.py during the installation. The default parameter values above can be changed if needed. 

# **Before running the system for non Docker version**

## **Activate your python environment and add the MULTICOM4 system path to PYTHONPATH**

```bash
conda activate multicom4

# Replace $MULTICOM4_INSTALL_DIR with your MULTICOM4 installation directory (absolute path)
export PYTHONPATH=$MULTICOM4_INSTALL_DIR

# e.g, 
# conda activate MULTICOM4
# export PYTHONPATH=/home/multicom4/MULTICOM4

```
Now MULTICOM4 is ready for you to make predictions.

# **Running the monomer/tertiary structure prediction pipeline**

Say we have a monomer with the sequence `<SEQUENCE>`. The input sequence file should be in the FASTA format as follows:

```fasta
>sequence_name
<SEQUENCE>
```

Note: It is recommended that the name of the sequence file in FASTA format should be the same as the sequence name.

Then run the following command:

```bash
# Please provide absolute path for the input parameters
python bin/monomer.py \
    --option_file=bin/db_option \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --output_dir=$OUTDIR
```

option_file is a file in the MULTICOM package to store some key parameter values for AlphaFold2 and AlphaFold-Multimer. fasta_path is the full path of the file storing the input protein sequence(s) in the FASTA format. output_dir specifies where the prediction results are stored. Please be aware that we have included a parameter (--run_img) that allows you to turn off the usage of the IMG database for faster prediction (--run_img=False). In the case of --run_img=True, the program will pause at the monomer prediction generation stage to wait for the IMG alignment to be created. Generating alignments from IMG may take a much longer time, potentially several days, because the database is very large. So run_img is set to false by default. It is advised that run_img is set to true only if other alignments cannot yield good results.

## **Output**

```
$OUTPUT_DIR/                                   # Your output directory
    N1_monomer_alignments_generation/          # Working directory for generating monomer MSAs
    N2_monomer_template_search/                # Working directory for searching monomer templates
    N3_monomer_structure_generation/           # Working directory for generating monomer structural predictions
    N4_monomer_structure_evaluation/           # Working directory for evaluating the monomer structural predictions
        - alphafold_ranking.csv    # AlphaFold2 pLDDT ranking
        - pairwise_ranking.tm      # Pairwise (APOLLO) ranking
        - pairwise_af_avg.ranking  # Average ranking of the two
```

* The predictions and ranking files are saved in the *N4_monomer_structure_evaluation* folder. You can check the AlphaFold2 pLDDT score ranking file (alphafold_ranking.csv) to look for the structure with the highest pLDDT score. The *pairwise_ranking.tm* and *pairwise_af_avg.ranking* are the other two ranking files. 

# **Running the multimer/quaternary structure prediction pipeline**

## **Folding a multimer**

Say we have a homomer with 4 copies of the same sequence
`<SEQUENCE>`. The input file should be in the format as follows:

```fasta
>sequence_1
<SEQUENCE>
>sequence_2
<SEQUENCE>
>sequence_3
<SEQUENCE>
>sequence_4
<SEQUENCE>
```

Then run the following command:

```bash
# Please provide absolute path for the input parameters
python bin/multimer.py \
    --option_file=bin/db_option \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --output_dir=$OUTDIR
```

## **Output**

```
$OUTPUT_DIR/                                   # Your output directory
    N1_monomer_alignments_generation/          # Working directory for generating monomer MSAs
        - Subunit A
        - Subunit B
        - ...
    N1_monomer_alignments_generation_img/      # Working directory for generating IMG MSA
        - Subunit A
        - Subunit B
        - ...
    N2_monomer_template_search/                # Working directory for searching monomer templates
        - Subunit A
        - Subunit B
        - ...
    N3_monomer_structure_generation/           # Working directory for generating monomer structural predictions
        - Subunit A
        - Subunit B
        - ...
    N4_monomer_alignments_concatenation/       # Working directory for concatenating the monomer MSAs
    N5_monomer_templates_search/               # Working directory for concatenating the monomer templates
    N6_multimer_structure_generation/          # Working directory for generating multimer structural predictions
    N7_monomer_structure_evaluation            # Working directory for evaluating monomer structural predictions
        - Subunit A
            # Rankings for all the predictions
            - alphafold_ranking.csv            # AlphaFold2 pLDDT ranking 
            - pairwise_ranking.tm              # Pairwise (APOLLO) ranking
            - pairwise_af_avg.ranking          # Average ranking of the two 

            # Rankings for the predictions generated by monomer structure prediction
            - alphafold_ranking_monomer.csv    # AlphaFold2 pLDDT ranking 
            - pairwise_af_avg_monomer.ranking  # Average ranking 

            # Rankings for the predictions extracted from multimer predictions
            - alphafold_ranking_multimer.csv   # AlphaFold2 pLDDT ranking 
            - pairwise_af_avg_multimer.ranking # Average ranking 

        - Subunit B
        - ...
    N8_multimer_structure_evaluation           # Working directory for evaluating multimer structural predictions
        - alphafold_ranking.csv                # AlphaFold2 pLDDT ranking
        - multieva.csv                         # Pairwise ranking using MMalign
        - pairwise_af_avg.ranking              # Average ranking of the two
```

* The predictions and ranking files are saved in *N8_multimer_structure_evaluation*, similarly, you can check the AlphaFold-Multimer confidence score ranking file (alphafold_ranking.csv) to look for the structure with the highest predicted confidence score generated by AlphaFold-Multimer. The *multieva.csv* and *pairwise_af_avg.ranking* are the other two ranking files.

* The monomer structures and ranking files are saved in *N7_monomer_structure_evaluation* if you want to check the predictions and rankings for the monomer structures.

# **Some CASP16 Prediction Examples**




