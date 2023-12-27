# MULTICOM4
This is the development repository for the MULTICOM protein structure prediction system for the CASP16 competition. It includes all the programs used in the CASP15 experiment and new programs developed for the CASP16 experiment. 

Before CASP16:

1. Add new predictors into MULTICOM4

     a. Pawan's predictor

     b. DeepFold (https://github.com/newtonjoo/deepfold), docker-based?

     c. MEGAFold (https://gitee.com/mindspore/mindscience/tree/master/MindSPONGE/applications/MEGAProtein)

     d. ESMFold
   
3. Update AlphaFold software (if available)
4. Update sequence/template databases

     a. AlphaFold2 databases (bfd, mgnify, uniref90, uniref30, uniprot, pdb_seqres, pdb70, pdb_mmcifs)

     b. Additional sequence databases from DeepMSA (TaraDB (https://zenodo.org/records/3380712), MetaSourceDB(https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/README.txt?))

     c. Inhouse template databases (pdb_sort90, pdb_complex)

     d. AlphaFoldDB template database (https://github.com/google-deepmind/alphafold/tree/main/afdb) ~ 23TB?
     
5. Update interaction databases (lower priority)

     a. String database
     b. Uniprot to pdb id mapping

6. Possible tools

     a. Improving deep learning protein monomer and complex structure prediction using DeepMSA2 with huge metagenomics data.  Wei Zheng,  Qiqige Wuyun, Yang Li, Chengxin Zhang, P.  Lydia Freddolino*, and Yang Zhang. Nature Methods, in press.
