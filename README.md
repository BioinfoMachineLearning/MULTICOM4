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

     a. AlphaFold2 databases (Updated on Hellbender on December 27)
   
     *   [BFD](https://bfd.mmseqs.com/),
     *   [MGnify v2023_02]([https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/),
     *   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
     *   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
     *   [PDB seqres](https://www.rcsb.org/)
     *   [UniRef30 v2023_02](https://gwdu111.gwdg.de/~compbiol/uniclust/2023_02/),
     *   [UniProt](https://www.uniprot.org/uniprot/),
     *   [UniRef90](https://www.uniprot.org/help/uniref).

     b. Additional sequence databases from DeepMSA
          
     *   [TaraDB](https://zenodo.org/records/3380712)
     *   [MetaSourceDB v2023_02](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/)
   
     c. Inhouse template databases (pdb_sort90, pdb_complex)

     d. AlphaFoldDB template database (https://github.com/google-deepmind/alphafold/tree/main/afdb) ~ 23TB?
     
6. Update interaction databases (lower priority)

     a. [String database v12.0](https://string-db.org/cgi/download?sessionId=bZgGNwyipdWy)
 
     b. Uniprot to pdb id mapping

7. Possible tools

     a. Improving deep learning protein monomer and complex structure prediction using DeepMSA2 with huge metagenomics data.  Wei Zheng,  Qiqige Wuyun, Yang Li, Chengxin Zhang, P.  Lydia Freddolino*, and Yang Zhang. Nature Methods, in press.
