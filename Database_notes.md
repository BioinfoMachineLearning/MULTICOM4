# MULTICOM4
This is database notes for MULTICOM4

- AlphaFold2 databases (Updated on Hellbender on December 27)
   
  *   [BFD](https://bfd.mmseqs.com/),
  *   [MGnify v2023_02]([https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/),
  *   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
  *   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
  *   [PDB seqres](https://www.rcsb.org/)
  *   [UniRef30 v2023_02](https://gwdu111.gwdg.de/~compbiol/uniclust/2023_02/),
  *   [UniProt](https://www.uniprot.org/uniprot/),
  *   [UniRef90](https://www.uniprot.org/help/uniref).

- DeepMSA2
  * Overlapped
     *   [BFD](https://bfd.mmseqs.com/),
     *   [UniRef30 v2023_02](https://gwdu111.gwdg.de/~compbiol/uniclust/2023_02/),
     *   [UniRef90](https://www.uniprot.org/help/uniref).
     *   [MGnify v2023_02]([https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/),
  * Non-overlapped
     *   [Metaclust](https://metaclust.mmseqs.org/current_release/)
     *   [UniClust30](https://gwdu111.gwdg.de/~compbiol/uniclust/2018_08/)
     *   JGI (curated by DeepMSA2) 

-  Inhouse template databases (pdb_sort90, pdb_complex)

     d. AlphaFoldDB template database (https://github.com/google-deepmind/alphafold/tree/main/afdb) ~ 23TB?

     e. ColabFold Database (2021_08) (https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh)
     
- Update interaction databases (lower priority)

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