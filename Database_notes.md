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

-  Inhouse template databases (pdb_sort90, pdb_complex): /bmlfast/bml_casp16/update_complex_dbs_newest

   - pdb_sort90
     ```
     # Need to change the pdb date cut off in bin/db_option_2024
     # e.g., set_pdb_release_date = 01-05-2024

     python scripts/P1_update_pdb_slim.py --option_file bin/db_option_2024
     python scripts/P2_update_mmcif.py --option_file bin/db_option_2024
     python scripts/P3_update_dataset.py --option_file bin/db_option_2024
     python scripts/P4_update_cm.py --option_file bin/db_option_2024

     # Note: Eventually we need to generate hmm files for each sequences in the database using the latest uniref30 database
     # The total number of the sequences are very large (~652719), however some sequences are identical
     # We use this script to filter the identical sequences and generate the a3m/hmm files for each sequence (~650000 to ~240000)
     python scripts/P5_generate_msas_fastas.py --fastafile databases/RCSB_PDB/cm_lib/pdb_cm.fasta --mmseq tools/mmseqs/bin/mmseqs --outputdir databases/RCSB_PDB/cm_lib/msa_fastas

     # Generate a3m/hmm for the filtered sequences
     python scripts/batch_run_hhblits.py --fasta_dir databases/RCSB_PDB/cm_lib/msa_fastas/fasta_rep --outputdir databases/RCSB_PDB/cm_lib/msa_fastas/msas_same --hhblits3_dir tools/hhsuite-3.2.0 --hhblits_db /bmlfast/bml_casp16/tools/alphafold_databases_multicom3/uniref30/UniRef30_2023_02
     
     # Generate (copy) a3m for all the sequences, the targets without generated a3ms will be saved in databases/RCSB_PDB/cm_lib/msa_fastas/missing
     mkdir databases/RCSB_PDB/cm_lib/msa_fastas/hmm_same
     cp databases/RCSB_PDB/cm_lib/msa_fastas/msas_same/*/*.hmm databases/RCSB_PDB/cm_lib/msa_fastas/hmm_same

     python scripts/P5_copy_same_msas_by_samelist.py --samelist databases/RCSB_PDB/cm_lib/msa_fastas/same_seq_rep.list --hmmdir databases/RCSB_PDB/cm_lib/msa_fastas/hmm_same --outdir databases/RCSB_PDB/msa/hmm/ --fastadir databases/RCSB_PDB/cm_lib/msa_fastas/fasta --missing databases/RCSB_PDB/cm_lib/msa_fastas/missing

     # update fr database
     python scripts/P5_update_fr_slim_v2.py --option_file bin/db_option_2024

     # update template database
     sh databases/hhsuite_pdb/add_tss_hmm.sh 
     ```

   - pdb_complex
     ```
     # Need to change the pdb date cut off in bin/complex_db_option_2023
     # e.g., set_pdb_release_date = 01-05-2024

     python scripts/C1_update_complex_pdb_slim.py --option_file bin/complex_db_option_2024
     python scripts/C2_update_complex_mmcif.py --option_file bin/complex_db_option_2024
     python scripts/C3_update_complex_dataset.py --option_file bin/complex_db_option_2024
     python scripts/C4_update_complex_cm.py --option_file bin/complex_db_option_2024

     # Note: similarly, we need to generate a3m/hmm files for each sequences in the database using the latest uniref30 database
     # We use this script to filter the identical sequences and generate the a3m/hmm files for each sequence (~650000 to ~240000)
     python scripts/C5_generate_msas_fastas.py --fastafile databases/Complex/cm_lib/pdb_cm.fasta --mmseq tools/mmseqs/bin/mmseqs --outputdir databases/Complex/cm_lib/msa_fastas/

     # However, we may have generated the a3m/hmm files for some sequences if we have updated the RCSB_PDB database
     python scripts/C5_copy_same_msas_in_pdb.py --samelist databases/Complex/cm_lib/msa_fastas/same_seq_rep.list --srcfastadir databases/RCSB_PDB/cm_lib/msa_fastas/fasta_rep --trgfastadir databases/Complex/cm_lib/msa_fastas/fasta/ --srca3mdir databases/RCSB_PDB/cm_lib/msa_fastas/hmm_same --outdir databases/Complex/msa/hmm --trgmissing databases/Complex/cm_lib/msa_fastas/missing

     # If missing hmm for some sequences, use the scripts above to generate the hmm files

     # update fr database
     python scripts/C5_update_complex_fr_slim_v2.py --option_file bin/complex_db_option_2024

     # update template database
     sh databases/hhsuite_complex/add_tss_hmm.sh 

     ```

-  AlphaFoldDB template database (https://github.com/google-deepmind/alphafold/tree/main/afdb) ~ 23TB?

-  ColabFold Database (2021_08) (https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh)
     
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