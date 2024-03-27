# MULTICOM4
This is the running note for the tools

- CombFold (https://github.com/dina-lab3D/CombFold)
     1. Extract all pair info from the input fasta (user-defined json file)

     2. Predict the structure of each pair (may need some coding to fit our prediction system)

     3. Rank the pair structure based on the interface plddt scores

     3. Predict the structure of larger groups, where the groups are built by the pairs ordered by their interface plddt scores.

     4. Given the predicted pair/group structures, select the representatives of each subunit by the plddt score in each predicted complex structures, and the transformation information between the subunits in the predicted structures. 


     - Defining subunits

          - Subunit is defined by 4 fields:

               - name: a unique name for the subunit
               - sequence: the amino acid sequence of the subunit
               - chain_names: a list of chain names representing also the stoichiometry of the subunit
               - start_res: the index of the start residue of the sequence on the chain. Needed to set constraints on other subunits on the same chains.

               ```
                    {
                         "A0": {"name": "A0", "chain_names": ["A", "B"], "start_res": 1, "sequence": "MKDILEKLEERRAQARLGGGEKRLEAQHKRGKLTARERIELLLDHGSFEE"},
                         "C0": {"name": "C0", "chain_names": ["C", "D"], "start_res": 1, "sequence": "MFDKILIANRGEIACRIIKTAQKMGIKTVAVYSDADRDAVHVAMADEAVH"},
                         "E0": {"name": "E0", "chain_names": ["E"], "start_res": 1, "sequence": "MGDKIESKKAAAAAEVSTVPGFLGVIESPEHAVTIADEIGYPVMIKASAGA"},
                         "E1": {"name": "E1", "chain_names": ["E"], "start_res": 51, "sequence": "GGGKGMRIAESADEVAEGFARAKSEASSSFGDDRVFVEKFITDPRHIEIQ"},
                    }
               ```

     - Predicting structures for pairs

          ``` 
          python3 scripts/prepare_fastas.py subunits.json --stage pairs --output-fasta-folder <path_to_output_folder> --max-af-size 1800

          # generate models using default_multimer in MULTICOM4

          ```

     - Predicting structures for larger groups

          ```
          python3 scripts/prepare_fastas.py subunits.json  --stage groups --output-fasta-folder <path_to_output_folder>--max-af-size 1800 --input-pairs-results <path_to_AFM_pairs_results>


          ```

     - Combinatorial Assembly

          ```

          python3 scripts/run_on_pdbs.py <path_to_subunits.json> <path_to_folder_of_pdbs> <path_to_empty_output_folder>


          ```