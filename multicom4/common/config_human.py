"""Model config."""

import copy
import ml_collections

MONOMER_PREDICTIONS_PER_MODEL_HUMAN = 200

MONOMER_HUMAN_CONFIG = ml_collections.ConfigDict({
    'common_config': {
        'num_ensemble': 1,
        'num_recycle': 12,
        'predictions_per_model': MONOMER_PREDICTIONS_PER_MODEL_HUMAN,
        'model_preset': 'monomer',
        'relax_topn_predictions': 5,
        'dropout': False,
        'dropout_structure_module': True,
        'msa_source': 'default',
        'template_source': 'pdb70',
        'model_ckpt': None,
    },
    'predictors':{
        'default': {
        },
        'default_pdb70_new': {
            'template_source': 'pdb70_newest',
        },
        'default_seq_temp': {
            'template_source': 'pdb_sort90'
        },
        'original': {
            'msa_source': 'original',
        },
        'ori_seq_temp': {
            'msa_source': 'original',
            'template_source': 'pdb_sort90'
        },
        # 'colabfold': {
        #     'msa_source': 'colabfold',
        #     #'predictions_per_model': 20,
        # },
        # 'colab_seq_temp': {
        #     'msa_source': 'colabfold',
        #     'template_source': 'pdb_sort90',
        #     #'predictions_per_model': 20,
        # },
        # 'img': {
        #     'msa_source': 'img',
        #     #'predictions_per_model': 20,
        # },
        # 'img_seq_temp': {
        #     'msa_source': 'img',
        #     'template_source': 'pdb_sort90',
        #     #'predictions_per_model': 20,
        # },
        'dhr': {
            'msa_source': 'dhr',
        },
        'def_drop_s': {
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_drop_nos': {
            'dropout': True,
            'dropout_structure_module': False,
        },
        'def_notemp': {
            'template_source': 'notemplate',
        },
        'def_notemp_drop_s': {
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_notemp_drop_nos': {
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': False,
        },

        'default_ptm': {
            'model_preset': 'monomer_ptm',
        },
        'def_ptm_drop_s': {
            'model_preset': 'monomer_ptm',
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_ptm_drop_nos': {
            'model_preset': 'monomer_ptm',
            'dropout': True,
            'dropout_structure_module': False,
        },
        'def_ptm_notemp': {
            'model_preset': 'monomer_ptm',
            'template_source': 'notemplate',
        },
        'def_ptm_not_drop_s': {
            'model_preset': 'monomer_ptm',
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_ptm_not_drop_nos': {
            'model_preset': 'monomer_ptm',
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': False,
        },
        
        'deepmsa_dMSA_hhb':{
            'msa_source': 'dMSA.hhb',
        },
        'deepmsa_dMSA_jac':{
            'msa_source': 'dMSA.jac',
        },
        'deepmsa_dMSA_hms':{
            'msa_source': 'dMSA.hms',
        },
        'deepmsa_dMSA':{
            'msa_source': 'dMSA',
        },
        'deepmsa_qMSA':{
            'msa_source': 'qMSA',
        },
        'deepmsa_aMSA':{
            'msa_source': 'aMSA',
        },
        'deepmsa_qMSA_hhb':{
            'msa_source': 'qMSA.hhb',
        },
        'deepmsa_qMSA_jac':{
            'msa_source': 'qMSA.jac',
        },
        'deepmsa_qMSA_hh3':{
            'msa_source': 'qMSA.hh3',
        },
        'deepmsa_qMSA_hms':{
            'msa_source': 'qMSA.hms',
        },
        'deepmsa_DeepJGI_hms':{
            'msa_source': 'DeepJGI.hms',
        },
        'deepmsa_DeepJGI':{
            'msa_source': 'DeepJGI',
        },
        'deepmsa_q3JGI':{
            'msa_source': 'q3JGI',
        },
        'deepmsa_q4JGI':{
            'msa_source': 'q4JGI',
        },
        'deepmsa_q3JGI_hms':{
            'msa_source': 'q3JGI.hms',
        },
        'deepmsa_q4JGI_hms':{
            'msa_source': 'q4JGI.hms',
        },
        # 'default_tmsearch': {
        #     'template_source': 'tmsearch',
        # },
        'def_esm_msa': {
            'input_msa_source': 'default',
            'msa_source': 'esm_msa',
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'def_esm_msa_ckpt5': {
            'input_msa_source': 'default',
            'msa_source': 'esm_msa',
            'model_ckpt': 'model_5',
            'predictions_per_model': MONOMER_PREDICTIONS_PER_MODEL_HUMAN * 5,
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'dom_hhsearch':{
            'msa_source': 'dom_hhsearch',
        },
        'dom_parser':{
            'msa_source': 'dom_parser',
        },
        'dom_unidoc':{
            'msa_source': 'dom_unidoc',
        },
        'dom_manual':{
            'msa_source': 'dom_manual',
        },
        # 'foldseek_refine': {
        #     'number_of_input_models': 5,
        #     'max_iteration': 5,
        #     'max_template_count': 50,
        #     'progressive_threshold': 2000,
        #     'number_of_output_models': 5,
        #     'relax_topn_predictions': 1,
        #     'foldseek_database': 'pdb+afdb',
        #     'msa_source': 'foldseek',
        #     'template_source': 'foldseek',
        #     'predictions_per_model': 5,
        # },
        # 'foldseek_refine_esm': {
        #     'number_of_input_models': 5,
        #     'max_iteration': 5,
        #     'max_template_count': 50,
        #     'progressive_threshold': 2000,
        #     'number_of_output_models': 5,
        #     'relax_topn_predictions': 1,
        #     'foldseek_database': 'esm_atlas',
        #     'plddt_threshold': 0.0,
        #     'ptm_threshold': 0.0,
        #     'msa_source': 'foldseek',
        #     'template_source': 'default',
        #     'predictions_per_model': 5,
        # },
        # 'foldseek_refine_esm_h': {
        #     'number_of_input_models': 5,
        #     'max_iteration': 5,
        #     'max_template_count': 50,
        #     'progressive_threshold': 2000,
        #     'number_of_output_models': 5,
        #     'relax_topn_predictions': 1,
        #     'foldseek_database': 'esm_atlas',
        #     'msa_source': 'foldseek',
        #     'template_source': 'default',
        #     'plddt_threshold': 0.7,
        #     'ptm_threshold': 0.7,
        #     'predictions_per_model': 5,
        # },
    }
})

HETEROMULTIMER_HUMAN_CONFIG = ml_collections.ConfigDict({
    'common_config': {
        'num_ensemble': 1,
        'num_recycle': 20,
        'predictions_per_model': 40,
        'model_preset': 'multimer',
        'relax_topn_predictions': 5,
        'dropout': False,
        'dropout_structure_module': True,
        'msa_unpaired_source': 'default',
        'msa_paired_source': 'default',
        'template_source': 'pdb_seqres',
    },
    'predictors':{
        'default_multimer': {
        },
        'def_mul_struct': {
            'template_source': 'foldseek_structure_based_template'
        },
        'def_mul_tmsearch': {
            'template_source': 'tmsearch_structure_based_template'
        },
        'def_mul_pdb70': {
            'template_source': 'sequence_based_template_pdb70'
        },
        'def_mul_pdb': {
            'template_source': 'sequence_based_template_pdb_sort90'
        },
        'def_mul_comp': {
            'template_source': 'sequence_based_template_pdb_complex'
        },
        # 'def_mul_af': {
        #     'template_source': 'alphafold_model_templates'
        # },
        'def_mul_drop_s': {
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_mul_drop_nos': {
            'dropout': True,
            'dropout_structure_module': False,
        },
        'def_mul_notemp': {
            'template_source': 'notemplate',
        },
        'def_mul_not_drop_s': {
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_mul_not_drop_nos': {
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': False,
        },
        'def_mul_nopair': {
            'msa_paired_source': 'None',
        },
        'uniclust_ox_a3m': {
            'msa_paired_source': 'uniclust_oxmatch_a3m',
        },
        'pdb_inter_ref_a3m': {
            'msa_paired_source': 'pdb_interact_uniref_a3m',
        },
        'pdb_inter_ref_sto': {
            'msa_paired_source': 'pdb_interact_uniref_sto',
        },
        'pdb_inter_prot_sto': {
            'msa_paired_source': 'pdb_interact_uniprot_sto',
        },
        'unidist_ref_a3m': {
            'msa_paired_source': 'uniprot_distance_uniref_a3m',
        },
        'unidist_ref_sto': {
            'msa_paired_source': 'uniprot_distance_uniref_sto',
        },
        'unidist_prot_sto': {
            'msa_paired_source': 'uniprot_distance_uniprot_sto',
        },
        'spec_inter_ref_a3m': {
            'msa_paired_source': 'species_interact_uniref_a3m',
        },
        'spec_struct': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'foldseek_structure_based_template'
        },
        'spec_pdb70': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'sequence_based_template_pdb70'
        },
        'spec_pdb': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'sequence_based_template_pdb_sort90'
        },
        'spec_comp': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'sequence_based_template_pdb_complex'
        },
        # 'spec_af': {
        #     'msa_paired_source': 'species_interact_uniref_a3m',
        #     'template_source': 'alphafold_model_templates'
        # },
        'spec_inter_ref_sto': {
            'msa_paired_source': 'species_interact_uniref_sto',
        },
        'spec_inter_prot_sto': {
            'msa_paired_source': 'species_interact_uniprot_sto',
        },
        'str_inter_ref_a3m': {
            'msa_paired_source': 'string_interact_uniref_a3m',
        },
        'str_inter_ref_sto': {
            'msa_paired_source': 'string_interact_uniref_sto',
        },
        'str_struct': {
            'msa_paired_source': 'string_interact_uniref_sto',
            'template_source': 'foldseek_structure_based_template'
        },
        'str_pdb70': {
            'msa_paired_source': 'string_interact_uniref_sto',
            'template_source': 'sequence_based_template_pdb70'
        },
        'str_pdb': {
            'msa_paired_source': 'string_interact_uniref_sto',
            'template_source': 'sequence_based_template_pdb_sort90'
        },
        'str_comp': {
            'msa_paired_source': 'string_interact_uniref_sto',
            'template_source': 'sequence_based_template_pdb_complex',
        },
        # 'str_af': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        #     'template_source': 'alphafold_model_templates',
        # },
        'str_inter_prot_sto': {
            'msa_paired_source': 'string_interact_uniprot_sto',
        },
        'deepmsa2': {   # common parameters for all deepmsa2 predictors
            'max_pairs': 20,
        },
        'def_mul_esm_msa': {
            'input_msa_source': 'default',
            'msa_paired_source': 'esm_msa',
        },
        'folds_iter': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'foldseek',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_nop': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'None',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'foldseek',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_not': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_notp': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'None',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_esm': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'pdb_seqres',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'folds_iter_esm_nop': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'None',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'pdb_seqres',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'folds_iter_esm_not': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'folds_iter_esm_notp': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'None',
            'msa_unpaired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        # 'AFProfile': 
        # {
        #     'confidence_threshold': 0.95,
        #     'max_iteration': 1000,
        #     'learning_rate': 0.0001,
        # },
        'def_mul_refine': {
            'number_of_input_models': 5,
            'max_iteration': 1,
            'max_template_count': 50,
            'progressive_threshold': 2000,
            'number_of_output_models': 5,
            'relax_topn_predictions': 1,
            'foldseek_database': 'pdb+afdb',
            'msa_source': 'foldseek',
            'template_source': 'foldseek',
        },
        'afsample_v1': {
            'model_preset': "multimer_v1",
            'dropout': True,
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v1_not': {
            'model_preset': "multimer_v1",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v1_r21_not': {
            'model_preset': "multimer_v1",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 21,
        },
        'afsample_v2': {
            'model_preset': "multimer_v2",
            'dropout': True,
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v2_not': {
            'model_preset': "multimer_v2",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v2_r21_not': {
            'model_preset': "multimer_v2",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 21,
        },
        'colabfold_casp16_web': {
            'msa_paired_source': 'colabfold_web',
            'msa_paired_source': 'colabfold_web',
        },
        'colabfold_casp16_web_not': {
            'msa_paired_source': 'colabfold_web',
            'msa_paired_source': 'colabfold_web',
            'template_source': 'notemplate',
        },
    }
})

HOMOMULTIMER_HUMAN_CONFIG = ml_collections.ConfigDict({
    'common_config': {
        'num_ensemble': 1,
        'num_recycle': 20,
        'predictions_per_model': 40,
        'model_preset': 'multimer',
        'relax_topn_predictions': 5,
        'dropout': False,
        'dropout_structure_module': True,
        'msa_paired_source': 'default',
        'template_source': 'pdb_seqres',
    },
    'predictors':{
        'default_multimer': {
        },
        'def_mul_struct': {
            'template_source': 'foldseek_structure_based_template'
        },
        'def_mul_tmsearch': {
            'template_source': 'tmsearch_structure_based_template'
        },
        'def_mul_pdb70': {
            'template_source': 'sequence_based_template_pdb70'
        },
        'def_mul_pdb': {
            'template_source': 'sequence_based_template_pdb_sort90'
        },
        'def_mul_comp': {
            'template_source': 'sequence_based_template_pdb_complex'
        },
        'def_mul_af': {
            'template_source': 'alphafold_model_templates'
        },
        'def_mul_drop_s': {
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_mul_drop_nos': {
            'dropout': True,
            'dropout_structure_module': False,
        },
        'def_mul_notemp': {
            'template_source': 'notemplate',
        },
        'def_mul_not_drop_s': {
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': True,
        },
        'def_mul_not_drop_nos': {
            'template_source': 'notemplate',
            'dropout': True,
            'dropout_structure_module': False,
        },
        'uniclust_ox_a3m': {
            'msa_paired_source': 'uniclust_oxmatch_a3m',
        },
        'pdb_inter_ref_a3m': {
            'msa_paired_source': 'pdb_interact_uniref_a3m',
        },
        'pdb_inter_ref_sto': {
            'msa_paired_source': 'pdb_interact_uniref_sto',
        },
        'pdb_inter_prot_sto': {
            'msa_paired_source': 'pdb_interact_uniprot_sto',
        },
        # 'uniprot_distance_uniref_a3m': {
        #     'msa_paired_source': 'uniprot_distance_uniref_a3m',
        # },
        # 'uniprot_distance_uniref_sto': {
        #     'msa_paired_source': 'uniprot_distance_uniref_sto',
        # },
        # 'uniprot_distance_uniprot_sto': {
        #     'msa_paired_source': 'uniprot_distance_uniprot_sto',
        # },
        'spec_inter_ref_a3m': {
            'msa_paired_source': 'species_interact_uniref_a3m',
        },
        'spec_struct': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'foldseek_structure_based_template'
        },
        'spec_pdb70': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'sequence_based_template_pdb70'
        },
        'spec_pdb': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'sequence_based_template_pdb_sort90'
        },
        'spec_comp': {
            'msa_paired_source': 'species_interact_uniref_a3m',
            'template_source': 'sequence_based_template_pdb_complex'
        },
        # 'spec_af': {
        #     'msa_paired_source': 'species_interact_uniref_a3m',
        #     'template_source': 'alphafold_model_templates'
        # },
        'spec_inter_ref_sto': {
            'msa_paired_source': 'species_interact_uniref_sto',
        },
        'spec_inter_prot_sto': {
            'msa_paired_source': 'species_interact_uniprot_sto',
        },
        # 'string_interact_uniref_a3m': {
        #     'msa_paired_source': 'string_interact_uniref_a3m',
        # },
        # 'string_interact_uniref_sto': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        # },
        # 'str_struct': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        #     'template_source': 'foldseek_structure_based_template'
        # },
        # 'str_pdb70': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        #     'template_source': 'sequence_based_template_pdb70'
        # },
        # 'str_pdb': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        #     'template_source': 'sequence_based_template_pdb_sort90'
        # },
        # 'str_comp': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        #     'template_source': 'sequence_based_template_pdb_complex',
        # },
        # 'str_af': {
        #     'msa_paired_source': 'string_interact_uniref_sto',
        #     'template_source': 'alphafold_model_templates',
        # },
        # 'string_interact_uniprot_sto': {
        #     'msa_paired_source': 'string_interact_uniprot_sto',
        # },
        'deepmsa2': {   # common parameters for all deepmsa2 predictors
            'max_pairs': 50,
        },
        'def_mul_esm_msa': {
            'input_msa_source': 'default',
            'msa_paired_source': 'esm_msa',
        },
        # 'AFProfile': 
        # {
        #     'confidence_threshold': 0.95,
        #     'max_iteration': 1000,
        #     'learning_rate': 0.0001,
        # },
        'folds_iter': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'foldseek',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_not': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_o': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'foldseek',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_o_not': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "pdb+afdb",
            'number_of_input': 2,
        },
        'folds_iter_esm': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'pdb_seqres',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'folds_iter_esm_not': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'folds_iter_esm_o': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'pdb_seqres',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'folds_iter_esm_o_not': {   # common parameters for all deepmsa2 predictors
            'msa_paired_source': 'foldseek',
            'template_source': 'notemplate',
            'foldseek_database': "esm_atlas",
            'number_of_input': 2,
        },
        'def_mul_refine': {
            'number_of_input_models': 5,
            'max_iteration': 5,
            'max_template_count': 50,
            'progressive_threshold': 2000,
            'number_of_output_models': 5,
            'relax_topn_predictions': 1,
            'foldseek_database': 'pdb+afdb',
            'msa_source': 'foldseek',
            'template_source': 'foldseek',
            'predictions_per_model': 40,
        },
        'afsample_v1': {
            'model_preset': "multimer_v1",
            'dropout': True,
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v1_not': {
            'model_preset': "multimer_v1",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v1_r21_not': {
            'model_preset': "multimer_v1",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 21,
        },
        'afsample_v2': {
            'model_preset': "multimer_v2",
            'dropout': True,
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v2_not': {
            'model_preset': "multimer_v2",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'afsample_v2_r21_not': {
            'model_preset': "multimer_v2",
            'dropout': True,
            'dropout_structure_module': False,
            'template_source': 'notemplate',
            'num_ensemble': 1,
            'num_recycle': 21,
        },
        'colabfold_casp16_web': {
            'msa_paired_source': 'colabfold_web',
        },
        'colabfold_casp16_web_not': {
            'msa_paired_source': 'colabfold_web',
            'template_source': 'notemplate',
        },
    }
})