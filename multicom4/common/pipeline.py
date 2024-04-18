import os, sys, argparse, time, copy, math
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.monomer_alignment_generation.alignment import read_fasta, write_fasta
from multicom4.monomer_alignment_generation.pipeline import *
from multicom4.monomer_structure_generation.pipeline_v2 import *
from multicom4.monomer_structure_evaluation.pipeline_sep import *
from multicom4.monomer_templates_search.sequence_based_pipeline_pdb import *
from multicom4.monomer_structure_refinement import iterative_refine_pipeline
from multicom4.multimer_structure_refinement import iterative_refine_pipeline_multimer
from multicom4.monomer_alignments_concatenation.pipeline_v3 import *
from multicom4.monomer_templates_concatenation import sequence_based_pipeline_complex_pdb, \
    sequence_based_pipeline_pdb, sequence_based_pipeline, structure_based_pipeline_v2, \
    structure_tmsearch_pipeline_v2
    
# from multicom4.multimer_structure_generation.pipeline import *
from multicom4.multimer_structure_generation.pipeline_v2 import *
from multicom4.multimer_structure_generation.pipeline_default import *
from multicom4.multimer_structure_generation.pipeline_homo_v2 import *
from multicom4.multimer_structure_generation.iterative_search_pipeline_v0_2 import *
from multicom4.multimer_structure_generation.iterative_search_pipeline_v0_2_old import *
from multicom4.multimer_structure_evaluation.pipeline import *
from multicom4.common.protein import *
import pandas as pd
import numpy as np
from multicom4.monomer_structure_refinement.util import cal_tmscore
from multicom4.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from multicom4.common import config

NUM_FINAL_MODELS = 10

def run_monomer_msa_pipeline(fasta, outdir, params, only_monomer=False):
    uniref30 = params['uniref_db']
    uniclust30 = params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']

    smallbfd = ""  # params['smallbfd_database']
    bfd = params['bfd_database']
    mgnify = params['mgnify_database']

    hhblits_binary = params['hhblits_program']
    hhfilter_binary = params['hhfilter_program']
    jackhmmer_binary = params['jackhmmer_program']

    colabfold_search_binary = params['colabfold_search_program']
    colabfold_split_msas_binary = params['colabfold_split_msas_program']
    colabfold_databases = params['colabfold_databases']
    mmseq_binary = params['mmseq_program']

    deepmsa2_path = params['deepmsa2_path']
    JGIclust_database_path = params['JGIclust_database']
    metaclust_database_path = params['metaclust_database']
    
    dhr_program_path = params['DHR_program_path']
    dhr_database_path = params['DHR_database_path']

    if only_monomer:
        uniprot_fasta = ""
    else:
        uniprot_fasta = params['uniprot_fasta']

    result = None
    try:
        pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary_path=jackhmmer_binary,
                                                         hhblits_binary_path=hhblits_binary,
                                                         hhfilter_binary_path=hhfilter_binary,
                                                         colabfold_search_binary=colabfold_search_binary,
                                                         colabfold_split_msas_binary=colabfold_split_msas_binary,
                                                         mmseq_binary=mmseq_binary,
                                                         deepmsa2_path=deepmsa2_path,
                                                         dhr_program_path=dhr_program_path,
                                                         uniref90_database_path=uniref90_fasta,
                                                         mgnify_database_path=mgnify,
                                                         small_bfd_database_path=smallbfd,
                                                         bfd_database_path=bfd,
                                                         uniref30_database_path=uniref30,
                                                         uniclust30_database_path=uniclust30,
                                                         uniprot_database_path=uniprot_fasta,
                                                         JGIclust_database_path=JGIclust_database_path,
                                                         metaclust_database_path=metaclust_database_path,
                                                         colabfold_databases=colabfold_databases,
                                                         dhr_database_path=dhr_database_path)
        result = pipeline.process(fasta, outdir)
    except Exception as e:
        print(e)
        return result

    # Run disorder prediction in the background
    fasta_name = os.path.basename(fasta)
    diso_out_file = os.path.join(outdir, fasta_name.replace('.fasta', '.diso'))
    logfile = os.path.join(outdir, 'run_diso.log')
    if not os.path.exists(diso_out_file) and not os.path.exists(logfile):
        print("Start to generate disorder prediction")
        cmd = f"perl {params['disopred_path']} {outdir}/{fasta_name} &> {outdir}/run_diso.log &"
        print(cmd)
        os.system(cmd)

    # Run DNCON4 and DeepDist in the background
    targetname = open(fasta).readlines()[0].rstrip('\n')[1:]
    contact_map_file = os.path.join(outdir, 'dncon4', f'{targetname}.dncon2.rr')
    logfile = os.path.join(outdir, 'dncon4.log')
    if not os.path.exists(contact_map_file) and not os.path.exists(logfile):
        print("Start to generate contact prediction using DNCON4")
        os.makedirs(os.path.join(outdir, 'dncon4'), exist_ok=True)
        cmd = f"sh {params['dncon4_program']} {fasta} {outdir}/dncon4 &> {logfile} &"
        print(cmd)
        os.system(cmd)

    dist_map_file = os.path.join(outdir, 'deepdist', f'{targetname}.txt')
    logfile = os.path.join(outdir, 'deepdist.log')
    if not os.path.exists(dist_map_file) and not os.path.exists(logfile):
        print("Start to generate distance map prediction using DeepDist")
        cmd = f"sh {params['deepdist_program']} {fasta} {outdir}/deepdist &> {outdir}/deepdist.log &"
        print(cmd)
        os.system(cmd)

    return result



def copy_same_sequence_msas(srcdir, trgdir, srcname, trgname, rename_prefix=True):
    for msa in os.listdir(srcdir):
        if rename_prefix and msa[0] != srcname:
            continue

        if msa.find('.a3m') > 0 or msa.find('.fasta') > 0:
            contents = open(os.path.join(srcdir, msa))
            new_contents = []
            for i, line in enumerate(contents):
                if i == 0:
                    new_contents += [f">{trgname}\n"]
                else:
                    new_contents += [line]
            
            outfile = os.path.join(trgdir, f"{trgname}{msa[1:]}")
            if not rename_prefix:
                outfile = os.path.join(trgdir, msa)
            fw = open(outfile, 'w')
            fw.writelines(new_contents)

        elif msa.find('.sto') > 0:
            new_contents = []
            for line in open(os.path.join(srcdir, msa)):
                if line[:7] == "#=GF ID":
                    line = line.replace(srcname, trgname)

                tmp = line.split()
                if len(tmp) > 0 and tmp[0] == srcname:
                    line = line[0].replace(srcname, trgname) + line[1:]
                new_contents += [line]
            
            outfile = os.path.join(trgdir, f"{trgname}{msa[1:]}")
            if not rename_prefix:
                outfile = os.path.join(trgdir, msa)
            fw = open(outfile, 'w')
            fw.writelines(new_contents)


def run_monomer_msa_pipeline_img(fasta, outdir, params):
    monomer_msa_pipeline_img = Monomer_alignment_generation_pipeline_img(
        deepmsa_binary_path=params['deepmsa2_program'],
        bfd_database_path=params['bfd_database'],
        img_database_path=params['img_database'],
        metaclust_database_path=params['metaclust_database'],
        mgnify_database_path=params['mgnify_database'],
        uniref90_database_path=params['uniref90_fasta'])

    img_msa = None
    try:
        img_msa = monomer_msa_pipeline_img.process(fasta, outdir)
    except Exception as e:
        print(e)
    return img_msa


def run_monomer_template_search_pipeline(params, targetname, sequence, a3m, outdir):
    template_file = None
    try:
        pipeline = monomer_sequence_based_template_search_pipeline(params)
        template_file = pipeline.search(targetname=targetname, sequence=sequence, a3m=a3m, outdir=outdir)
    except Exception as e:
        print(e)
    return template_file

def run_monomer_structure_generation_pipeline_v2(params, targetname, fasta_path, alndir, img_alndir, 
                                                 templatedir, outdir, run_methods=None, 
                                                 run_script=False):
    try:
        pipeline = Monomer_structure_prediction_pipeline_v2(params, run_methods=run_methods)
        pipeline.process_single(targetname=targetname,
                                fasta_path=fasta_path,
                                alndir=alndir,
                                img_alndir=img_alndir,
                                template_dir=templatedir,
                                outdir=outdir,
                                run_script=run_script)
    except Exception as e:
        print(e)
        return False
    return True

def cal_tmscores(tmscore_program, src_pdbs, trg_pdb, outputdir):
    tmscores = []
    for src_pdb in src_pdbs:
        tmscore, gdtscore = cal_tmscore(tmscore_program,
                                        src_pdb,
                                        trg_pdb,
                                        os.path.join(outputdir, 'tmp'))
        tmscores += [tmscore]
    return tmscores

def select_models_by_tmscore(tmscore_program, ranking_df_file, outputdir, prefix, tmscore_threshold):
    selected_models = []
    selected_models_path = []
    ranking_df = pd.read_csv(ranking_df_file)
    for i in range(30):
        model = ranking_df.loc[i, 'model']
        tmscores = cal_tmscores(tmscore_program, selected_models_path, os.path.join(outputdir, 'pdb', model), outputdir)
        if len(tmscores) == 0 or np.max(np.array(tmscores)) < tmscore_threshold:
            selected_models += [model]
            selected_models_path += [os.path.join(outputdir, 'pdb', model)]

    for i in range(30):
        model = ranking_df.loc[i, 'model']
        if model in selected_models:
            continue
        selected_models += [model]
        selected_models_path += [os.path.join(outputdir, 'pdb', model)]
        if len(selected_models) >= NUM_FINAL_MODELS:
            break

    for i in range(NUM_FINAL_MODELS):
        final_pdb = os.path.join(outputdir, f'{prefix}{i+1}.pdb')
        os.system("cp " + selected_models_path[i] + " " + final_pdb)
    
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(os.path.join(outputdir, f'{prefix}_selected.csv'))

def select_models_by_cluster(ranking_df_file, cluster_result_file, outputdir, prefix):

    clusters = {}
    for line in open(cluster_result_file):
        modelname, cluster_idx = line.rstrip('\n').split()
        clusters[modelname + '.pdb'] = cluster_idx
    
    selected_models = []
    selected_clusters = []
    ranking_df = pd.read_csv(ranking_df_file)
    for i in range(30):
        model = ranking_df.loc[i, 'model']
        if clusters[model] not in selected_clusters:
            selected_models += [model]
            selected_clusters += [clusters[model]]

    for i in range(30):
        model = ranking_df.loc[i, 'model']
        if model in selected_models:
            continue
        selected_models += [model]
        if len(selected_models) >= NUM_FINAL_MODELS:
            break

    for i in range(NUM_FINAL_MODELS):
        final_pdb = os.path.join(outputdir, f'{prefix}{i+1}.pdb')
        os.system("cp " + os.path.join(outputdir, 'pdb', selected_models[i]) + " " + final_pdb)
    
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(os.path.join(outputdir, f'{prefix}_selected.csv'))
    return selected_models

def select_models_monomer_only(qa_result, outputdir, params):

    select_models_by_tmscore(tmscore_program=params['tmscore_program'],
                             ranking_df_file=qa_result["alphafold"],
                             outputdir=outputdir,
                             prefix="ai",
                             tmscore_threshold=0.8)
    
    if "gate" not in qa_result:
        return

    select_models_by_cluster(ranking_df_file=qa_result["gate"],
                             cluster_result_file=qa_result["gate_cluster"],
                             outputdir=outputdir,
                             prefix="gate")
    
    select_models_by_cluster(ranking_df_file=qa_result["gate_af_avg"],
                             cluster_result_file=qa_result["gate_cluster"],
                             outputdir=outputdir,
                             prefix="llm")


def select_models_with_multimer(qa_result, outputdir, params):
    
    select_models_by_tmscore(tmscore_program=params['tmscore_program'],
                             ranking_df_file=qa_result["alphafold_multimer"],
                             outputdir=outputdir,
                             prefix="ai",
                             tmscore_threshold=0.8)

    if "gate" not in qa_result:
        return

    select_models_by_cluster(ranking_df_file=qa_result["gate_multimer"],
                             cluster_result_file=qa_result["gate_cluster"],
                             outputdir=outputdir,
                             prefix="gate")
    
    select_models_by_cluster(ranking_df_file=qa_result["gate_af_avg_multimer"],
                             cluster_result_file=qa_result["gate_cluster"],
                             outputdir=outputdir,
                             prefix="llm")


def run_monomer_evaluation_pipeline(params, targetname, fasta_file, input_monomer_dir, outputdir, 
                                    contact_map_file, dist_map_file,
                                    run_methods=None,
                                    input_multimer_dir="",
                                    generate_final_models=False, model_count=5):

    makedir_if_not_exists(outputdir)
    qa_result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     use_gpu=True,
                                                     run_methods=run_methods)
    try:
        qa_result = pipeline.process(targetname=targetname, fasta_file=fasta_file,
                                     monomer_model_dir=input_monomer_dir, multimer_model_dir=input_multimer_dir,
                                     output_dir=outputdir, model_count=model_count,
                                     contact_map_file=contact_map_file, dist_map_file=dist_map_file)
    except Exception as e:
        print(e)

    if generate_final_models:
        if input_multimer_dir == "" or not os.path.exists(input_multimer_dir):
            select_models_monomer_only(qa_result=qa_result, outputdir=outputdir, params=params)
        else:
            select_models_with_multimer(qa_result=qa_result, outputdir=outputdir, params=params)

    return qa_result


def rerun_monomer_evaluation_pipeline(params, targetname, fasta_file, outputdir):
    makedir_if_not_exists(outputdir)
    result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     use_gpu=True)
    try:
        result = pipeline.reprocess(targetname=targetname, fasta_file=fasta_file,
                                    output_dir=outputdir)
    except Exception as e:
        print(e)

    return result


# def run_monomer_refinement_pipeline(params, refinement_inputs, outdir, finaldir, prefix):
#     pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params)
#     pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir)

#     makedir_if_not_exists(finaldir)

#     pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection(params)
#     pipeline.select_v1(indir=outdir, outdir=finaldir, prefix=prefix)


def getTopN(maxP, chainN, is_homomers):
    topN = 10
    if is_homomers:
        topN = 50
    else:
        topN = math.floor(math.pow(maxP, 1.0 / chainN))
    return topN

def run_monomer_msas_concatenation_pipeline(multimer, run_methods, monomer_aln_dir, monomer_model_dir,
                                            outputdir, params, is_homomers=False):
    chains = multimer.split(',')
    alignment = {'outdir': outputdir}

    config_deepmsa2 = config.HETEROMULTIMER_CONFIG.predictors.deepmsa2
    if is_homomers:
        config_deepmsa2 = config.HOMOMULTIMER_CONFIG.predictors.deepmsa2

    topN = getTopN(config_deepmsa2.max_pairs, len(chains), is_homomers)

    for i in range(len(chains)):
        chain = chains[i]
        chain_aln_dir = os.path.join(monomer_aln_dir, chain)
        if os.path.exists(chain_aln_dir):
            chain_a3ms = {'name': chain,
                          'colabfold_a3m': os.path.join(chain_aln_dir, f"{chain}_colabfold.a3m"),
                          'uniref30_a3m': os.path.join(chain_aln_dir, f"{chain}_uniref30.a3m"),
                          'uniref90_sto': os.path.join(chain_aln_dir, f"{chain}_uniref90.sto"),
                          'uniprot_sto': os.path.join(chain_aln_dir, f"{chain}_uniprot.sto"),
                          'uniclust30_a3m': os.path.join(chain_aln_dir, f"{chain}_uniclust30.a3m")}

            deepmsa_ranking_file = os.path.join(monomer_model_dir, chain, 'deepmsa.rank')
            contents = open(deepmsa_ranking_file).readlines()
            for index, line in enumerate(contents):
                if index >= topN:
                    break
                line = line.rstrip('\n')
                msa_name, plddt = line.split()
                msa_name = msa_name.replace('_', '.')
                chain_a3ms['deepmsa_' + msa_name] = os.path.join(chain_aln_dir, 'DeepMSA2_a3m', 'finalMSAs', msa_name + '.a3m')
            chain_a3ms['deepmsa_ranking_file'] = deepmsa_ranking_file
        else:
            chain_a3ms = {'name': chain}
        alignment[f"chain{i + 1}"] = chain_a3ms

    complete = True
    for name in alignment:
        if name == 'outdir':
            continue
        a3ms = alignment[name]
        if len(a3ms) == 1:
            complete = False
        for key in a3ms:
            if key.find('uni') >= 0:
                if not os.path.exists(a3ms[key]):
                    complete = False
                else:
                    contents = open(a3ms[key]).readlines()
                    if len(contents) == 0:
                        print(f"File: {a3ms[key]} is empty!")
                        complete = False
                        os.system(f"rm {a3ms[key]}")

    if complete:
        alignments = [alignment]
        print(f"Total {len(alignments)} pairs can be concatenated")
        print("Start to concatenate alignments for dimers")

        if not os.path.exists(os.path.join(alignment['outdir'], 'DONE')):
            monomer_alignments_concatenation_pipeline = Monomer_alignments_concatenation_pipeline(params=params,
                                                                                                run_methods=run_methods)
            alignments = monomer_alignments_concatenation_pipeline.concatenate(alignments, is_homomers=is_homomers)
        else:
            print("The multimer alignments have been generated!")
    else:
        print("The a3ms for dimers are not complete!")


def run_monomer_templates_concatenation_pipeline(multimers, monomer_aln_dir, monomer_model_dir, outdir, params):
    monomer_template_inputs = []
    monomer_sequences = []
    for chain in multimers:
        chain_template_a3m = os.path.join(monomer_aln_dir, chain, f"{chain}_uniref90.sto")
        if not os.path.exists(chain_template_a3m):
            raise Exception(f"Cannot find uniref90 alignment for {chain}")
        seq = open(os.path.join(monomer_aln_dir, chain, f"{chain}.fasta")).readlines()[1].rstrip('\n')
        chain_template_input = sequence_based_pipeline.monomer_template_input(name=chain,
                                                                              msa_path=chain_template_a3m,
                                                                              hmm_path="", seq=seq)
        monomer_template_inputs += [chain_template_input]
        monomer_sequences += [seq]

    print("searching complex sequence based template search pipelineL RCSB_PDB")
    pdb_seq_dir = os.path.join(outdir, 'pdb_seq')
    makedir_if_not_exists(pdb_seq_dir)
    if not os.path.exists(os.path.join(pdb_seq_dir, 'sequence_templates.csv')):
        pipeline = sequence_based_pipeline_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(copy.deepcopy(monomer_template_inputs), pdb_seq_dir)

    print("searching complex sequence based template search pipeline: Complex")
    complex_pdb_seq_dir = os.path.join(outdir, 'complex_pdb_seq')
    makedir_if_not_exists(complex_pdb_seq_dir)
    if not os.path.exists(os.path.join(complex_pdb_seq_dir, 'sequence_templates.csv')):
        pipeline = sequence_based_pipeline_complex_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(copy.deepcopy(monomer_template_inputs), complex_pdb_seq_dir)

    print("searching complex sequence based template search pipeline: pdb70")
    pdb70_seq_dir = os.path.join(outdir, 'pdb70_seq')
    makedir_if_not_exists(pdb70_seq_dir)
    if not os.path.exists(os.path.join(pdb70_seq_dir, 'sequence_templates.csv')):
        pipeline = sequence_based_pipeline.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(copy.deepcopy(monomer_template_inputs), pdb70_seq_dir)

    print("searching complex structure based template search pipeline")
    struct_temp_dir = os.path.join(outdir, 'struct_temp')
    makedir_if_not_exists(struct_temp_dir)
    if not os.path.exists(os.path.join(struct_temp_dir, 'structure_templates.csv')):
        monomer_pdbs = []
        for chain in multimers:
            monomer_pdb = os.path.join(monomer_model_dir, chain, 'default', 'ranked_0.pdb')
            if not os.path.exists(monomer_pdb):
                print(f"Cannot find teritary structure for {chain}: {monomer_pdb}")
                continue
            monomer_trg_pdb = os.path.join(struct_temp_dir, f"{chain}.pdb")
            os.system(f"cp {monomer_pdb} {monomer_trg_pdb}")
            monomer_pdbs += [monomer_trg_pdb]
        pipeline = structure_based_pipeline_v2.Complex_structure_based_template_search_pipeline(params)
        pipeline.search(monomer_sequences, monomer_pdbs, struct_temp_dir)

    print("searching complex structure based tmsearch pipeline")
    struct_temp_dir = os.path.join(outdir, 'tmsearch')
    makedir_if_not_exists(struct_temp_dir)
    if not os.path.exists(os.path.join(struct_temp_dir, 'structure_templates.csv')):
        monomer_pdbs = []
        monomer_tmsearch_result_dirs = []
        for chain in multimers:
            monomer_pdb = os.path.join(monomer_model_dir, chain, 'default', 'ranked_0.pdb')
            if not os.path.exists(monomer_pdb):
                print(f"Cannot find teritary structure for {chain}: {monomer_pdb}")
                continue
            monomer_trg_pdb = os.path.join(struct_temp_dir, f"{chain}.pdb")
            os.system(f"cp {monomer_pdb} {monomer_trg_pdb}")
            monomer_pdbs += [monomer_trg_pdb]

            monomer_tmsearch_result_dir = os.path.join(monomer_model_dir, chain, 'default_tmsearch', 'tmsearch')
            monomer_tmsearch_result_dirs += [monomer_tmsearch_result_dir]

        pipeline = structure_tmsearch_pipeline_v2.Complex_structure_tmsearch_based_template_search_pipeline(params)
        pipeline.search(monomer_sequences, monomer_pdbs, monomer_tmsearch_result_dirs, struct_temp_dir)


# def run_multimer_structure_generation_pipeline(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
#                                                  template_dir,
#                                                  monomer_model_dir, output_dir):
#     try:
#         pipeline = multimer_structure_prediction_pipeline(params)
#         result = pipeline.process(fasta_path=fasta_path,
#                                   chain_id_map=chain_id_map,
#                                   aln_dir=aln_dir,
#                                   complex_aln_dir=complex_aln_dir,
#                                   template_dir=template_dir,
#                                   monomer_model_dir=monomer_model_dir,
#                                   output_dir=output_dir)
#     except Exception as e:
#         print(e)
#         return False
#     return True


def run_multimer_structure_generation_pipeline_v2(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
                                                    template_dir,
                                                    monomer_model_dir, 
                                                    output_dir, 
                                                    run_methods=None, 
                                                    run_script=False,
                                                    run_deepmsa=True):
    try:
        pipeline = Multimer_structure_prediction_pipeline_v2(params, run_methods)
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  complex_aln_dir=complex_aln_dir,
                                  template_dir=template_dir,
                                  monomer_model_dir=monomer_model_dir,
                                  output_dir=output_dir,
                                  run_script=run_script)

        if run_deepmsa:
            print(f"Start to generate DeepMSA models")
            result = pipeline.process_deepmsa(fasta_path=fasta_path,
                                            chain_id_map=chain_id_map,
                                            aln_dir=aln_dir,
                                            complex_aln_dir=complex_aln_dir,
                                            output_dir=output_dir,
                                            run_script=run_script)
    except Exception as e:
        print(e)
        return False
    return True


def run_multimer_structure_generation_pipeline_default(params, fasta_path, chain_id_map, aln_dir, output_dir):
    try:
        pipeline = Multimer_structure_prediction_pipeline_default(params)
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  output_dir=output_dir)
    except Exception as e:
        print(e)
        return False
    return True

# def run_multimer_structure_generation_homo_pipeline_img_v2(params, fasta_path, chain_id_map, aln_dir, output_dir):
#     try:
#         pipeline = Multimer_structure_prediction_homo_pipeline_v2(params)
#         result = pipeline.process_img(fasta_path=fasta_path,
#                                       chain_id_map=chain_id_map,
#                                       aln_dir=aln_dir,
#                                       output_dir=output_dir)
#     except Exception as e:
#         print(e)
#         return False
#     return True


def run_multimer_structure_generation_homo_pipeline_v2(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
                                                       template_dir,
                                                       monomer_model_dir, output_dir, 
                                                       run_methods=None, 
                                                       run_script=False,
                                                       run_deepmsa=True):
    try:
        pipeline = Multimer_structure_prediction_homo_pipeline_v2(params, run_methods)
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  complex_aln_dir=complex_aln_dir,
                                  template_dir=template_dir,
                                  monomer_model_dir=monomer_model_dir,
                                  output_dir=output_dir,
                                  run_script=run_script)
        
        if run_deepmsa:
            print(f"Start to generate DeepMSA models")
            result = pipeline.process_deepmsa(fasta_path=fasta_path,
                                            chain_id_map=chain_id_map,
                                            aln_dir=aln_dir,
                                            complex_aln_dir=complex_aln_dir,
                                            output_dir=output_dir,
                                            run_script=run_script)
    except Exception as e:
        print(e)
        return False
    return True


class foldseek_iterative_monomer_input:
    def __init__(self, monomer_pdb_dirs, monomer_alphafold_a3ms):
        self.monomer_pdb_dirs = monomer_pdb_dirs
        self.monomer_alphafold_a3ms = monomer_alphafold_a3ms


def run_multimer_structure_generation_pipeline_foldseek(params, fasta_path, chain_id_map, pipeline_inputs, 
                                                        outdir, config_name, monomer_template_stos=[], start=0, is_homomers=False):
    pipeline = Multimer_iterative_generation_pipeline_monomer(params=params, config_name=config_name, is_homomers=is_homomers)
    try:
        for i, pipeline_input in enumerate(pipeline_inputs):
            if is_homomers:
                pipeline.search_single_homo(fasta_file=fasta_path,
                                            chain_id_map=chain_id_map,
                                            monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                            monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                            monomer_template_stos=monomer_template_stos,
                                            outdir=os.path.join(outdir, f"{config_name}_{i + 1 + start}"))
            else:
                pipeline.search_single(fasta_file=fasta_path,
                                       chain_id_map=chain_id_map,
                                       monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                       monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                       monomer_template_stos=monomer_template_stos,
                                       outdir=os.path.join(outdir, f"{config_name}_{i + 1 + start}"))
    except Exception as e:
        print(e)
        return False
    return True


def run_multimer_structure_generation_pipeline_foldseek_old(params, fasta_path, chain_id_map, pipeline_inputs, outdir,
                                                            config_name, monomer_template_stos=[], start=0, is_homomers=False):
    pipeline = Multimer_iterative_generation_pipeline_monomer_old(params=params, config_name=config_name, is_homomers=is_homomers)
    try:
        for i, pipeline_input in enumerate(pipeline_inputs):
            if is_homomers:
                pipeline.search_single_homo(fasta_file=fasta_path,
                                            chain_id_map=chain_id_map,
                                            monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                            monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                            monomer_template_stos=monomer_template_stos,
                                            outdir=os.path.join(outdir, f"{config_name}_{i + 1 + start}"))
            else:
                pipeline.search_single(fasta_file=fasta_path,
                                       chain_id_map=chain_id_map,
                                       monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                       monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                       monomer_template_stos=monomer_template_stos,
                                       outdir=os.path.join(outdir, f"{config_name}_{i + 1 + start}"))
    except Exception as e:
        print(e)
        return False
    return True


def extract_monomer_models_from_complex(complex_pdb, complex_pkl, chain_id_map, workdir):
    makedir_if_not_exists(workdir)

    chain_group = {}
    for chain_id in chain_id_map:
        find = False
        for chain_id_seq in chain_group:
            if chain_id_map[chain_id_seq].sequence == chain_id_map[chain_id].sequence:
                chain_group[chain_id_seq] += [chain_id]
                find = True
        if not find:
            chain_group[chain_id] = [chain_id]

    pkldir = os.path.join(workdir, 'pkl')
    makedir_if_not_exists(pkldir)

    alphafold_qa = Alphafold_pkl_qa(ranking_methods=['plddt_avg'])
    for chain_id in chain_group:

        chain_pdb_dict = extract_monomer_pdbs(complex_pdb=complex_pdb,
                                              sequence=chain_id_map[chain_id].sequence,
                                              output_prefix=os.path.join(workdir, 'chain'))
        print(chain_pdb_dict)

        same_seq_chains = chain_group[chain_id]
        same_seq_pkldir = os.path.join(pkldir, chain_id)
        makedir_if_not_exists(same_seq_pkldir)
        for same_seq_chain in same_seq_chains:
            pdbname = chain_pdb_dict[same_seq_chain]['pdbname']
            extract_pkl(src_pkl=complex_pkl,
                        residue_start=chain_pdb_dict[same_seq_chain]['chain_start'],
                        residue_end=chain_pdb_dict[same_seq_chain]['chain_end'],
                        output_pkl=os.path.join(same_seq_pkldir, pdbname.replace('.pdb', '.pkl')))
        alphafold_ranking = alphafold_qa.run(same_seq_pkldir)
        alphafold_ranking.to_csv(same_seq_pkldir + '_alphafold_ranking.csv')
        chain_top1_pdb = os.path.join(workdir, chain_id + "_top1.pdb")
        os.system("cp " + os.path.join(workdir, alphafold_ranking.loc[0, 'model']) + " " + chain_top1_pdb)

    return chain_group


def rerun_multimer_evaluation_pipeline(params, fasta_path, chain_id_map, outdir, is_homomer=False):
    makedir_if_not_exists(outdir)
    pipeline = Multimer_structure_evaluation_pipeline(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.reprocess(fasta_path=fasta_path,
                                                chain_id_map=chain_id_map, output_dir=outdir, 
                                                is_homomer=is_homomer)
    except Exception as e:
        print(e)

def cal_tmscores_usalign(usalign_program, src_pdbs, trg_pdb, outputdir):
    tmscores = []
    for src_pdb in src_pdbs:
        cmd = usalign_program + ' ' + src_pdb + ' ' + trg_pdb + " -TMscore 6 -ter 1 | grep TM-score | awk '{print $2}' "
        tmscore_contents = os.popen(cmd).read().split('\n')
        tmscore = float(tmscore_contents[1].rstrip('\n'))
        tmscores += [tmscore]

    return tmscores

def select_models_by_usalign(usalign_program, ranking_df_file, outputdir, prefix, tmscore_threshold):
    selected_models = []
    selected_model_paths = []
    ranking_df = pd.read_csv(ranking_df_file)
    for i in range(30):
        model = ranking_df.loc[i, 'model']
        tmscores = cal_tmscores_usalign(usalign_program, selected_model_paths, os.path.join(outputdir, 'pdb', model), outputdir)
        if len(tmscores) == 0 or np.max(np.array(tmscores)) < tmscore_threshold:
            selected_models += [model]
            selected_model_paths += [os.path.join(outputdir, 'pdb', model)]

    for i in range(30):
        model = ranking_df.loc[i, 'model']
        if model in selected_models:
            continue
        selected_models += [model]
        selected_model_paths += [os.path.join(outputdir, 'pdb', model)]
        if len(selected_models) >= NUM_FINAL_MODELS:
            break

    for i in range(NUM_FINAL_MODELS):
        final_pdb = os.path.join(outputdir, f'{prefix}{i+1}.pdb')
        os.system("cp " + selected_model_paths[i] + " " + final_pdb)
    
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(os.path.join(outputdir, f'{prefix}_selected.csv'))
    return selected_models


def run_multimer_evaluation_pipeline(params, fasta_path, chain_id_map,
                                     indir, outdir, is_homomer=False, model_count=5):
    makedir_if_not_exists(outdir)
    pipeline = Multimer_structure_evaluation_pipeline(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.process(fasta_path=fasta_path,
                                              chain_id_map=chain_id_map,
                                              model_dir=indir,
                                              output_dir=outdir, 
                                              model_count=model_count,
                                              is_homomer=is_homomer)
    except Exception as e:
        print(e)

    select_models_multimer(chain_id_map=chain_id_map, qa_result=multimer_qa_result, outputdir=outdir, params=params, is_homomer=is_homomer)

    return multimer_qa_result

def extract_final_monomer_models_from_complex(chain_id_map, selected_models, outputdir, prefix, is_homomer):
    for i in range(NUM_FINAL_MODELS):
        final_pdb = os.path.join(outputdir, f"{prefix}{i + 1}.pdb")
        if not is_homomer:
            chain_group = extract_monomer_models_from_complex(complex_pdb=final_pdb,
                                                              complex_pkl=os.path.join(outputdir, 'pkl', selected_models[i].replace('.pdb', '.pkl')),
                                                              chain_id_map=chain_id_map, 
                                                              workdir=os.path.join(outputdir, f"{prefix}{i + 1}"))
            for chain_id in chain_group:
                chain_outdir = os.path.join(outputdir, f"{prefix}_{chain_id}")
                makedir_if_not_exists(chain_outdir)
                srcpdb = os.path.join(outputdir, f"{prefix}{i + 1}", f"{chain_id}_top1.pdb")
                trgpdb = os.path.join(chain_outdir, f'{prefix}{i+1}.pdb')
                os.system(f"cp {srcpdb} {trgpdb}")

def select_models_multimer(chain_id_map, qa_result, outputdir, params, is_homomer):
    
    ai_selected_models = select_models_by_usalign(usalign_program=params['usalign_program'],
                                                ranking_df_file=qa_result["alphafold"],
                                                outputdir=outputdir,
                                                prefix="ai",
                                                tmscore_threshold=0.6)

    extract_final_monomer_models_from_complex(chain_id_map=chain_id_map,
                                              selected_models=ai_selected_models, 
                                              outputdir=outputdir, 
                                              prefix="ai", 
                                              is_homomer=is_homomer)

    gate_selected_models = select_models_by_cluster(ranking_df_file=qa_result["gate"],
                             cluster_result_file=qa_result["gate_cluster"],
                             outputdir=outputdir,
                             prefix="gate")
    

    extract_final_monomer_models_from_complex(chain_id_map=chain_id_map,
                                              selected_models=gate_selected_models, 
                                              outputdir=outputdir, 
                                              prefix="gate", 
                                              is_homomer=is_homomer)

    llm_selected_models = select_models_by_cluster(ranking_df_file=qa_result["gate_af_avg"],
                             cluster_result_file=qa_result["gate_cluster"],
                             outputdir=outputdir,
                             prefix="llm")

    extract_final_monomer_models_from_complex(chain_id_map=chain_id_map,
                                              selected_models=llm_selected_models, 
                                              outputdir=outputdir, 
                                              prefix="llm", 
                                              is_homomer=is_homomer)
