# Modified from Alphafold2 codes

"""Library to run HHblits from Python."""

import glob
import os
import subprocess
from typing import Any, Mapping, Optional, Sequence
from absl import logging
from multicom4.tool import utils
import shlex
from concurrent.futures import ProcessPoolExecutor

def run_command(inparams):
    cmdcontent, title = inparams
    logging.info('Launching subprocess "%s"', cmdcontent)
    cmd = shlex.split(cmdcontent)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with utils.timing(title):
        stdout, stderr = process.communicate()
        retcode = process.wait()
    

class DeepMSA2_runner:
    """Python wrapper of the HHblits binary."""

    def __init__(self,
                 tool_path,
                 bfd_database_path,
                 metaclust_database_path,
                 mgnify_database_path,
                 uniref90_database_path,
                 uniref30_database_path,
                 uniclust30_database_path,
                 JGIclust_database_path,
                 ):
        
        self.bfd_database_path = bfd_database_path
        self.metaclust_database_path = metaclust_database_path
        self.mgnify_database_path = mgnify_database_path
        self.uniref90_database_path = uniref90_database_path
        self.uniref30_database_path = uniref30_database_path
        self.uniclust30_database_path = uniclust30_database_path
        self.JGIclust_database_path = JGIclust_database_path

        self.tool_path = tool_path
        self.qMSApkg=os.path.join(tool_path, "bin/qMSA")
        self.dMSApkg=os.path.join(tool_path, "bin/dMSA")

    def query(self, input_fasta_path: str, outpath: str, cpu: int=10) -> Mapping[str, Any]:
        """Queries the database using HHblits."""
        
        targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        # Modified from DeepMSA2_noIMG.py
        dMSAhhblitsdb = self.uniclust30_database_path
        dMSAjackhmmerdb = self.uniref90_database_path
        dMSAhmmsearchdb = self.metaclust_database_path

        qMSAhhblitsdb = self.uniref30_database_path
        qMSAjackhmmerdb = self.uniref90_database_path
        qMSAhhblits3db = self.bfd_database_path
        qMSAhmmsearchdb = self.mgnify_database_path

        mMSAJGI = self.JGIclust_database_path

        # qMSA
        MSAdir = os.path.join(outpath, 'MSA')
        os.makedirs(MSAdir, exist_ok=True)
        qseq_file = os.path.join(MSAdir, 'qMSA.fasta')
        os.system(f"cp {input_fasta_path} {qseq_file}")
        
        process_list = []
        if not os.path.exists(os.path.join(MSAdir, 'qMSA.aln')):
            cmdcontent = f"python {self.qMSApkg}/scripts/qMSA2.py " \
                        f"-hhblitsdb={qMSAhhblitsdb} " \
                        f"-jackhmmerdb={qMSAjackhmmerdb} " \
                        f"-hhblits3db={qMSAhhblits3db} " \
                        f"-hmmsearchdb={qMSAhmmsearchdb} " \
                        f"-ncpu={cpu} " \
                        f"-outdir={MSAdir} " \
                        f"{qseq_file}"
            process_list.append([cmdcontent, 'qMSA query'])
            
        # dMSA
        if not os.path.exists(os.path.join(MSAdir, 'dMSA.aln')):
            dseq_file = os.path.join(MSAdir, 'dMSA.fasta')
            os.system(f"cp {input_fasta_path} {dseq_file}")
            cmdcontent = f"python {self.dMSApkg}/scripts/build_MSA.py " \
                        f"-hhblitsdb={dMSAhhblitsdb} " \
                        f"-jackhmmerdb={dMSAjackhmmerdb} " \
                        f"-hmmsearchdb={dMSAhmmsearchdb} " \
                        f"-ncpu={cpu} " \
                        f"-outdir={MSAdir} " \
                        f"{dseq_file}"        
            process_list.append([cmdcontent, 'dMSA query'])

        if len(process_list) > 0:
            with ProcessPoolExecutor(max_workers=len(process_list)) as executor:
                results = executor.map(run_command, process_list)

        if os.path.isfile(f"{MSAdir}/dMSA.a3m") and os.path.isfile(f"{MSAdir}/qMSA.a3m"):
            print("DeepMSA2_noIMG has been complete!")
       
        JGIdir = os.path.join(outpath, 'JGI')
        os.makedirs(JGIdir, exist_ok=True)
        if not os.path.isfile(f"{MSAdir}/qMSA.jaca3m") and not os.path.isfile(f"{MSAdir}/dMSA.jaca3m"):
            print(f"{targetname} does not need additional JGI search due to no jack result. skip")
        else:

            process_list = []

            recorddir = f"{outpath}/record"
            os.makedirs(recorddir, exist_ok=True)

            logging.info('1st step is starting!')
            print("5555555555555555555555")
            user_id = os.environ["USER"]

            content = open(f"{mMSAJGI}/list").readlines()
            for j in range(len(content)):
                DBfasta = content[j].strip('\n')
                if DBfasta == "":
                    break
                if os.path.exists(f"{JGIdir}/{DBfasta}.cdhit"):
                    continue
                tag = DBfasta
                tmpdir = f"{JGIdir}/tmp/{user_id}/{tag}"
                
                # cmdcontent = f"mkdir -p {tmpdir}\n" \
                # f"cd {tmpdir}\n\n" \
                # f"echo hostname: `hostname`  >>{recorddir}/ware_{tag}\n" \
                # f"echo starting time: `date` >>{recorddir}/ware_{tag}\n" \
                # f"echo pwd `pwd`             >>{recorddir}/ware_{tag}\n\n" \
                # f"cp {MSAdir}/qMSA.hh3aln seq.hh3aln\n" \
                # f"if [ ! -s 'seq.hh3aln' ];then\n" \
                # f"    cp {MSAdir}/dMSA.jacaln seq.hh3aln\n" \
                # f"fi\n" \
                # f"if [ ! -s 'seq.hh3aln' ];then\n" \

                # f"    cp {MSAdir}/dMSA.hhbaln seq.hh3aln\n" \
                # f"fi\n\n" \
                # f"sed = seq.hh3aln |sed 'N;s/\\n/\\t/'|sed 's/^/>/g'|sed 's/\\t/\\n/g'| {self.qMSApkg}/bin/qhmmbuild -n aln --amino -O seq.afq --informat afa seq.hmm -\n\n" \
                # f"{self.qMSApkg}/bin/qhmmsearch --cpu 1 -E 10 --incE 1e-3 -A {DBfasta}.match --tblout {DBfasta}.tbl -o {DBfasta}.out seq.hmm {JGI}/{DBfasta}\n" \
                # f"{self.qMSApkg}/bin/esl-sfetch -f {JGI}/{DBfasta} {DBfasta}.tbl|sed 's/*//g' > {DBfasta}.fseqs\n" \
                # f"{self.qMSApkg}/bin/cd-hit -i {DBfasta}.fseqs -o {datadir}/JGI/{DBfasta}.cdhit -c 1 -M 3000\n\n" \
                # f"sync\n" \
                # f"rm -rf {tmpdir}\n" \
                # f"echo ending time: `date`   >>{recorddir}/ware_{tag}\n"

                bashfile = os.path.join(recorddir, DBfasta + ".sh")
                with open(bashfile, 'w') as fw:
                    fw.write(f"mkdir -p {tmpdir}\n")
                    fw.write(f"cd {tmpdir}\n\n")
                    fw.write(f"echo hostname: `hostname`  >>{recorddir}/ware_{tag}\n")
                    fw.write(f"echo starting time: `date` >>{recorddir}/ware_{tag}\n")
                    fw.write(f"echo pwd `pwd`             >>{recorddir}/ware_{tag}\n\n")
                    fw.write(f"cp {MSAdir}/qMSA.hh3aln seq.hh3aln\n")
                    fw.write(f"if [ ! -s 'seq.hh3aln' ];then\n")
                    fw.write(f"    cp {MSAdir}/dMSA.jacaln seq.hh3aln\n")
                    fw.write(f"fi\n")
                    fw.write(f"if [ ! -s 'seq.hh3aln' ];then\n")
                    fw.write(f"    cp {MSAdir}/dMSA.hhbaln seq.hh3aln\n")
                    fw.write(f"fi\n")
                    fw.write(f"sed = seq.hh3aln |sed 'N;s/\\n/\\t/'|sed 's/^/>/g'|sed 's/\\t/\\n/g'| {self.qMSApkg}/bin/qhmmbuild -n aln --amino -O seq.afq --informat afa seq.hmm -\n\n")
                    fw.write(f"{self.qMSApkg}/bin/qhmmsearch --cpu 1 -E 10 --incE 1e-3 -A {DBfasta}.match --tblout {DBfasta}.tbl -o {DBfasta}.out seq.hmm {mMSAJGI}/{DBfasta}\n")
                    fw.write(f"{self.qMSApkg}/bin/esl-sfetch -f {mMSAJGI}/{DBfasta} {DBfasta}.tbl|sed 's/*//g' > {DBfasta}.fseqs\n")
                    fw.write(f"{self.qMSApkg}/bin/cd-hit -i {DBfasta}.fseqs -o {JGIdir}/{DBfasta}.cdhit -c 1 -M 3000\n\n")
                    fw.write(f"sync\n")
                    fw.write(f"rm -rf {tmpdir}\n")
                    fw.write(f"echo ending time: `date`   >>{recorddir}/ware_{tag}\n")
  
                cmd = 'bash ' + bashfile
                process_list.append([cmd, 'JGI query'])
            
            if len(process_list) > 0:
                with ProcessPoolExecutor(max_workers=len(process_list)) as executor:
                    results = executor.map(run_command, process_list)

            print("2rd step/JGI combination is starting!\n")
            tag = "sJGI"

            tmpdir = f"{JGIdir}/tmp/{user_id}/{tag}"
            os.makedirs(tmpdir, exist_ok=True)
            cmd = ['python',
                os.path.join(self.tool_path, 'bin/JGImod.py'),
                outpath,
                tmpdir,
                qMSAhhblitsdb,
                qMSAjackhmmerdb,
                qMSAhhblits3db
                ]

            logging.info('Launching subprocess "%s"', ' '.join(cmd))
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with utils.timing('JGI combination'):
                stdout, stderr = process.communicate()
                retcode = process.wait()

            a3m_out = os.path.join(JGIdir, 'DeepJGI.a3m')
            if not os.path.exists(a3m_out):
                print(f"Cannot find the generated a3m file for {targetname} by DeepMS2: {a3m_out}")  

        deepmsa_noimg_tags = ['dMSA.hhb', 'dMSA.jac', 'dMSA.hms', 'dMSA', 'qMSA', 'aMSA', 'qMSA.hhb', 'qMSA.jac', 'qMSA.hh3', 'qMSA.hms']
        deepmsa_deepjgi_tags = ['DeepJGI.hms', 'DeepJGI']
        deepmsa_jgi_tags = ['q3JGI', 'q4JGI', 'q3JGI.hms', 'q4JGI.hms']

        finalMSAdir = os.path.join(outpath, 'finalMSAs')
        os.makedirs(finalMSAdir, exist_ok=True)
        
        for tag in deepmsa_noimg_tags:
            deepmsa_noimg_a3m = os.path.join(outpath, 'MSA', tag + 'a3m')
            if not os.path.exists(deepmsa_noimg_a3m):
                deepmsa_noimg_a3m = os.path.join(outpath, 'MSA', tag + '.a3m')
                if not os.path.exists(deepmsa_noimg_a3m):
                    continue
            contents = open(deepmsa_noimg_a3m).readlines()
            with open(f"{finalMSAdir}/{tag}.a3m", 'w') as fw:
                fw.write(f'>{targetname}\n')
                i = 1
                while i < len(contents):
                    if contents[i].find('>seq') == 0:
                        i += 1
                    else:
                        fw.write(contents[i])
                    i += 1
        
        for tag in deepmsa_deepjgi_tags:
            deepmsa_img_a3m = os.path.join(outpath, 'JGI', tag + 'a3m')
            if not os.path.exists(deepmsa_img_a3m):
                deepmsa_img_a3m = os.path.join(outpath, 'JGI', tag + '.a3m')
                if not os.path.exists(deepmsa_img_a3m):
                    continue
            contents = open(deepmsa_img_a3m).readlines()
            with open(f"{finalMSAdir}/{tag}.a3m", 'w') as fw:
                fw.write(f'>{targetname}\n')
                i = 1
                while i < len(contents):
                    if contents[i].find('>seq') == 0:
                        i += 1
                    else:
                        fw.write(contents[i])
                    i += 1
        
        for tag in deepmsa_jgi_tags:
            deepmsa_img_a3m = os.path.join(outpath, 'JGI', tag + 'a3m')
            if not os.path.exists(deepmsa_img_a3m):
                deepmsa_img_a3m = os.path.join(outpath, 'JGI', tag + '.a3m')
                if not os.path.exists(deepmsa_img_a3m):
                    continue
            contents = open(deepmsa_img_a3m).readlines()
            with open(f"{finalMSAdir}/{tag}.a3m", 'w') as fw:
                fw.write(f'>{targetname}\n')
                i = 1
                while i < len(contents):
                    if contents[i].find('>seq') == 0:
                        i += 1
                    else:
                        fw.write(contents[i])
                    i += 1
                
        return os.path.join(finalMSAdir, "dMSA.a3m")



