import sys
import os
import re
import glob
import subprocess, argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bash_file', type=str, required=True)
    parser.add_argument('--refdir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    contents = open(args.bash_file).readlines()

    cmd_info = contents[-1]
    if cmd_info.find('alphafold') < 0:
        raise Exception(f"Cannot find the cmd info: {contents[-1]}")

    # Regular expression to find the model_preset value
    match = re.search(r'--model_preset=([\w-]+)', cmd_info)
    if match:
        model_preset_value = match.group(1)
        print("model_preset value:", model_preset_value)
    else:
        raise Exception(f"Cannot find the cmd info: {contents[-1]}")
    
    if model_preset_value == "monomer":
        match = re.search(r'--num_monomer_predictions_per_model ([\w-]+)', cmd_info)
        if match:
            num_monomer_predictions_per_model = match.group(1)
            print("num_monomer_predictions_per_model value:", num_monomer_predictions_per_model)
        else:
            raise Exception(f"Cannot find the cmd info: {contents[-1]}")

        task_name = os.path.basename(args.bash_file).replace('.sh', '')
        os.makedirs(args.outdir, exist_ok=True)
        monomer_ckpts = ['model_1', 'model_2', 'model_3', 'model_4', 'model_5']
        for monomer_ckpt in monomer_ckpts:
            for i in range(1, int(num_monomer_predictions_per_model), 10):
                single_task_cmd_info = cmd_info + f"--model_ckpt={monomer_ckpt} --nstruct_start={i}"
                single_task_cmd_info.replace('TOPN', 'NONE')
                single_task_file = os.path.join(args.outdir, f"{task_name}_{monomer_ckpt}_{i}.sh")

                new_contents = copy.deepcopy(contents)
                new_contents[-1] = single_task_cmd_info
                with open(single_task_file, 'w') as fw:
                    fw.writelines(new_contents)

                bfinished = True
                for j in range(i, 10):
                    model_name = f'unrelaxed_{monomer_ckpt}_pred_{j}.pdb'
                    if not os.path.exists(args.refdir + '/' + model_name):
                        bfinished = False
                        break

                if not bfinished:
                    os.system('sbatch ' + single_task_file)
    else:
        match = re.search(r'--num_multimer_predictions_per_model=([\w-]+)', cmd_info)
        if match:
            num_multimer_predictions_per_model = match.group(1)
            print("num_multimer_predictions_per_model value:", num_multimer_predictions_per_model)
        else:
            raise Exception(f"Cannot find the cmd info: {contents[-1]}")

        task_name = os.path.basename(args.bash_file).replace('.sh', '')
        os.makedirs(args.outdir, exist_ok=True)
        
        #if model_preset_value == "multimer":
        multimer_ckpts = ['model_1_multimer_v3', 'model_2_multimer_v3', 
                          'model_3_multimer_v3', 'model_4_multimer_v3',
                          'model_5_multimer_v3']

        for multimer_ckpt in multimer_ckpts:
            for i in range(0, int(num_multimer_predictions_per_model), 2):
                single_task_cmd_info = cmd_info + f"--model_ckpt={multimer_ckpt} --nstruct_start={i}"
                single_task_cmd_info.replace('TOPN', 'NONE')
                single_task_file = os.path.join(args.outdir, f"{task_name}_{multimer_ckpt}_{i}.sh")

                new_contents = copy.deepcopy(contents)
                new_contents[-1] = single_task_cmd_info
                with open(single_task_file, 'w') as fw:
                    fw.writelines(new_contents)

                bfinished = True
                for j in range(i, 10):
                    model_name = f'unrelaxed_{multimer_ckpt}_pred_{j}.pdb'
                    if not os.path.exists(args.refdir + '/' + model_name):
                        bfinished = False
                        break

                # if not bfinished:
                #     os.system('sbatch ' + single_task_file)


