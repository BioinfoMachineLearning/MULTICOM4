import os
import argparse
import concurrent.futures
from functools import partial

def makedir_if_not_exists(directory):
    """Creates a directory if it does not exist."""
    os.makedirs(directory, exist_ok=True)

def process_single_structure(native_pdb, native_dir, pred_dir, out_dir, pre2zhang_script):
    """Processes a single predicted structure using pre2zhang.pl."""
    if not native_pdb.endswith('.pdb'):
        return
    
    target_name = native_pdb.replace('.pdb', '')
    target_fullname = target_name.split('-')[0]  # Extract target name (e.g., T1201 from T1201-D1.pdb)

    pred_target_dir = os.path.join(pred_dir, target_fullname)
    if not os.path.isdir(pred_target_dir):
        print(f"Predicted structures for {target_fullname} not found.")
        return

    target_out_dir = os.path.join(out_dir, target_name)
    makedir_if_not_exists(target_out_dir)

    target_mol2_out_dir = os.path.join(out_dir, 'MOL2', target_name)
    makedir_if_not_exists(target_mol2_out_dir)

    for pred_pdb in os.listdir(pred_target_dir):
        input_pdb_path = os.path.join(pred_target_dir, pred_pdb)
        output_pdb_path = os.path.join(target_out_dir, f"{pred_pdb}_filtered.pdb")
        combined_pdb_path = os.path.join(target_mol2_out_dir, f"{pred_pdb}_combined.pdb")
        
        if os.path.exists(output_pdb_path) and os.path.exists(combined_pdb_path):
            continue

        os.system(f"perl {pre2zhang_script} {input_pdb_path} {os.path.join(native_dir, native_pdb)} {output_pdb_path}")
        print(f"Processed {pred_pdb} -> {output_pdb_path}")

        # Create a combined output file with both native and filtered predicted structure
        with open(combined_pdb_path, 'w') as combined_file:
            combined_file.write("MOLECULE\n")
            with open(os.path.join(native_dir, native_pdb), 'r') as native_file:
                combined_file.writelines(native_file.readlines())
            combined_file.write("END\n")
            combined_file.write("MOLECULE\n")
            with open(output_pdb_path, 'r') as pred_file:
                combined_file.writelines(pred_file.readlines())
            combined_file.write("END")
            print(f"Created combined file {combined_pdb_path}")

def process_predicted_structures(native_dir, pred_dir, out_dir, pre2zhang_script, num_workers=4):
    """Runs pre2zhang.pl on predicted structures in parallel."""
    makedir_if_not_exists(out_dir)

    native_pdbs = [pdb for pdb in os.listdir(native_dir) if pdb.endswith('.pdb')]

    # Use multiprocessing to speed up processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        executor.map(partial(process_single_structure, native_dir=native_dir, pred_dir=pred_dir, 
                             out_dir=out_dir, pre2zhang_script=pre2zhang_script), native_pdbs)

def main():
    parser = argparse.ArgumentParser(description="Process predicted structures using pre2zhang.pl in parallel.")
    parser.add_argument("-n", "--nativedir", help="Directory containing native PDB structures.", type=str, required=True)
    parser.add_argument("-p", "--preddir", help="Directory containing predicted structures.", type=str, required=True)
    parser.add_argument("-o", "--outdir", help="Output directory for processed PDBs.", type=str, required=True)
    parser.add_argument("-s", "--script", help="Path to pre2zhang.pl script.", type=str, required=True)
    parser.add_argument("-w", "--workers", help="Number of parallel workers.", type=int, default=150)

    args = parser.parse_args()
    process_predicted_structures(args.nativedir, args.preddir, args.outdir, args.script, args.workers)

if __name__ == "__main__":
    main()
