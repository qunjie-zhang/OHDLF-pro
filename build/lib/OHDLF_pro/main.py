import shutil
import time
import os
import sys
import subprocess
import uuid
import random
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

# Check Python dependencies
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
    from tqdm import tqdm
except ImportError as e:
    print("\033[91m[Error] Missing Python libraries.\033[0m")
    print(f"Details: {e}")
    print("Please install them using: pip install biopython tqdm")
    sys.exit(1)

# Check DISCO dependencies
try:
    import treeswift
except ImportError:
    pass


def check_external_dependencies(process_type):
    """
    Verify the presence of external software dependencies in the environment variables.
    """
    # 基础依赖
    required_tools = ["mafft", "makeblastdb", "blastp"]

    # Check specific tools according to the process type
    if process_type == 1:
        # Type 1 need RAxML
        required_tools.append("raxmlHPC-PTHREADS")

    elif process_type == 2:
        # Type 2 need IQ-TREE, ASTRAL
        required_tools.append("iqtree")
        required_tools.append("astral-pro3")

        # check treeswift dependencies
        try:
            import treeswift
        except ImportError:
            print("\033[91m[CRITICAL ERROR] Python library 'treeswift' is missing.\033[0m")
            print("It is required for DISCO (Process Type 2).")
            print("Please install it: \033[92mpip install treeswift\033[0m")
            sys.exit(1)


    missing_tools = []
    print("--- Checking external dependencies ---")
    for tool in required_tools:
        path = shutil.which(tool)
        if path:
            print(f"[{tool}]: Found at {path}")
        else:
            print(f"\033[91m[{tool}]: NOT FOUND\033[0m")
            missing_tools.append(tool)

    if missing_tools:
        print("\n\033[91m[CRITICAL ERROR] The following required tools are missing:\033[0m")
        print(f"  {', '.join(missing_tools)}")
        print("Please install them via conda or your package manager before running this script.")
        print("Example: \033[92mconda install -c bioconda mafft iqtree blast astral-pro3 raxml\033[0m")
        sys.exit(1)
    print("--------------------------------------\n")


def fasta_load(fasta_path):
    """Generator to read fasta file"""
    if not os.path.exists(fasta_path):
        return
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    seq_name = None
    seq = None
    for line in lines:
        if not line: break
        if line.startswith('>'):
            if seq:
                yield seq_name, seq
            seq = ''
            seq_name = line[1:].strip('\n')
        else:
            if seq is not None:
                seq += line.strip('\n')
    if seq:
        yield seq_name, seq


def run_command(command, quiet=True):

    try:
        stdout = subprocess.DEVNULL if quiet else None
        stderr = subprocess.DEVNULL if quiet else None
        subprocess.run(command, shell=True, check=True, stdout=stdout, stderr=stderr)
        return True
    except subprocess.CalledProcessError:
        return False


# ------------------------------------------------------------------
# 1. Filter CNV (Filter Orthologous)
# ------------------------------------------------------------------
def filter_orthologous(loss, duplication, output_dir):
    print(f"--- Step 1: Filtering Orthologous Groups (Max Loss: {loss}, Max Dup: {duplication}) ---")
    loss_num = float(loss)
    max_duplication_num = int(duplication)

    if not os.path.exists('Orthogroups') or not os.path.exists('Orthogroup_Sequences'):
        print("\033[91mError: Please place this script in the Orthofinder output directory (Results_*)\033[0m")
        print("Missing directories: Orthogroups/ or Orthogroup_Sequences/")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    tsv_file_path = 'Orthogroups/Orthogroups.GeneCount.tsv'
    with open(tsv_file_path, 'r', encoding='utf8') as f:
        tsv_lines = f.read().splitlines()

    species_name = [i for i in tsv_lines[0].split()[1:-1]]
    max_loss_allowed = round(loss_num * len(species_name))

    del tsv_lines[0]
    target_dict = {}

    for content in tsv_lines:
        parts = content.split()
        og_name = parts.pop(0)
        _ = parts.pop()
        counts = [int(x) for x in parts]

        if any(c > max_duplication_num for c in counts): continue
        if counts.count(0) > max_loss_allowed: continue

        if og_name not in target_dict:
            target_dict[og_name] = []

        for idx, count in enumerate(counts):
            if count == 0:
                target_dict[og_name].append(species_name[idx])

    if not target_dict:
        print("\033[93mWarning: No Orthogroups matched your criteria. Check -l and -d parameters.\033[0m")
        sys.exit(0)

    print(f"Found {len(target_dict)} orthogroups matching criteria.")
    for fasta_name, missing_species in tqdm(target_dict.items(), desc="Extracting Sequences"):
        source_file = f'Orthogroup_Sequences/{fasta_name}.fa'
        if not os.path.exists(source_file): continue

        with open(f'{output_dir}/{fasta_name}_repeat.fa', 'w', encoding='utf8') as file:
            for name, seq in fasta_load(source_file):
                file.write(f'>{name}\n{seq}\n')
            for sp in missing_species:
                file.write(f'>{sp}|EMPTY_DATA_{random.random()}\n')
                file.write(100 * '-' + '\n')


# ------------------------------------------------------------------
# 2. BLAST Filter (Parallel Worker)
# ------------------------------------------------------------------
def worker_blast_filter(file_info):
    filename, input_dir, output_dir, threshold, temp_root_dir = file_info
    input_path = os.path.join(input_dir, filename)
    uid = str(uuid.uuid4())

    temp1 = os.path.join(temp_root_dir, f"temp1_{uid}.fa")
    temp2 = os.path.join(temp_root_dir, f"temp2_{uid}.fa")
    temp_db = os.path.join(temp_root_dir, f"temp_db_{uid}")
    blast_out = os.path.join(temp_root_dir, f"blast_{uid}.out")

    def clean_temp_files():
        try:
            for f in [temp1, temp2, blast_out]:
                if os.path.exists(f): os.remove(f)
            for ext in [".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto"]:
                db_file = temp_db + ext
                if os.path.exists(db_file): os.remove(db_file)
        except OSError:
            pass

    try:
        species_sequences = defaultdict(list)
        for record in SeqIO.parse(input_path, "fasta"):
            sp_name = record.id.split('|')[0]
            species_sequences[sp_name].append(record)

        all_species_ok = True

        for sp, seqs in species_sequences.items():
            if len(seqs) > 1:
                SeqIO.write(seqs[0], temp1, "fasta")
                SeqIO.write(seqs, temp2, "fasta")

                cmd_db = f"makeblastdb -in {temp1} -dbtype prot -out {temp_db} -logfile {os.devnull}"
                if not run_command(cmd_db):
                    clean_temp_files()
                    raise RuntimeError("makeblastdb failed")

                cmd_blast = f"blastp -query {temp2} -db {temp_db} -out {blast_out} -outfmt 6"
                if not run_command(cmd_blast):
                    clean_temp_files()
                    raise RuntimeError("blastp failed")

                if not os.path.exists(blast_out) or os.path.getsize(blast_out) == 0:
                    all_species_ok = False
                    clean_temp_files()
                    break

                max_sims = defaultdict(float)
                with open(blast_out, 'r') as b_res:
                    for line in b_res:
                        cols = line.strip().split('\t')
                        if len(cols) < 3: continue
                        q, s, ident = cols[0], cols[1], float(cols[2])
                        max_sims[(q, s)] = max(max_sims[(q, s)], ident)

                current_sp_ok = True
                if not max_sims:
                    current_sp_ok = False
                else:
                    for val in max_sims.values():
                        if val < threshold:
                            current_sp_ok = False
                            break

                clean_temp_files()

                if not current_sp_ok:
                    all_species_ok = False
                    break

        clean_temp_files()

        if all_species_ok:
            shutil.copy(input_path, os.path.join(output_dir, filename))
            return 1
        return 0

    except Exception:
        clean_temp_files()
        return 0


def blast_identity_filter(sim, input_dir, output_dir, threads):
    print(f"--- Step 2: BLAST Identity Filtering (Threshold: {sim}%, Threads: {threads}) ---")
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    temp_workspace = "OHDLF_temp_workspace"
    if not os.path.exists(temp_workspace):
        os.makedirs(temp_workspace)
    else:
        shutil.rmtree(temp_workspace)
        os.makedirs(temp_workspace)

    print(f"Temporary files will be stored in: {temp_workspace}/")

    files = [f for f in os.listdir(input_dir) if f.endswith(('.fa', '.fasta'))]
    if not files:
        print("\033[93mWarning: No files found in filter step. Skipping BLAST.\033[0m")
        os.rmdir(temp_workspace)
        return

    tasks = [(f, input_dir, output_dir, sim, temp_workspace) for f in files]

    kept_count = 0
    try:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(worker_blast_filter, t) for t in tasks]
            for f in tqdm(as_completed(futures), total=len(tasks), desc="BLAST Filtering"):
                kept_count += f.result()
    finally:
        print("Cleaning up temporary workspace...")
        try:
            shutil.rmtree(temp_workspace)
        except Exception as e:
            print(f"\033[93mWarning: Could not remove temp dir {temp_workspace}: {e}\033[0m")

    print(f"BLAST Filter finished. Kept {kept_count}/{len(files)} files.")


# ------------------------------------------------------------------
# 3. MAFFT Align (Parallel Worker)
# ------------------------------------------------------------------
def worker_mafft(file_info):
    input_file, output_file = file_info
    if os.path.getsize(input_file) == 0:
        return False

    cmd = f"mafft --quiet --auto {input_file} > {output_file}"
    success = run_command(cmd)

    if success and (not os.path.exists(output_file) or os.path.getsize(output_file) == 0):
        return False
    return success


def align_by_mafft(input_dir, output_dir, threads):
    print(f"--- Step 3: MAFFT Alignment (Threads: {threads}) ---")
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    files = [f for f in os.listdir(input_dir) if f.endswith('.fa') or f.endswith('.fasta')]
    if not files:
        print("\033[93mWarning: No files to align.\033[0m")
        return

    tasks = []
    for f in files:
        in_path = os.path.join(input_dir, f)
        out_path = os.path.join(output_dir, f"{f}_align.fa")
        tasks.append((in_path, out_path))

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(worker_mafft, t) for t in tasks]
        for _ in tqdm(as_completed(futures), total=len(tasks), desc="Running MAFFT"):
            pass


# ------------------------------------------------------------------
# 4. multiple_seq_fill (Type 1 only)
# ------------------------------------------------------------------
def worker_fill(file_info):
    input_path, output_path = file_info
    if not os.path.exists(input_path) or os.path.getsize(input_path) == 0:
        return False

    try:
        records = list(SeqIO.parse(input_path, "fasta"))
        if not records: return False

        id_seq_dict = {}
        for record in records:
            seq_id = record.id.split("|")[0]
            seq_str = str(record.seq)
            if seq_id in id_seq_dict:
                existing_seq = id_seq_dict[seq_id]
                new_seq_list = list(existing_seq)
                length = min(len(seq_str), len(existing_seq))
                for i in range(length):
                    if seq_str[i] != existing_seq[i]:
                        if seq_str[i] == '-' or existing_seq[i] == '-':
                            new_seq_list[i] = '-'
                        else:
                            new_seq_list[i] = 'X'
                id_seq_dict[seq_id] = "".join(new_seq_list)
            else:
                id_seq_dict[seq_id] = seq_str

        with open(output_path, 'w', encoding='utf8') as f:
            for s_id, s_seq in id_seq_dict.items():
                f.write(f">{s_id}\n{s_seq}\n")
        return True
    except Exception:
        return False


def align_and_fill(input_dir, output_dir, threads):
    print(f"--- Step 4: Align and Fill (Threads: {threads}) ---")
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith('.fa')]

    tasks = [(os.path.join(input_dir, f), os.path.join(output_dir, f"{f}_fill.fa")) for f in files]
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(worker_fill, t) for t in tasks]
        for _ in tqdm(as_completed(futures), total=len(tasks), desc="Processing Seqs"):
            pass


# ------------------------------------------------------------------
# 5. Type1 (Merge -> RAxML)
# ------------------------------------------------------------------
def seq_merge(input_dir, output_file):
    print(f"--- Step 5 (Type 1): Merging Sequences into {output_file} ---")
    seq_dict = defaultdict(str)
    files = [f for f in os.listdir(input_dir) if f.endswith(".fa")]

    if not files:
        print("\033[93mWarning: No files to merge.\033[0m")
        return

    for file_name in tqdm(files, desc="Merging"):
        file_path = os.path.join(input_dir, file_name)
        if os.path.getsize(file_path) == 0: continue
        for name, seq in fasta_load(file_path):
            seq_dict[name] += seq

    if not seq_dict:
        print("\033[91mError: Merged sequence is empty.\033[0m")
        return

    with open(output_file, "w") as output:
        for s_id, s_seq in seq_dict.items():
            output.write(f">{s_id}\n{s_seq}\n")


def fas_to_phy_single(input_file, output_file=None):
    if not os.path.exists(input_file): return False
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}.phy"

    names, seqs = [], []
    max_name_len = 0
    try:
        with open(input_file, 'r') as f:
            content = f.read().replace(' ', '')
        lines = content.splitlines()
        current_seq, current_name = '', ''
        for line in lines:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_name:
                    names.append(current_name)
                    seqs.append(current_seq)
                current_name = line.split("|")[0][1:]
                current_seq = ''
                if len(current_name) > max_name_len: max_name_len = len(current_name)
            else:
                current_seq += line
        if current_name:
            names.append(current_name)
            seqs.append(current_seq)

        if not names: return False

        with open(output_file, 'w') as out:
            out.write(f"{len(names)} {len(seqs[0])}\n")
            for n, s in zip(names, seqs):
                space = ' ' * (max_name_len + 10 - len(n))
                out.write(f"{n}{space}{s}\n")
        return True
    except Exception:
        return False


def run_raxml_step(input_phy, threads):
    """
    Step 6 (Type 1): Running RAxML
    """
    print(f"--- Step 6 (Type 1): Running RAxML (Threads: {threads}) ---")

    if not os.path.exists(input_phy):
        print(f"\033[91m[Error] Input file '{input_phy}' is missing. Cannot run RAxML.\033[0m")
        return

    # cmd：raxmlHPC-PTHREADS -x 12345 -p 12345 -# 1000 -T {threads} -m PROTGAMMAIJTTF -s final_OrthologsAlign_GDL.phy -f a -n OHDLF_tree
    output_suffix = "OHDLF_tree"

    # remove old files of RAxML (if they exist)
    possible_outputs = [f"RAxML_info.{output_suffix}", f"RAxML_bestTree.{output_suffix}",
                        f"RAxML_bootstrap.{output_suffix}", f"RAxML_bipartitions.{output_suffix}",
                        f"RAxML_bipartitionsBranchLabels.{output_suffix}"]

    for f in possible_outputs:
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError:
                pass

    cmd = f"raxmlHPC-PTHREADS -x 12345 -p 12345 -# 1000 -T {threads} -m PROTGAMMAIJTTF -s {input_phy} -f a -n {output_suffix}"

    print(f"Executing: {cmd}")
    #  quiet=False help user see the Bootstrap process of RAxML
    if run_command(cmd, quiet=False):
        print(f"RAxML finished. Output files with suffix '{output_suffix}' generated.")
    else:
        print("\033[91m[Error] RAxML failed.\033[0m")


# ------------------------------------------------------------------
# 6. Type 2 : IQ-TREE  -> DISCO  -> ASTRAL
# ------------------------------------------------------------------

def worker_iqtree_direct(args):
    """
    args: (input_fasta_path, output_prefix)
    """
    input_path, output_prefix = args
    cmd = f"iqtree -s {input_path} -pre {output_prefix} -B 1000 --bnni -T 1 --quiet -redo"
    return run_command(cmd)


def run_iqtree_step(input_dir, output_dir, threads):
    print(f"--- Step 4 (Type 2): Running IQ-TREE on Alignments (Threads: {threads}) ---")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    files = [f for f in os.listdir(input_dir) if f.endswith('.fa') or f.endswith('.fasta')]

    if not files:
        print("\033[93mWarning: No alignment files found for IQ-TREE.\033[0m")
        return

    tasks = []
    for f in files:
        in_path = os.path.join(input_dir, f)
        base_name = f.replace(".fa", "").replace(".fasta", "")
        out_prefix = os.path.join(output_dir, base_name)
        tasks.append((in_path, out_prefix))

    print(f"Submitting {len(tasks)} IQ-TREE jobs...")
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(worker_iqtree_direct, t) for t in tasks]
        for _ in tqdm(as_completed(futures), total=len(tasks), desc="IQ-TREE Building"):
            pass

    print(f"IQ-TREE finished. Trees saved in {output_dir}")



def worker_disco(args):
    """
    args: (input_tree_file, output_disco_file)
    调用包内部的 utils/disco.py 脚本
    """
    input_tree, output_tree = args

    # 1. 动态获取 utils/disco.py 的绝对路径
    # __file__ 是当前 main.py 的路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # 拼接路径
    disco_script = os.path.join(current_dir, "utils", "disco.py")

    if not os.path.exists(disco_script):
        print(f"\033[91m[Error] Internal script missing: {disco_script}\033[0m")
        return False

    # 2. 使用 sys.executable 确保环境一致
    cmd = f"{sys.executable} {disco_script} -i {input_tree} -d \"|\" -o {output_tree}"

    if run_command(cmd):
        return True
    return False


def run_disco_step(input_dir, output_dir, threads):
    print(f"--- Step 5 (Type 2): Running DISCO decomposition (Threads: {threads}) ---")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    files = [f for f in os.listdir(input_dir) if f.endswith('.treefile')]

    if not files:
        print("\033[91mError: No .treefile found in IQ-TREE output directory.\033[0m")
        return

    tasks = []
    for f in files:
        in_path = os.path.join(input_dir, f)
        out_name = f"DISCO_{f}"
        out_path = os.path.join(output_dir, out_name)
        tasks.append((in_path, out_path))

    print(f"Submitting {len(tasks)} DISCO jobs...")
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(worker_disco, t) for t in tasks]
        for _ in tqdm(as_completed(futures), total=len(tasks), desc="Running DISCO"):
            pass

    print(f"DISCO finished. Decomposed trees saved in {output_dir}")

    all_trees_path = "all_disco.trees"
    count = 0
    with open(all_trees_path, 'w') as outfile:
        for f in os.listdir(output_dir):
            if f.startswith("DISCO_") and (f.endswith(".treefile") or ".treefile" in f):
                path = os.path.join(output_dir, f)
                if os.path.getsize(path) > 0:
                    with open(path, 'r') as infile:
                        outfile.write(infile.read())
                        count += 1
    if count > 0:
        print(f"Summary: Concatenated {count} DISCO trees into '{all_trees_path}'")


def run_astral_step(input_file, output_file, threads):
    """
    Step 6 (Type 2): Running ASTRAL-Pro3
    """
    print(f"--- Step 6 (Type 2): Running ASTRAL-Pro3 (Threads: {threads}) ---")

    if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
        print(f"\033[91m[Error] Input file '{input_file}' is missing or empty. Cannot run ASTRAL.\033[0m")
        return

    cmd = f"astral-pro3 -t {threads} -R -i {input_file} -o {output_file}"

    print(f"Executing: {cmd}")
    if run_command(cmd, quiet=False):
        print(f"ASTRAL-Pro3 finished. Output saved to: {output_file}")
    else:
        print("\033[91m[Error] ASTRAL-Pro3 failed.\033[0m")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    arg = argparse.ArgumentParser(description='OHDLF-pro: a pipeline designed to filter and address orthologous gene heterogeneity, duplication, and loss')
    arg.add_argument("-l", "--loss", default=0, required=True, help="Allowable loss value")
    arg.add_argument("-d", "--duplication", default=3, required=True, help="Allowable max copy value")
    arg.add_argument("-s", "--similarity", default=97, required=False, help="Similarity value")
    arg.add_argument("-p", "--process_type", type=int, choices=[1, 2], required=True,
                     help="1: Concatenation (RAxML), 2: Coalescence (ASTRAL)")
    arg.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")

    args = arg.parse_args()


    check_external_dependencies(args.process_type)

    start_time = time.time()

    suffix = "_parallel" if args.process_type == 2 else ""
    dir_all = f"all_GDL_Orthologue_Sequences{suffix}"
    dir_blast = f"GDL_Orthologue_Sequences{suffix}"
    dir_mafft = f"GDL_Orthologue_Sequences_mafft{suffix}"
    # Type 1
    dir_fill = f"GDL_Orthologue_Sequences_mafft_fill{suffix}"
    # Type 2
    dir_iqtree = "GDL_Orthologue_Sequences_iqtree"
    dir_disco = "GDL_Orthologue_Sequences_DISCO"

    # 1. Filter
    filter_orthologous(args.loss, args.duplication, dir_all)

    # 2. BLAST
    blast_identity_filter(float(args.similarity), dir_all, dir_blast, args.threads)

    # 3. MAFFT
    align_by_mafft(dir_blast, dir_mafft, args.threads)

    if args.process_type == 1:
        # Type 1 : Fill -> Merge -> Phy -> RAxML
        # 4. Fill
        align_and_fill(dir_mafft, dir_fill, args.threads)

        # 5. Result - Concatenation
        final_out = "final_OrthologsAlign_GDL.fasta"
        seq_merge(dir_fill, final_out)

        # 6. Convert to Phy and Run RAxML
        fas_to_phy_single(final_out)

        #  phy filename (default final_OrthologsAlign_GDL.phy)
        final_phy = final_out.replace(".fasta", ".phy")
        run_raxml_step(final_phy, args.threads)

        print(f"Analysis Finished (Type 1). Check RAxML outputs.")

    elif args.process_type == 2:
        # Type 2 : Skip Fill -> IQ-TREE (Parallel) -> DISCO -> ASTRAL

        # 4. IQ-TREE Direct (Input: MAFFT dir, Output: IQ-TREE dir)
        run_iqtree_step(dir_mafft, dir_iqtree, args.threads)

        # 5. DISCO (Input: IQ-TREE dir, Output: DISCO dir)
        run_disco_step(dir_iqtree, dir_disco, args.threads)

        # 6. ASTRAL (Input: all_disco.trees from prev step, Output: Root dir)
        run_astral_step("all_disco.trees", "OHDLF_DISCO_ASTRAL.nwk", args.threads)

        print(f"Analysis Finished (Type 2). Check 'OHDLF_DISCO_ASTRAL.nwk'")

    print(f"Total Runtime: {(time.time() - start_time) / 60:.2f} minutes")

def run():
    arg = argparse.ArgumentParser(
        description='OHDLF-pro: a pipeline designed to filter and address orthologous gene heterogeneity, duplication, and loss')
    arg.add_argument("-l", "--loss", default=0, required=True, help="Allowable loss value")
    arg.add_argument("-d", "--duplication", default=3, required=True, help="Allowable max copy value")
    arg.add_argument("-s", "--similarity", default=97, required=False, help="Similarity value")
    arg.add_argument("-p", "--process_type", type=int, choices=[1, 2], required=True,
                     help="1: Concatenation (RAxML), 2: Coalescence (ASTRAL)")
    arg.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")

    args = arg.parse_args()
    main(args)


if __name__ == "__main__":
    run()