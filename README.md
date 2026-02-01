# OHDLF-pro
## Introduction
**OHDLF-pro** is the advanced evolution of the OHDLF pipeline, designed for high-performance phylogenomic analysis. Building upon OrthoFinder outputs, it automates the screening, alignment, and phylogenetic tree construction process with enhanced speed and accuracy.

Compared to the standard OHDLF, the **Pro** version introduces critical upgrades:

1. **⚡ Multi-threading Support**: New parallel computing architecture (via the `-t` parameter) significantly accelerates BLAST filtering, alignment, and tree building steps.
    
2. **✨ Visual Progress Tracking**: Integrated **dynamic animations** (progress bars) provide real-time feedback on current tasks and estimated completion time.
    
3. **🛡️ Robust Coalescence (Type 2)**: Now integrates **DISCO** (Decomposition of Species Trees) into the Coalescence workflow. This allows for the robust handling of **multi-copy gene families**, decomposing complex gene trees into single-copy orthologs to yield more statistically reliable species trees compared to simple filtering.

## dependencies
- Biopython
- tqdm
- treeswift
External Bioinformatics Tools (conda Recommended)
```
conda install -c bioconda iqtree aster blast raxml mafft
```
users can download the OHDLF.yaml file to directly configure the environment.
## Install
```
pip install OHDLF-pro
```

## Quick Start
```
Concatenation:
OHDLF-pro -l 0.05 -d 6 -s 97 -p 1 -t 8
Coalescence:
OHDLF-pro -l 0 -d 6 -s 97 -p 2 -t 8
```
## Usage
```
OHDLF-pro.py -l [LOSS] -d [DUP] -s [SIM] -p [TYPE] -t [THREADS]


Commands:
  -l / --loss Allowable the max missing rate of gene. This option is required.
  -d / --duplication Allowable the max duplication number of gene. This option is required.
  -s / --similarity Allowable the similarity threshold of gene. If you do not set this parameter, the program will use '97' by default
  -p / --process_type process_type: 1 for Concatenation, 2 for Coalescence
  -t / --threads number of allowed threads
```

## Input
**Input** :You need to 'cd' to the Orthofinder output directory named 'Results_XXX'. The software depends on two directories: 'Orthogroup_Sequences' and 'Orthogroups'.

## Output

### Type 1 (Concatenation)

- `final_OrthologsAlign_GDL.phy`: The concatenated alignment (phy file).
    
- `RAxML_bestTree.OHDLF_tree`: The final Maximum Likelihood tree.

### Type 2 (Coalescence with DISCO)

- `GDL_Orthologue_Sequences_iqtree`: Individual gene trees inferred by IQ-TREE.
    
- `GDL_Orthologue_Sequences_DISCO`: **Decomposed gene trees** processed by DISCO (Multi-copy -> Single-copy).
    
- `all_disco.trees`: The combined input file for ASTRAL.
    
- **`OHDLF_DISCO_ASTRAL.nwk`**: The final, robust species tree.