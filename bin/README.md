# busco_multigene_tree.py

This Python script automates extraction of BUSCO single-copy orthologs across multiple genomes, builds alignments at a user-specified completeness threshold, and infers a concatenated phylogeny with IQ-TREE using partition models.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Input directory structure](#input-directory-structure)
- [Output](#output)
- [Authors](#authors)


## Features

- **collect**: gather per-gene FASTA files from BUSCO outputs
- **select**: choose genes present in at least X% of samples
- **align**: align & trim genes with MAFFT & trimAl
- **infer**: concatenate alignments (AMAS) and infer a tree (IQ-TREE)
- **all**: run the full pipeline (collect → select → align → infer)


## Prerequisites

- Linux or macOS
- [miniforge](https://github.com/conda-forge/miniforge)
- Python 3.10+


## Installation

1. Clone the repository
   ```bash
   git clone <your-repo-url>
   cd <your-repo-name>
   ```

2. Create and activate a conda environment
	```bash
	conda create -n busco_phylogenomics 'python>=3.10' \
		biopython tqdm mafft trimal iqtree amas \
		-c conda-forge -c bioconda
	conda activate busco_phylogenomics
	```

3. Make the script executable
	```bash
	chmod +x busco_multigene_tree.py
	```


## Usage
```bash
busco_multigene_tree.py <subcommand> [options]
busco_multigene_tree.py -h
busco_multigene_tree.py collect --help
```

## Fraction
- `--fraction 0.8` builds a phylogenomic tree with genes that 80% of the samples possess
- Multiple fractions can be selected with comma (`--fraction 0.8,0.9,1.0`)

## Input directory structure

The script expects the input directory below
```
input/
├── sample_A
│   ├── busco_output
│   │   ├── logs/ …
│   │   ├── prodigal_output/ …
│   │   ├── run_<busco_lineage>			# e.g. bacillota_odb12
│   │   │   ├── busco_sequences
│   │   │   │   ├── fragmented_busco_sequences
│   │   │   │   ├── multi_copy_busco_sequences
│   │   │   │   └── single_copy_busco_sequences
│   │   │   │       ├── geneA.faa
│   │   │   │       ├── geneA.fna
│   │   │   │       ├── geneB.faa
│   │   │   │       ├── geneB.fna
│   │   │   │       └── …
│   │   │   ├── full_table.tsv
│   │   │   ├── hmmer_output/ …
│   │   │   ├── missing_busco_list.tsv
│   │   │   ├── short_summary.json
│   │   │   └── short_summary.txt
│   │   ├── short_summary.specific.<busco_lineage>.busco.json
│   │   └── short_summary.specific.<busco_lineage>.busco.txt
│   └── other_dirs/ …
├── sample_B/ …
├── sample_C/ …
└── …
```


## Output
```
output/
├── seqs/
│   ├── raw/           # raw per-gene FASTA files
│   ├── aligned/       # MAFFT alignments
│   └── trimmed/       # trimAl-trimmed alignments
├── frac90pct_results/ # results for 90% completeness
│   ├── frac90pct_genes.txt
│   ├── concat.faa
│   ├── partitions.nex
│   └── frac90pct.iqtree.*  # IQ-TREE outputs
└── frac80pct_results/ …    # if used 0.8 threshold
```

if `seqs/raw`, `seqs/aligned`, ir `seqs/trimmed` already contains files, the software will abort.

## Authors
- Akito Shima (ASUQ) — akito.shima@oist.jp
