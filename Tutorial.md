# MeH v2 tutorial

MeH profiles DNA methylation heterogeneity from bisulfite-sequencing BAM
files. Version 2 requires Python 3.8 or newer and supports current pandas
releases.

## Installation

Clone the repository and switch to the v2 branch:

```bash
git clone https://github.com/britishcoffee/MeH.git
cd MeH
git switch v2
```

Install MeH and its Python dependencies:

```bash
python3 -m pip install -e .
```

After installation, either `mehscr` or `python3 MeHscr.py` can be used. Display
all available options with:

```bash
mehscr --help
```

## Input files

Place the following files in one input directory:

```text
input/
├── sample1.bam
├── sample1.bam.bai
├── sample2.bam
├── sample2.bam.bai
├── genome.fa
└── genome.fa.fai
```

Requirements:

- BAM files must be coordinate-sorted and indexed.
- Each BAM must use the same reference genome.
- The reference FASTA must be indexed.
- PCR duplicates should preferably be removed before running MeH.
- Sample names are taken from the complete BAM filename stem. For example,
  `treated_rep1.bam` becomes sample `treated_rep1`.

MeH v2 reads BAM files directly. It does not split, modify, or delete input
files.

## Run the included example

The repository includes four small BAM files and one reference in `MeHdata`.
From the repository root, run:

```bash
python3 MeHscr.py -w 4 -c 4 --CG
```

Because `MeHdata` is the default input and output directory, the equivalent
explicit command is:

```bash
python3 MeHscr.py \
  --CG \
  --windowsize 4 \
  --cores 4 \
  --data-dir ./MeHdata \
  --output-dir ./MeHdata
```

## Analyse your own data

Keeping results separate from input data is recommended:

```bash
python3 MeHscr.py \
  --CG \
  --data-dir /path/to/input \
  --output-dir /path/to/results \
  -w 4 \
  -c 4
```

If the input directory contains multiple FASTA files, select one explicitly:

```bash
python3 MeHscr.py \
  --CG \
  --data-dir /path/to/input \
  --reference /path/to/genome.fa \
  --output-dir /path/to/results
```

At least one sequence context must be selected. Multiple contexts can be
analysed in one command:

```bash
python3 MeHscr.py --CG --CHG --CHH \
  --data-dir /path/to/input \
  --output-dir /path/to/results
```

## Important options

| Option | Description |
|---|---|
| `-w`, `--windowsize` | Cytosines per window; allowed values are 2–6 (default: 4). |
| `-c`, `--cores` | Maximum number of BAM/context jobs run in parallel (default: 4). |
| `-m`, `--MeH` | Score: 1 abundance, 2 pairwise, 3 phylogeny, 4 entropy, or 5 epipolymorphism. |
| `-d`, `--dist` | Pattern distance: 1 Hamming or 2 weighted-degree kernel. |
| `--CG`, `--CHG`, `--CHH` | Sequence contexts to analyse. |
| `--opt` | Write methylation-pattern counts for each eligible window. |
| `--mlv` | Calculate and write methylation levels. |
| `--imp` | Enable BSImp imputation for reads with one missing value. |
| `--seed` | Seed used by imputation, allowing reproducible runs (default: 0). |
| `--progress-every` | Print progress after this percentage interval (default: 5). |
| `--quiet` | Suppress per-BAM progress messages. |

The phylogeny score (`-m 3`) supports window sizes 3 and 4 only. Invalid
combinations are rejected before analysis starts.

### Pattern counts and methylation levels

```bash
python3 MeHscr.py --CG --opt --mlv \
  --data-dir /path/to/input \
  --output-dir /path/to/results
```

### Reproducible imputation

```bash
python3 MeHscr.py --CG --imp --seed 42 \
  --data-dir /path/to/input \
  --output-dir /path/to/results
```

Using the same input, options, and seed produces the same imputation choices.

## Output files

For two samples named `sample1` and `sample2`, a CG analysis can produce:

```text
results/
├── MeHscreening.log
├── CG_sample1_0.csv
├── CG_sample2_0.csv
└── CG_Results.csv
```

The per-sample files contain:

| Column | Description |
|---|---|
| `chrom` | Reference sequence name. |
| `pos` | One-based position of the first cytosine in the window. |
| `MeH` | Methylation heterogeneity score. |
| `dis` | Genomic distance from the first to last cytosine. |
| `strand` | `f` for forward or `r` for reverse. |

`CG_Results.csv` contains mean MeH values in 400 bp bins, with one column per
sample. CHG and CHH analyses use the same naming scheme.

Additional options create these files:

- `--opt`: `CG_opt_sample1.csv`, containing `p01`, `p02`, etc. pattern counts.
- `--mlv`: `CG_ML_sample1_0.csv` and the binned `CG_ML_Results.csv` summary.

## Convert results for IGV

The included R script converts a summary file to bedGraph format:

```bash
Rscript tobed.R -m ./results/CG_Results.csv
```

Use `Rscript tobed.R -h` for all options.

## Find differential heterogeneous regions

`finddhr.R` requires a MeH summary, a sample-condition definition, and a gene
list:

```bash
Rscript finddhr.R \
  -m ./results/CG_Results.csv \
  -g ./MeHdata/genelist.txt \
  -o ./results/CG \
  -s W,W,D,D \
  -p 1000
```

Typical outputs are:

- `CG_DHR_Result.csv`: differential heterogeneous regions associated with
  genes or promoters.
- `CG_MeH_Result.csv`: MeH differences used in the DHR analysis.

Use `Rscript finddhr.R -h` for filtering and statistical options.

## Migrating from the original version

The main command-line flags remain compatible, but v2 fixes duplicate read
counting caused by many-to-many DataFrame merges in the original script.
Consequently, the number of eligible windows and some MeH values can differ
from results produced by the old implementation. Re-run all samples in a
comparison with the same MeH version rather than mixing old and v2 outputs.
