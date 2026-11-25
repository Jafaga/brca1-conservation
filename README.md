# Sequence conservation across species

## What this project investigates

**How conserved is the same protein-coding gene across species?** I am investigating whether human DNA has greater similarity to orthologous sequences from more closely related species. An ortholog is a corresponding gene in different species inherited from a common ancestor.

The program implements **Smith–Waterman local alignment from scratch**. It compares coding DNA sequences (CDS), traces back a highest-scoring matching region, and records alignment score, percent identity, alignment length, mismatches, gaps, coverage, runtime, and memory use.

The terminal demo has three examples, **BRCA1, TP53, and BRCA2**. Each uses human, chimpanzee, mouse, and cat sequences from NCBI. The exact transcripts and lengths for each gene are in its manifest. BRCA1 remains the original checkpoint dataset. TP53 and BRCA2 show how the results vary by gene.

### Hypothesis and what the results tell us

Under fixed alignment scoring, a closely related species is expected to show higher similarity to human in corresponding gene sequences. The program provides preliminary evidence about conservation in the selected regions and measures the computational cost of obtaining those results.

Based on the proposal feedback, the analysis uses **specific orthologous genes with documented sequences and lengths**, a manually checkable toy example, and real public data. It measures scores, identity, and coverage. A count of local alignments would be misleading because Smith–Waterman can find local matches even between unrelated inputs.

Interpret a high identity together with coverage and score. A short local match can have high identity while covering little of the input. These similarity measures are descriptive evidence of conservation. Estimating evolutionary distance would require a defined evolutionary model and further analysis.

### Connection to distance-based phylogeny

The project uses the first two steps from lecture. It starts with orthologous coding sequences and aligns them to identify matching positions. This pilot compares each species with human using local alignment. A phylogeny extension would align comparable homologous positions across all four species, calculate distances for all six species pairs under a consistent rule, and then use the distance matrix to build a tree. The current alignment scores, identity, and coverage describe the selected windows.

### Current scope

Each selected gene has three non-overlapping pieces per species. The program runs nine real comparisons and nine composition-preserving shuffled controls. The piece length is the smaller of 600 bases and one third of the shortest selected CDS. The terminal reports that length for each run. This keeps a live demonstration manageable while ensuring all species contribute complete pieces.

Same-offset pieces may shift out of correspondence after insertions or deletions. Windows within one gene are not independent biological replicates. One shuffle per piece gives a preliminary control. Full-length comparisons, broader ortholog and gene sampling, and HPC benchmarking are the next steps toward the planned large-scale analysis.

## Set up on your Mac

Open **Terminal** or **Terminal → New Terminal** in VS Code. Run commands from the project folder. If you already have a clone, open that folder and run `git pull --ff-only origin main`. For a new copy, run:

```bash
git clone https://github.com/Jafaga/brca1-conservation.git
cd brca1-conservation
```

Use Python 3.10 or newer. You can reuse the virtual environment you already set up, or create one with `python3 -m venv .venv`. Activate it before running the demo.

```bash
source .venv/bin/activate
python --version
```

The terminal alignment and summaries use the Python standard library. For optional PNG/PDF plots or the validation tests, install dependencies with `python -m pip install -r requirements.txt`.

## Download official NCBI files

Install the [NCBI Datasets command-line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) if needed. Each command below downloads one gene's ortholog package with coding DNA sequences.

| Gene | Official human gene record | Downloaded file |
|---|---|---|
| BRCA1 | [NCBI Gene 672](https://www.ncbi.nlm.nih.gov/gene/672) | `brca1_orthologs.zip` |
| TP53 | [NCBI Gene 7157](https://www.ncbi.nlm.nih.gov/gene/7157) | `tp53_orthologs.zip` |
| BRCA2 | [NCBI Gene 675](https://www.ncbi.nlm.nih.gov/gene/675) | `brca2_orthologs.zip` |

```bash
datasets download gene gene-id 672 --ortholog all --include cds --filename brca1_orthologs.zip
datasets download gene gene-id 7157 --ortholog all --include cds --filename tp53_orthologs.zip
datasets download gene gene-id 675 --ortholog all --include cds --filename brca2_orthologs.zip
```

Use your existing BRCA1 ZIP if you already downloaded it. Supply its actual path after `--package`. Keep the ZIP intact. The program reads `ncbi_dataset/data/cds.fna` for the coding DNA and `ncbi_dataset/data/data_report.jsonl` for the gene metadata. See the [NCBI gene package documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/gene-package/).

## Run the live demo in Terminal

Start with a hand-checkable example. The best local match is `ACGT`, giving a score of 8, 100% identity within the match, and 50% coverage of each input.

```bash
python demo.py --query TTACGTAA --target GGACGTCC
```

Then run one command per official ZIP. Each command imports the downloaded data and computes new alignments.

```bash
python demo.py --gene BRCA1 --package brca1_orthologs.zip --output runs/brca1
python demo.py --gene TP53 --package tp53_orthologs.zip --output runs/tp53
python demo.py --gene BRCA2 --package brca2_orthologs.zip --output runs/brca2
```

If you need to rehearse offline, use the bundled NCBI snapshots. For example, run `python demo.py --data-dir data/TP53 --output runs/tp53-rehearsal`. BRCA1 is in `data` and BRCA2 is in `data/BRCA2`.

### Reading the output

Each row averages the three real windows for one comparison to human and includes the mean shuffled score for context. Percent identity uses all alignment columns, including gaps. Human coverage shows how much of the human input window appears in the alignment. For comparisons across genes at a common window length, add `--window 300` to each command.

### Saved outputs

| File | What it contains |
|---|---|
| `alignments.csv` | One row per comparison with score, identity, coverage, and timing |
| `alignments.json` | Per-comparison results and the aligned DNA strings |
| `benchmark.json` | Input size, runtime, memory, scoring parameters, and window sizes |
| `input/manifest.json` | Package name, accessions, lengths, and checksums for an imported ZIP |
| `input/source_package.zip` | A copy of the official ZIP used for the run |

Add `--plots` to a run to save PNG and PDF figures. See [Methods and reproducibility](docs/METHODS.md) for the scoring recurrence, transcript selection, and scaling discussion. The original BRCA1 [conservation figure](results/conservation.png) and [benchmark measurements](results/benchmark.json) are included.
