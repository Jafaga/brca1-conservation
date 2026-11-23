# Methods and reproducibility

## Research question and professor feedback

Does human BRCA1 have greater sequence similarity to a closely related mammalian ortholog than to more distant mammalian orthologs under fixed alignment scoring?

The broader hypothesis predicts higher similarity for more closely related species. This checkpoint tests the computational workflow with one gene and three comparison species. Conservation is quantified using alignment score, identity and coverage. Transcript variants are not treated as independent species. The pilot includes both a manually checkable toy case and public sequence data with exact accessions and lengths. These choices address the feedback recovered from the prior project discussion.

## Selected real data

| Species | NCBI Gene | RefSeq transcript | CDS bases |
|---|---|---|---|
| Human, Homo sapiens | 672 | NM_007294.4 | 5592 |
| Chimpanzee, Pan troglodytes | 449497 | NM_001045493.1 | 5592 |
| Mouse, Mus musculus | 12189 | NM_009764.3 | 5439 |
| Cat, Felis catus | 101081937 | XM_019817933.3 | 5622 |

Cat is a predicted XM_ transcript. It was included as an additional mammalian comparison. Zebrafish was not found in the queried human BRCA1 ortholog set. The pilot therefore uses the four named mammalian genes above. Reassess transcript selection before the broader biological analysis.

The bundled manifest records the actual retrieval timestamp, source URL and SHA256 checksums. Sequences were downloaded afresh for this prototype. The earlier 562-gene download remains on your Mac. This new package contains four genes and is not a copy of that local ZIP.

Refresh from NCBI, or use your own downloaded ortholog ZIP.

```bash
python3 prepare_data.py
python3 prepare_data.py --package /path/to/brca1_orthologs.zip
```

Refreshing can change accession versions. Check `data/manifest.json` after refresh. The included snapshots and their checksums preserve this run.

## Pilot design

Default input is three consecutive, non-overlapping 600-base windows (CDS coordinates 0–600, 600–1200, 1200–1800) per sequence. Human is the query for each comparison. There are nine real comparisons and nine shuffled controls, using seed 691. Each shuffled target preserves its window's nucleotide composition.

```bash
python3 demo.py --window 300 --windows 3 --output results_300
```

Windows at matching offsets are not guaranteed to represent precisely homologous regions after insertions or deletions. Full-length comparisons and broader ortholog sampling are next steps. Windows within one BRCA1 sequence are not independent biological replicates. One shuffle per window is a diagnostic control and cannot establish a p-value. This pilot does not estimate a phylogeny or evolutionary distance.

## Algorithm and metrics

For query position i and target position j,

`H[i,j] = max(0, H[i-1,j-1] + substitution, H[i-1,j] - 2, H[i,j-1] - 2)`

Match +2, mismatch -1, linear gap -2 per base. N always receives mismatch. First row and column are zero. Traceback starts at the first maximum-scoring cell and ends at a zero cell. Ties resolve diagonal, then up, then left.

- Identity = matching unambiguous bases / all alignment columns, including gaps.
- Query/target coverage = bases consumed / respective input window length.
- Normalized score = raw score / (2 × shorter input length).
- Coordinates are zero-based and half-open, relative to each window. Add `cds_start_0based` to obtain CDS-relative coordinates.

Time and space are O(mn). Scores use an int array and traceback uses one byte per cell. On this machine the DP buffers occupy 5(m+1)(n+1) bytes. The program caps a pair at 50 million cells to keep this traceback prototype practical. Biopython is used only for independent validation.

## Measured outputs

`results/alignments.csv` contains every observation and runtime. `alignments.json` also includes alignment strings. `benchmark.json` describes input volume, cell updates, scoring, hardware/OS context, runtime and memory measurement scope. `conservation.png` and `benchmark.png` have corresponding PDF versions.

Runtime uses `time.perf_counter`. Memory uses the whole-process peak resident set size via `resource.ru_maxrss`, normalized to MiB on macOS/Linux. That is a high-water mark across the process, not a fresh allocation measurement for each pair. Benchmark metadata is captured before plotting imports. Network download and plotting time are excluded. Your Mac and HPC timings will differ.

## Scaling toward the final project

The checkpoint computes 18 pairs at 600 × 600 bases each. For one million pairs at the same lengths, multiply the measured mean seconds per pair by one million for a rough serial estimate. With p workers, divide that estimate by p only as an ideal lower bound. Scheduling, serialization, I/O and load imbalance add cost.

The per-pair DP buffers scale quadratically with sequence length. Doubling both lengths requires four times as many cells. Whole-process memory includes interpreter overhead, and each worker has its own buffers. Full CDS pairs of roughly 5,600 bases require about 150 MiB in DP buffers per active worker. Benchmark those before submitting large jobs.

A million comparisons is still a planned computational workload. The present results contain 18 measured comparisons. A scientifically meaningful expansion requires more orthologous genes and a defined sampling plan. Arbitrary all-against-all windows from a single gene do not provide a million independent biological observations.

## References

- NCBI gene ortholog packages: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/
- NCBI gene data packages: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-gene-data-package/
- Exact API download URL and sequence identifiers: `data/manifest.json`.
- Smith TF, Waterman MS. Identification of common molecular subsequences. J Mol Biol. 1981;147:195–197. https://doi.org/10.1016/0022-2836(81)90087-5

