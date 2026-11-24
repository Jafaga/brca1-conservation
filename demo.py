"""Reproduce the checkpoint pilot and its measured results."""
import argparse
import csv
import json
import platform
import random
import resource
import sys
import time
from pathlib import Path
from alignment import read_fasta, smith_waterman

ROOT = Path(__file__).resolve().parent


def peak_rss_mib():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024 if sys.platform == 'darwin' else 1024)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT / 'data')
    parser.add_argument('--window', type=int, default=600)
    parser.add_argument('--windows', type=int, default=3)
    parser.add_argument('--output', type=Path, default=ROOT / 'results')
    parser.add_argument('--query', help='Align a custom DNA sequence with --target.')
    parser.add_argument('--target')
    args = parser.parse_args()
    if args.query is not None or args.target is not None:
        if args.query is None or args.target is None:
            parser.error('--query and --target must be supplied together')
        print(json.dumps(smith_waterman(args.query, args.target), indent=2))
        return
    if args.window < 1 or args.windows < 1:
        parser.error('Window size and count must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.data_dir / 'manifest.json').read_text())
    records = read_fasta(args.data_dir / manifest.get('fasta_file', 'brca1_selected.fna'))
    metadata = manifest['sequences']
    required = args.window * args.windows
    if required > min(len(sequence) for _, sequence in records):
        parser.error('Requested windows exceed the shortest selected CDS')
    started = time.perf_counter()
    rows, alignments = [], []
    rng = random.Random(691)
    human = records[0][1]
    for index in range(args.windows):
        start = index * args.window
        query = human[start:start + args.window]
        for species_index in range(1, len(records)):
            target = records[species_index][1][start:start + args.window]
            shuffled = list(target)
            rng.shuffle(shuffled)
            for condition, sequence in [('real', target), ('shuffled', ''.join(shuffled))]:
                tick = time.perf_counter()
                result = smith_waterman(query, sequence)
                elapsed = time.perf_counter() - tick
                row = dict(species=metadata[species_index]['species'],
                           accession=metadata[species_index]['accession'], condition=condition,
                           window=index + 1, cds_start_0based=start, query_length=len(query),
                           target_length=len(sequence), runtime_seconds=elapsed, process_peak_rss_mib=peak_rss_mib())
                row.update({k: v for k, v in result.items() if not k.startswith('aligned_')})
                row['normalized_score'] = result['score'] / (2 * min(len(query), len(sequence)))
                rows.append(row)
                alignments.append(dict(**row, aligned_query=result['aligned_query'], aligned_target=result['aligned_target']))
                print(f"{condition:8} window {index+1}  {row['species']:17} score {result['score']:4}  identity {result['percent_identity']:5.1f}%  {elapsed:.3f}s", flush=True)
    benchmark = dict(platform=platform.platform(), python=sys.version.split()[0],
        scoring=dict(match=2, mismatch=-1, gap=-2, ambiguous_N='mismatch'),
        window_length=args.window, windows_per_species=args.windows, comparisons=len(rows),
        real_comparisons=sum(r['condition'] == 'real' for r in rows),
        selected_sequences=len(records), selected_cds_bases=sum(len(s) for _, s in records),
        unique_real_input_bases=len(records) * required,
        pairwise_input_bases=sum(r['query_length'] + r['target_length'] for r in rows),
        dp_cells_total=sum(r['dp_cells'] for r in rows),
        alignment_runtime_seconds=sum(r['runtime_seconds'] for r in rows),
        pipeline_runtime_seconds=time.perf_counter() - started,
        peak_process_rss_mib=peak_rss_mib(),
        memory_method='resource.ru_maxrss, whole-process high-water mark; measured before importing plotting libraries',
        max_dp_buffer_bytes=max(r['dp_buffer_bytes'] for r in rows),
        scope=f'{args.windows} position-matched, non-overlapping CDS windows per species; within-gene windows are not biological replicates.',
        download_included=False, plotting_included=False, random_seed=691)
    with (args.output / 'alignments.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (args.output / 'alignments.json').write_text(json.dumps(alignments, indent=2) + '\n')
    (args.output / 'benchmark.json').write_text(json.dumps(benchmark, indent=2) + '\n')
    from figures import create_figures
    create_figures(args.output, rows, benchmark, manifest)
    report = args.output / 'report.html'
    report.write_text(report.read_text().replace('BRCA1', manifest.get('gene', 'BRCA1')))
    print(f"\nResults written to {args.output}\nOpen {args.output / 'report.html'} in your browser.")


if __name__ == '__main__':
    main()
