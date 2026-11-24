"""Figures use measured CSV results. No invented observations."""
import html
from pathlib import Path
import statistics


def create_figures(output, rows, benchmark, manifest):
    gene = manifest.get('gene', 'BRCA1')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    species = list(dict.fromkeys(row['species'] for row in rows))
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
    metrics = [('normalized_score', 'Score / maximum possible score'),
               ('percent_identity', 'Identity (% of alignment columns)'),
               ('query_coverage', 'Human window coverage (%)')]
    for axis, (metric, label) in zip(axes, metrics):
        for condition, color, delta in [('real', '#6542b5', -.09), ('shuffled', '#b37724', .09)]:
            means = []
            for i, name in enumerate(species):
                values = [r[metric] for r in rows if r['species'] == name and r['condition'] == condition]
                axis.scatter([i + delta] * len(values), values, color=color, alpha=.6, s=28)
                means.append(statistics.mean(values))
            axis.plot([i+delta for i in range(len(species))], means, marker='_', markersize=18,
                      linestyle='none', color=color, label=condition)
        axis.set_xticks(range(len(species)), [s.replace(' ', '\n') for s in species])
        axis.set_ylabel(label)
        axis.grid(axis='y', alpha=.15)
    axes[0].legend(frameon=False)
    fig.suptitle(f'{gene} pilot: human CDS windows against three mammalian orthologs')
    fig.text(.5, .015, 'Each dot is one window. Bars mark means. Windows from one gene are not biological replicates.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0,.05,1,.94))
    fig.savefig(output / 'conservation.png', dpi=190)
    fig.savefig(output / 'conservation.pdf')
    plt.close(fig)
    real = [r for r in rows if r['condition'] == 'real']
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    labels = [f"{r['species'].split()[-1]} W{r['window']}" for r in real]
    axes[0].bar(range(len(real)), [r['runtime_seconds'] for r in real], color='#6542b5')
    axes[0].set_ylabel('Alignment runtime (seconds)')
    axes[0].set_xticks(range(len(real)), labels, rotation=50, ha='right')
    axes[1].bar(['DP buffers', 'Peak process RSS'], [benchmark['max_dp_buffer_bytes']/1048576, benchmark['peak_process_rss_mib']], color=['#6542b5','#b37724'])
    axes[1].set_ylabel('Memory (MiB)')
    fig.suptitle('Measured preliminary computational cost')
    fig.tight_layout()
    fig.savefig(output / 'benchmark.png', dpi=190)
    fig.savefig(output / 'benchmark.pdf')
    plt.close(fig)
    summary = []
    for name in species:
        records = [r for r in real if r['species'] == name]
        summary.append(f"<tr><td>{html.escape(name)}</td><td>{statistics.mean(r['score'] for r in records):.1f}</td><td>{statistics.mean(r['percent_identity'] for r in records):.1f}%</td><td>{statistics.mean(r['query_coverage'] for r in records):.1f}%</td></tr>")
    inventory = ''.join(f"<tr><td>{s['species']}</td><td>{s['accession']}</td><td>{s['gene_id']}</td><td>{s['cds_length']:,}</td></tr>" for s in manifest['sequences'])
    (output / 'report.html').write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>BRCA1 conservation pilot</title>
<style>body{{font:17px/1.6 system-ui;max-width:1100px;margin:50px auto;padding:0 24px;color:#272333;background:#faf9fc}}h1{{font-size:42px;line-height:1.1}}h2{{margin-top:38px}}table{{border-collapse:collapse;width:100%;background:white}}td,th{{padding:12px;text-align:left;border-bottom:1px solid #ddd}}img{{max-width:100%;background:white}}a{{color:#6542b5}}.meta{{color:#61596c}}code{{background:#ece7f5;padding:3px}}</style>
<p class="meta">ICS 691G · Working computational prototype</p><h1>How conserved is BRCA1 across species?</h1>
<p>Hypothesis for the broader study: more closely related species will tend to have higher BRCA1 sequence similarity to human under fixed alignment scoring. This pilot compares three mammalian orthologs using Smith–Waterman implemented from scratch.</p>
<h2>Real sequence input</h2><table><tr><th>Species</th><th>RefSeq transcript</th><th>NCBI Gene</th><th>Full CDS bases</th></tr>{inventory}</table>
<p>The pilot uses the first {benchmark['windows_per_species']} non-overlapping {benchmark['window_length']}-base windows from each CDS. Same-offset windows can miss homologous regions shifted by insertions or deletions. Broader biological conclusions need full-length alignment and more orthologs.</p>
<h2>Preliminary sequence results</h2><img src="conservation.png" alt="Measured BRCA1 alignment similarity, identity and coverage">
<table><tr><th>Species vs human</th><th>Mean score</th><th>Mean identity</th><th>Mean human window coverage</th></tr>{''.join(summary)}</table>
<p>Identity is matches divided by all alignment columns, including gaps. Coverage is the consumed portion of the input window. Shuffled controls preserve each target window's base composition, with seed 691. One shuffle per window provides a demonstration control. These observations do not estimate an evolutionary distance or statistical significance.</p>
<h2>Measured benchmark</h2><p>{benchmark['comparisons']} alignments ({benchmark['real_comparisons']} real and {benchmark['comparisons']-benchmark['real_comparisons']} shuffled), {benchmark['unique_real_input_bases']:,} unique real input bases, {benchmark['dp_cells_total']:,} DP cells. Alignment runtime {benchmark['alignment_runtime_seconds']:.3f} seconds. Peak process memory {benchmark['peak_process_rss_mib']:.2f} MiB. Measurements come from {html.escape(benchmark['platform'])}, Python {benchmark['python']}. Download and plotting are excluded.</p><img src="benchmark.png" alt="Measured runtime and memory">
<p><a href="alignments.csv">Download all observations (CSV)</a> · <a href="benchmark.json">Benchmark metadata</a> · <a href="alignments.json">Traceback alignments</a></p>
<h2>Run the demo again</h2><p>From the project directory run <code>python3 demo.py</code>. For a hand-checkable example run <code>python3 demo.py --query TTACGTAA --target GGACGTCC</code>.</p>
</html>''')
