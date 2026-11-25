"""Plot measured alignment and benchmark results."""
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
