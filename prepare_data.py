"""Download four named BRCA1 genes and select reproducible CDS records."""
import argparse
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from pathlib import Path
from datetime import datetime, timezone
from alignment import read_fasta

URL = 'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/672,449497,12189,101081937/download?include_annotation_type=FASTA_CDS'
GENES = {'BRCA1': 672, 'TP53': 7157, 'BRCA2': 675}
SPECIES = ['Homo sapiens', 'Pan troglodytes', 'Mus musculus', 'Felis catus']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene', choices=GENES, default='BRCA1')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--package', type=Path, help='Use an existing NCBI Datasets ZIP, including your full ortholog package.')
    args = parser.parse_args()
    root = args.output or Path(__file__).resolve().parent / 'data'
    root.mkdir(parents=True, exist_ok=True)
    url = URL
    if not args.package and args.gene != 'BRCA1':
        ortholog_url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{GENES[args.gene]}/orthologs?taxon_filter=all'
        response = json.load(urllib.request.urlopen(ortholog_url, timeout=60))
        genes_by_species = {report['gene']['taxname']: report['gene']['gene_id'] for report in response['reports']}
        ids = ','.join(str(genes_by_species[species]) for species in SPECIES)
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{ids}/download?include_annotation_type=FASTA_CDS'
    root.mkdir(exist_ok=True)
    raw = args.package.read_bytes() if args.package else urllib.request.urlopen(url, timeout=60).read()
    archive = zipfile.ZipFile(io.BytesIO(raw))
    cds = archive.read('ncbi_dataset/data/cds.fna')
    reports = archive.read('ncbi_dataset/data/data_report.jsonl')
    (root / 'source_cds.fna').write_bytes(cds)
    records = read_fasta(root / 'source_cds.fna')
    genes = {str(row['geneId']): row for row in map(json.loads, reports.decode().splitlines())}
    human_genes = [row for row in genes.values() if row['taxname'] == 'Homo sapiens']
    if len(human_genes) != 1 or str(human_genes[0]['geneId']) != str(GENES[args.gene]):
        raise ValueError('Use the single-gene ortholog package matching --gene.')
    manifest, selected = [], []
    for species in SPECIES:
        candidates = [(header, seq) for header, seq in records if f'[organism={species}]' in header]
        if not candidates:
            raise ValueError(f'No coding sequence found for {species}')
        canonical = [r for r in candidates if r[0].startswith('NM_007294.4:')] if species == 'Homo sapiens' and args.gene == 'BRCA1' else []
        # Prefer curated NM records, then longest CDS, then header for deterministic ties.
        header, seq = sorted(canonical or candidates, key=lambda r: (not r[0].startswith('NM_'), -len(r[1]), r[0]))[0]
        gene_id = re.search(r'\[GeneID=(\d+)\]', header).group(1)
        accession = header.split()[0].split(':')[0]
        gene = genes[gene_id]
        selected.append(f'>{header}\n' + '\n'.join(seq[i:i+70] for i in range(0, len(seq), 70)))
        manifest.append(dict(species=species, accession=accession, gene_id=gene_id,
                             cds_length=len(seq), header=header, chromosomes=gene.get('chromosomes', []),
                             sha256=hashlib.sha256(seq.encode()).hexdigest(),
                             selection='Human NM_007294.4 when present; otherwise prefer NM_, longest CDS, lexical header tie-break.'))
    fasta_file = 'brca1_selected.fna' if args.gene == 'BRCA1' else 'selected_cds.fna'
    (root / fasta_file).write_text('\n'.join(selected) + '\n')
    (root / 'manifest.json').write_text(json.dumps(dict(gene=args.gene, fasta_file=fasta_file, source_url=url if not args.package else None,
        source_package_name=args.package.name if args.package else None,
        retrieved_utc=datetime.now(timezone.utc).isoformat(), package_sha256=hashlib.sha256(raw).hexdigest(),
        source_gene_records=len(genes), source_cds_records=len(records), sequences=manifest), indent=2) + '\n')
    (root / 'source_cds.fna').unlink()
    (root / 'source_package.zip').write_bytes(raw)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
