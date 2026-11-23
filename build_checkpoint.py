"""Generate a checkpoint reflection and report from measured results."""
import argparse,csv,json,statistics
import reportlab
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak

ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--repository',default='https://github.com/Jafaga/brca1-conservation')
args=p.parse_args()
b=json.loads((ROOT/'results/benchmark.json').read_text())
m=json.loads((ROOT/'data/manifest.json').read_text())
rows=list(csv.DictReader((ROOT/'results/alignments.csv').open()))
pdfmetrics.registerFont(TTFont('DejaVu', str(Path(reportlab.__file__).parent / 'fonts/Vera.ttf')))
pdfmetrics.registerFont(TTFont('DejaVu-Bold', str(Path(reportlab.__file__).parent / 'fonts/VeraBd.ttf')))
styles=getSampleStyleSheet()
for style in styles.byName.values():
    style.fontName = 'DejaVu-Bold' if style.name in ['Title','Heading1','Heading2'] else 'DejaVu'
styles['BodyText'].fontSize=10
styles['BodyText'].leading=14
styles['BodyText'].spaceAfter=10
styles['Title'].textColor=colors.HexColor('#55338d')
story=[]
def para(text,style='BodyText'):
    story.append(Paragraph(text,styles[style]))
def table(data,widths):
    t=Table([[Paragraph(str(v),styles['BodyText']) for v in row] for row in data],colWidths=widths,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#ece6f5')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,-1),.4,colors.HexColor('#dad4e2')),('TOPPADDING',(0,0),(-1,-1),6)]))
    story.extend([t,Spacer(1,10)])
para('BRCA1 Sequence Conservation','Title')
para('ICS 691G - Mid-term Checkpoint | Working prototype and preliminary results')
para('Research question','Heading2')
para('Does human BRCA1 have higher sequence similarity to a closely related mammalian ortholog than to more distant mammalian orthologs under fixed Smith-Waterman scoring? The broader hypothesis predicts higher similarity for more closely related species. This pilot evaluates the computational workflow on one gene.')
para('Progress reflection','Heading2')
para('For this checkpoint, I have a working prototype that uses Smith-Waterman to compare BRCA1 coding sequences across species. The program finds a local alignment and reports its score, percent identity, length, matches, mismatches, gaps, and coverage. I included a short toy example that gives the expected ACGT alignment with a score of 8. The code also passed checks against Biopython for 100 short sequence pairs and all nine real-data comparisons. The demo uses four public NCBI sequences, with their accessions and lengths recorded so the input can be checked.')
para('One challenge is choosing comparable sequences because a gene can have several transcript variants and different lengths across species. For this first run, I used one transcript per species and three 600-base windows from each sequence. Based on the proposal feedback, I am measuring similarity using alignment scores, percent identity, and coverage. The windows provide a manageable starting point, although insertions and deletions can shift matching regions. My next step is to compare the full coding sequences and expand the set of orthologs.')
para('This small run gives me measured runtime and memory usage before moving to larger inputs. I still need to define the larger biological dataset, benchmark full-length sequences, and divide the work across HPC jobs. The million-comparison run is planned for the final stage. For now, the prototype can read real sequences, calculate local alignments, and generate figures that help me examine sequence conservation.')
para('Code repository','Heading2')
para(f'<link href="{escape(args.repository)}" color="#55338d">{escape(args.repository)}</link>')
para('Run python3 demo.py from the project folder, then open results/report.html. README.md includes setup, data preparation, a toy example, metric definitions, and a live demo walkthrough.')
story.append(PageBreak())
para('Data and computational benchmark','Title')
table([['Species','Gene','RefSeq transcript','CDS bases']]+[[s['species'],s['gene_id'],s['accession'],s['cds_length']] for s in m['sequences']],[145,85,155,70])
para('Human NM_007294.4 was selected explicitly. For other species, selection preferred curated NM_ transcripts, then longest CDS, with a lexical tie-break. Cat uses a predicted XM_ transcript. Data were downloaded afresh for this prototype. The source URL, actual timestamp, headers, and checksums are in data/manifest.json.')
para('The input subset uses CDS intervals [0,600), [600,1200), and [1200,1800) from each sequence. Human is compared with each of the other three species. Each target window is also shuffled once with seed 691, preserving base composition. Scoring is match +2, mismatch -1, gap -2 per base. N receives a mismatch score.')
table([['Measurement','Observed value'],['Selected complete CDS data',f"4 sequences; {b['selected_cds_bases']:,} bases"],['Unique real subset used',f"{b['unique_real_input_bases']:,} bases; 1,800 per species"],['Comparisons and DP work',f"18 pairs (9 real + 9 shuffled); {b['dp_cells_total']:,} cells"],['Alignment runtime',f"{b['alignment_runtime_seconds']:.3f} seconds"],['Pipeline runtime before plotting',f"{b['pipeline_runtime_seconds']:.3f} seconds"],['Peak process memory',f"{b['peak_process_rss_mib']:.2f} MiB"],['Largest per-pair DP buffers',f"{b['max_dp_buffer_bytes']/1048576:.2f} MiB"]],[180,250])
para(f"Environment: {b['platform']}, Python {b['python']}. Runtime uses time.perf_counter. Memory is process peak resident set size from resource.ru_maxrss, captured before plotting imports. DP buffer size is computed from allocated arrays. Download and plotting are excluded. These are execution-environment measurements; Mac and HPC costs still need measurement.")
mean=b['alignment_runtime_seconds']/b['comparisons']
para('Scale toward the final analysis','Heading2')
para(f"This subset covers 7,200 of 22,245 selected CDS bases (32.4%). The 18 pairs are 0.0018% of a planned million-comparison workload. At the same input lengths, the observed mean {mean:.4f} seconds per pair extrapolates to about {mean*1e6/3600:.1f} serial hours for one million pairs. This estimate is unmeasured and excludes scheduling and I/O. Both time and traceback storage scale as O(mn). Doubling both input lengths quadruples DP work. A roughly 5,600-by-5,600 pair needs about 150 MiB of DP buffers plus process overhead, per concurrent worker.")
story.append(PageBreak())
para('Preliminary sequence results','Title')
story.append(Image(str(ROOT/'results/conservation.png'),width=468,height=158.4))
para('Figure 1. Dots show individual windows; short bars show means. Normalized score divides score by 2 times the shorter input length. Identity divides matches by all alignment columns, including gaps. Coverage is the fraction of the human input window consumed.')
data=[['Comparison to human','Mean score','Identity','Coverage']]
for s in m['sequences'][1:]:
    a=[r for r in rows if r['species']==s['species'] and r['condition']=='real']
    data.append([s['species'],f"{statistics.mean(float(r['score']) for r in a):.1f}",f"{statistics.mean(float(r['percent_identity']) for r in a):.1f}%",f"{statistics.mean(float(r['query_coverage']) for r in a):.1f}%"])
table(data,[155,85,90,90])
para('Chimpanzee has the highest similarity in the selected windows, followed by cat and mouse. Each real comparison scores above its shuffled counterpart. This preliminary pattern supports using the prototype for a larger conservation analysis. Score, identity, and coverage must be interpreted together because a short local match can have high identity.')
para('Windows from one gene are not independent biological replicates. Matching offsets may miss regions shifted by indels. One shuffle per window is a diagnostic control. These results do not establish statistical significance, a general species relationship, or an evolutionary distance. Full-CDS comparisons, additional genes, and broader ortholog sampling remain next steps.')
story.append(PageBreak())
para('Preliminary performance figure','Title')
story.append(Image(str(ROOT/'results/benchmark.png'),width=468,height=187.2))
para('Figure 2. Measured runtime for the nine real comparisons, with the maximum allocated DP buffers and whole-process peak resident memory. All individual values are saved in results/alignments.csv. Peak RSS is a high-water mark across this process, rather than an isolated per-pair allocation.')
para('Reproducibility and live demonstration','Heading2')
para('Create a Python virtual environment and install requirements.txt. Run python3 -m unittest discover -s tests -v, then python3 demo.py. Open results/report.html. For the toy case run python3 demo.py --query TTACGTAA --target GGACGTCC. Show the recurrence and traceback in alignment.py. The exact source snapshot, fixed scoring, seed, checksums, result CSV, and traceback strings are included.')
para('Methodology changes based on earlier feedback','Heading2')
para('The study now uses specified orthologous coding sequences with exact lengths and public accessions. A hand-checkable toy dataset accompanies the real sequences. Conservation is evaluated using score, identity, and coverage. The planned large run will need additional biologically defined inputs and HPC benchmarking. Cat is an added mammalian comparison. Zebrafish was absent from the queried human BRCA1 ortholog set, so it was excluded from this pilot.')
para('References','Heading2')
para('NCBI Datasets. Gene and ortholog data packages. https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/')
para('Exact accessions, source download URL, retrieval timestamp, and checksums are provided in data/manifest.json.')
para('Smith TF, Waterman MS. Identification of common molecular subsequences. J Mol Biol. 1981;147:195-197. https://doi.org/10.1016/0022-2836(81)90087-5')
def footer(c,d):
    c.setFont('DejaVu',8)
    c.drawString(46,28,'ICS 691G | BRCA1 sequence conservation | Preliminary checkpoint results')
    c.drawRightString(566,28,str(d.page))
SimpleDocTemplate(str(ROOT/'checkpoint.pdf'),pagesize=(612,792),leftMargin=46,rightMargin=46,topMargin=36,bottomMargin=45,title='ICS 691G - BRCA1 Mid-term Checkpoint',author='Justine Afaga').build(story,onFirstPage=footer,onLaterPages=footer)
