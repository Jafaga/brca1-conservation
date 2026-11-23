const toyButton = document.querySelector('#toy-button');
const realButton = document.querySelector('#real-button');
const status = document.querySelector('#status');
const names = {'Pan troglodytes': 'Chimpanzee', 'Felis catus': 'Cat', 'Mus musculus': 'Mouse'};

async function request(path) {
  const response = await fetch(path, {method: 'POST'});
  const result = await response.json();
  if (!response.ok) throw new Error(result.error || 'The run could not finish. See Terminal for details.');
  return result;
}

toyButton.addEventListener('click', async () => {
  toyButton.disabled = true;
  const output = document.querySelector('#toy-result');
  try {
    const result = await request('/api/toy');
    output.className = 'toy-output';
    output.innerHTML = `<h3>The matching piece is <code>${result.aligned_query}</code></h3><p>4 matches × 2 points = <strong>${result.score} points</strong>. All four aligned letters match, so identity is <strong>${result.percent_identity}%</strong>. The alignment uses 4 of the 8 letters in each input, so coverage is <strong>${result.query_coverage}%</strong>.</p><p>Smith–Waterman finds the best matching local region and can leave the ends out.</p>`;
  } catch (error) {
    output.textContent = error.message;
  } finally {
    output.hidden = false;
    toyButton.disabled = false;
  }
});

realButton.addEventListener('click', async () => {
  realButton.disabled = true;
  realButton.textContent = 'Calculating alignments…';
  status.textContent = 'Running Python on the bundled NCBI sequences. This usually takes a few seconds. Keep Terminal open.';
  try {
    const result = await request('/api/pilot');
    const benchmark = result.benchmark;
    document.querySelector('#benchmark').innerHTML = `<div class="metric"><strong>${benchmark.comparisons}</strong><span>alignments calculated</span></div><div class="metric"><strong>${benchmark.alignment_runtime_seconds.toFixed(3)} s</strong><span>alignment runtime on this device</span></div><div class="metric"><strong>${benchmark.peak_process_rss_mib.toFixed(2)} MiB</strong><span>peak Python process memory</span></div>`;
    const mean = (rows, key) => rows.reduce((sum, row) => sum + row[key], 0) / rows.length;
    document.querySelector('#summary').innerHTML = Object.entries(names).map(([species, name]) => {
      const rows = result.alignments.filter(row => row.species === species && row.condition === 'real');
      return `<tr><td><strong>${name}</strong><small>${species}</small></td><td>${mean(rows, 'score').toFixed(1)}</td><td>${mean(rows, 'percent_identity').toFixed(1)}%</td><td>${mean(rows, 'query_coverage').toFixed(1)}%</td></tr>`;
    }).join('');
    document.querySelector('#manifest').innerHTML = result.manifest.sequences.map(sequence => `<tr><td>${sequence.species}</td><td>${sequence.accession}</td><td>${sequence.cds_length.toLocaleString()} bases</td></tr>`).join('');
    const alignment = result.alignments[0];
    const query = alignment.aligned_query.slice(0, 80);
    const target = alignment.aligned_target.slice(0, 80);
    const marks = [...query].map((base, index) => base === target[index] && base !== '-' ? '|' : ' ').join('');
    document.querySelector('#alignment').textContent = `Human   ${query}\n        ${marks}\nChimp   ${target}`;
    document.querySelector('#figure').src = `data:image/png;base64,${result.figure}`;
    document.querySelector('#real-result').hidden = false;
    status.textContent = 'Run complete. Scores were recalculated with the project’s Python implementation. Runtime and memory are measured before plotting.';
  } catch (error) {
    status.textContent = error.message;
  } finally {
    realButton.disabled = false;
    realButton.textContent = 'Run the real DNA demo';
  }
});
