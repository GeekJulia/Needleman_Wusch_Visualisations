#!/usr/bin/env python3
"""
Needleman-Wunsch Global Sequence Alignment — Interactive Visualization
Requires: flask (pip install flask)
Run:      python needleman_wunsch.py
Then open: http://localhost:5000
"""

from flask import Flask, render_template_string, request, jsonify
import json

app = Flask(__name__)

# ──────────────────────────────────────────────
#  ALGORITHM CORE
# ──────────────────────────────────────────────

def needleman_wunsch(seq1: str, seq2: str, match: int, mismatch: int, gap: int):
    """
    Returns:
        matrix   – 2-D list of scores  (rows = seq1+1, cols = seq2+1)
        traceback– 2-D list of direction sets per cell {'D','U','L'}
        path     – list of (i,j) tuples on the optimal traceback path
        aligned  – (aligned_seq1, aligned_seq2) strings
    """
    m, n = len(seq1), len(seq2)
    INF = float('-inf')

    # Initialise
    matrix = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[set() for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        matrix[i][0] = i * gap
        traceback[i][0].add('U')
    for j in range(1, n + 1):
        matrix[0][j] = j * gap
        traceback[0][j].add('L')

    # Fill
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = matrix[i - 1][j - 1] + s
            up   = matrix[i - 1][j]     + gap
            left = matrix[i][j - 1]     + gap
            best = max(diag, up, left)
            matrix[i][j] = best
            if diag == best: traceback[i][j].add('D')
            if up   == best: traceback[i][j].add('U')
            if left == best: traceback[i][j].add('L')

    # Traceback
    path = []
    i, j = m, n
    while i > 0 or j > 0:
        path.append((i, j))
        dirs = traceback[i][j]
        if 'D' in dirs and i > 0 and j > 0:
            i -= 1; j -= 1
        elif 'U' in dirs and i > 0:
            i -= 1
        else:
            j -= 1
    path.append((0, 0))
    path_set = set(path)

    # Reconstruct alignment
    a1, a2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        dirs = traceback[i][j]
        if 'D' in dirs and i > 0 and j > 0:
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            i -= 1; j -= 1
        elif 'U' in dirs and i > 0:
            a1.append(seq1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(seq2[j - 1])
            j -= 1

    aligned1 = ''.join(reversed(a1))
    aligned2 = ''.join(reversed(a2))

    # Serialise traceback as lists for JSON
    tb_serial = [[list(cell) for cell in row] for row in traceback]

    return matrix, tb_serial, list(path_set), (aligned1, aligned2)


# ──────────────────────────────────────────────
#  HTML TEMPLATE
# ──────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>Needleman–Wunsch Visualizer</title>
<link href="https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap" rel="stylesheet"/>
<style>
  /* ── DESIGN TOKENS ── */
  :root {
    --bg:       #0d1117;
    --surface:  #161b22;
    --panel:    #1c2330;
    --border:   #30363d;
    --text:     #e6edf3;
    --muted:    #8b949e;
    --accent:   #58a6ff;
    --accent2:  #f78166;
    --match:    #3fb950;
    --path:     #d2a8ff;
    --path-bg:  rgba(210,168,255,0.18);
    --hover-bg: rgba(88,166,255,0.12);
    --diag:     #ffa657;
    --mono:     'Space Mono', monospace;
    --sans:     'DM Sans', sans-serif;
    --radius:   8px;
    --cell:     52px;
  }

  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    font-family: var(--sans);
    background: var(--bg);
    color: var(--text);
    min-height: 100vh;
    padding: 2rem 1.5rem 4rem;
  }

  /* ── HEADER ── */
  header {
    text-align: center;
    margin-bottom: 2.5rem;
    padding-bottom: 1.5rem;
    border-bottom: 1px solid var(--border);
  }
  header h1 {
    font-family: var(--mono);
    font-size: clamp(1.4rem, 3vw, 2.2rem);
    color: var(--accent);
    letter-spacing: -0.03em;
  }
  header p {
    color: var(--muted);
    font-size: 0.9rem;
    margin-top: 0.4rem;
  }

  /* ── LAYOUT ── */
  .layout {
    display: grid;
    grid-template-columns: 320px 1fr;
    gap: 1.5rem;
    max-width: 1400px;
    margin: 0 auto;
    align-items: start;
  }
  @media (max-width: 900px) {
    .layout { grid-template-columns: 1fr; }
  }

  /* ── PANEL ── */
  .panel {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 1.25rem 1.25rem 1.5rem;
  }
  .panel h2 {
    font-family: var(--mono);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    color: var(--muted);
    margin-bottom: 1rem;
  }

  /* ── FORM ── */
  .field { margin-bottom: 0.9rem; }
  label {
    display: block;
    font-size: 0.78rem;
    font-weight: 500;
    color: var(--muted);
    margin-bottom: 0.3rem;
    text-transform: uppercase;
    letter-spacing: 0.06em;
  }
  input[type=text], input[type=number] {
    width: 100%;
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--text);
    font-family: var(--mono);
    font-size: 0.95rem;
    padding: 0.55rem 0.75rem;
    outline: none;
    transition: border-color 0.2s;
  }
  input:focus { border-color: var(--accent); }
  input[type=text] { text-transform: uppercase; letter-spacing: 0.15em; }

  .score-row {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 0.6rem;
  }

  .btn-run {
    width: 100%;
    margin-top: 1.2rem;
    padding: 0.7rem;
    background: var(--accent);
    color: #000;
    border: none;
    border-radius: 6px;
    font-family: var(--mono);
    font-size: 0.9rem;
    font-weight: 700;
    cursor: pointer;
    letter-spacing: 0.05em;
    transition: opacity 0.2s, transform 0.1s;
  }
  .btn-run:hover { opacity: 0.88; }
  .btn-run:active { transform: scale(0.97); }

  /* ── LEGEND ── */
  .legend {
    margin-top: 1.4rem;
    display: flex;
    flex-direction: column;
    gap: 0.45rem;
  }
  .legend-item {
    display: flex;
    align-items: center;
    gap: 0.55rem;
    font-size: 0.78rem;
    color: var(--muted);
  }
  .legend-swatch {
    width: 18px; height: 18px;
    border-radius: 4px;
    flex-shrink: 0;
  }

  /* ── SCORE DISPLAY ── */
  .score-display {
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.7rem 1rem;
    margin-top: 1rem;
    font-family: var(--mono);
    font-size: 0.85rem;
  }
  .score-display span { color: var(--path); font-weight: 700; }

  /* ── MAIN AREA ── */
  .main-area { display: flex; flex-direction: column; gap: 1.5rem; }

  /* ── MATRIX WRAPPER ── */
  .matrix-wrapper {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 1.25rem;
    overflow-x: auto;
  }
  .matrix-wrapper h2 {
    font-family: var(--mono);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    color: var(--muted);
    margin-bottom: 1.2rem;
  }

  /* ── TABLE ── */
  table {
    border-collapse: collapse;
    font-family: var(--mono);
    font-size: 0.82rem;
  }
  th, td {
    width: var(--cell);
    height: var(--cell);
    text-align: center;
    vertical-align: middle;
    border: 1px solid var(--border);
    position: relative;
    user-select: none;
  }
  th {
    background: var(--panel);
    color: var(--accent);
    font-weight: 700;
    font-size: 0.9rem;
    letter-spacing: 0.08em;
    border-color: #3d444d;
  }
  th.seq-label { color: var(--accent2); font-size: 1rem; }
  th.empty { background: var(--bg); border-color: transparent; }

  td { background: var(--surface); cursor: pointer; transition: background 0.15s; }
  td:hover { background: var(--hover-bg) !important; }

  td.path-cell { background: var(--path-bg) !important; }
  td.path-cell .score-val { color: var(--path); font-weight: 700; }

  td.init-cell { background: rgba(88,166,255,0.06); }
  td.init-cell .score-val { color: var(--accent); }

  td.highlight-parent { background: rgba(255,166,87,0.2) !important; outline: 2px solid var(--diag); outline-offset: -2px; }

  /* ── CELL CONTENT ── */
  .cell-inner {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    width: 100%; height: 100%;
    gap: 2px;
  }
  .score-val {
    font-size: 0.88rem;
    font-weight: 700;
    line-height: 1;
  }
  .arrows {
    font-size: 0.7rem;
    color: var(--muted);
    line-height: 1;
    letter-spacing: -0.02em;
  }
  .arrow-D { color: var(--diag); }
  .arrow-U { color: var(--match); }
  .arrow-L { color: var(--accent2); }

  /* ── ALIGNMENT DISPLAY ── */
  .alignment-panel {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 1.25rem;
  }
  .alignment-panel h2 {
    font-family: var(--mono);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    color: var(--muted);
    margin-bottom: 1.1rem;
  }
  .alignment-seqs {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
  }
  .seq-row {
    display: flex;
    gap: 0.7rem;
    align-items: center;
    flex-wrap: wrap;
  }
  .seq-label-tag {
    font-family: var(--mono);
    font-size: 0.7rem;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.1em;
    min-width: 28px;
  }
  .seq-chars { display: flex; gap: 2px; flex-wrap: wrap; }
  .seq-char {
    width: 28px; height: 28px;
    display: flex; align-items: center; justify-content: center;
    font-family: var(--mono);
    font-size: 0.85rem;
    font-weight: 700;
    border-radius: 4px;
    border: 1px solid var(--border);
  }
  .seq-char.match   { background: rgba(63,185,80,0.2);  color: var(--match);  border-color: var(--match); }
  .seq-char.gap     { background: rgba(247,129,102,0.2); color: var(--accent2); border-color: var(--accent2); }
  .seq-char.mismatch{ background: rgba(255,166,87,0.2);  color: var(--diag);   border-color: var(--diag); }

  .match-bar { display: flex; gap: 2px; flex-wrap: wrap; }
  .match-sym {
    width: 28px; height: 16px;
    display: flex; align-items: center; justify-content: center;
    font-family: var(--mono);
    font-size: 0.75rem;
    color: var(--match);
  }
  .match-sym.mismatch-sym { color: var(--accent2); }

  /* ── TOOLTIP / STATUS ── */
  #cell-info {
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.6rem 0.9rem;
    font-size: 0.78rem;
    color: var(--muted);
    margin-top: 0.8rem;
    min-height: 2.2rem;
    font-family: var(--mono);
    display: none;
  }
  #cell-info.visible { display: block; }
  #cell-info em { color: var(--text); font-style: normal; }

  /* ── PLACEHOLDER ── */
  .placeholder {
    display: flex;
    align-items: center;
    justify-content: center;
    height: 200px;
    color: var(--muted);
    font-size: 0.85rem;
    flex-direction: column;
    gap: 0.5rem;
  }
  .placeholder svg { opacity: 0.3; }

  /* ── LOADING ── */
  #loading {
    display: none;
    text-align: center;
    padding: 2rem;
    color: var(--muted);
    font-family: var(--mono);
    font-size: 0.85rem;
  }
  .spinner {
    display: inline-block;
    width: 20px; height: 20px;
    border: 2px solid var(--border);
    border-top-color: var(--accent);
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
    margin-right: 0.5rem;
    vertical-align: middle;
  }
  @keyframes spin { to { transform: rotate(360deg); } }

  /* ── FADE-IN ── */
  .fade-in { animation: fadeIn 0.35s ease forwards; }
  @keyframes fadeIn { from { opacity:0; transform: translateY(8px); } to { opacity:1; transform: none; } }
</style>
</head>
<body>

<header>
  <h1>Needleman–Wunsch</h1>
  <p>Global Sequence Alignment — Interactive Matrix Visualization</p>
</header>

<div class="layout">

  <!-- ── CONTROLS PANEL ── -->
  <aside class="panel">
    <h2>Parameters</h2>

    <div class="field">
      <label>Sequence 1 (vertical)</label>
      <input type="text" id="seq1" value="CTTAACT" placeholder="e.g. CTTAACT"/>
    </div>
    <div class="field">
      <label>Sequence 2 (horizontal)</label>
      <input type="text" id="seq2" value="CGGATCAT" placeholder="e.g. CGGATCAT"/>
    </div>

    <div class="score-row">
      <div class="field">
        <label>Match</label>
        <input type="number" id="match" value="1"/>
      </div>
      <div class="field">
        <label>Mismatch</label>
        <input type="number" id="mismatch" value="-1"/>
      </div>
      <div class="field">
        <label>Gap</label>
        <input type="number" id="gap" value="-3"/>
      </div>
    </div>

    <button class="btn-run" onclick="runAlignment()">▶ Run Alignment</button>

    <div id="score-result" class="score-display" style="display:none">
      Optimal score: <span id="optimal-score">—</span>
    </div>

    <!-- Legend -->
    <div class="legend">
      <div class="legend-item">
        <div class="legend-swatch" style="background:var(--path-bg); border:1px solid var(--path)"></div>
        Traceback path
      </div>
      <div class="legend-item">
        <div class="legend-swatch" style="background:rgba(88,166,255,0.06); border:1px solid var(--accent)"></div>
        Initialisation row/col
      </div>
      <div class="legend-item">
        <div class="legend-swatch" style="background:rgba(255,166,87,0.2); border:1px solid var(--diag)"></div>
        Highlighted parent cells
      </div>
      <div class="legend-item">
        <div style="font-family:var(--mono); color:var(--diag); font-size:1rem">↖</div>
        Diagonal (match/mismatch)
      </div>
      <div class="legend-item">
        <div style="font-family:var(--mono); color:var(--match); font-size:1rem">↑</div>
        Up (gap in seq2)
      </div>
      <div class="legend-item">
        <div style="font-family:var(--mono); color:var(--accent2); font-size:1rem">←</div>
        Left (gap in seq1)
      </div>
    </div>
  </aside>

  <!-- ── MAIN AREA ── -->
  <main class="main-area">

    <!-- Matrix -->
    <div class="matrix-wrapper">
      <h2>Scoring Matrix</h2>
      <div id="loading"><span class="spinner"></span>Computing…</div>
      <div id="matrix-container">
        <div class="placeholder">
          <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1">
            <rect x="3" y="3" width="7" height="7"/><rect x="14" y="3" width="7" height="7"/>
            <rect x="14" y="14" width="7" height="7"/><rect x="3" y="14" width="7" height="7"/>
          </svg>
          Enter sequences and click Run
        </div>
      </div>
      <div id="cell-info"></div>
    </div>

    <!-- Alignment result -->
    <div class="alignment-panel" id="alignment-panel" style="display:none">
      <h2>Optimal Alignment</h2>
      <div class="alignment-seqs" id="alignment-display"></div>
    </div>

  </main>
</div>

<script>
let lastData = null;

async function runAlignment() {
  const seq1    = document.getElementById('seq1').value.trim().toUpperCase();
  const seq2    = document.getElementById('seq2').value.trim().toUpperCase();
  const match   = parseInt(document.getElementById('match').value);
  const mismatch= parseInt(document.getElementById('mismatch').value);
  const gap     = parseInt(document.getElementById('gap').value);

  if (!seq1 || !seq2) { alert('Please enter both sequences.'); return; }
  if (seq1.length > 30 || seq2.length > 30) { alert('For readability, keep sequences ≤ 30 characters.'); return; }

  document.getElementById('loading').style.display = 'block';
  document.getElementById('matrix-container').innerHTML = '';
  document.getElementById('alignment-panel').style.display = 'none';
  document.getElementById('score-result').style.display = 'none';

  try {
    const resp = await fetch('/align', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ seq1, seq2, match, mismatch, gap })
    });
    const data = await resp.json();
    lastData = data;
    renderMatrix(data);
    renderAlignment(data);
    document.getElementById('optimal-score').textContent = data.optimal_score;
    document.getElementById('score-result').style.display = 'block';
  } catch(e) {
    document.getElementById('matrix-container').innerHTML =
      '<div class="placeholder">Error: ' + e.message + '</div>';
  } finally {
    document.getElementById('loading').style.display = 'none';
  }
}

function renderMatrix(data) {
  const { seq1, seq2, matrix, traceback, path } = data;
  const pathSet = new Set(path.map(([r,c]) => `${r},${c}`));

  const m = seq1.length, n = seq2.length;

  let html = '<table class="fade-in">';

  // Header row: empty + empty + seq2 chars
  html += '<tr>';
  html += '<th class="empty"></th>';
  html += '<th class="empty"></th>';
  html += '<th class="empty">—</th>';
  for (let j = 0; j < n; j++) {
    html += `<th class="seq-label">${seq2[j]}</th>`;
  }
  html += '</tr>';

  // Data rows
  for (let i = 0; i <= m; i++) {
    html += '<tr>';
    // Seq1 label
    if (i === 0) {
      html += '<th class="empty"></th>';
      html += '<th class="empty">—</th>';
    } else {
      html += '<th class="empty"></th>';
      html += `<th class="seq-label">${seq1[i-1]}</th>`;
    }

    for (let j = 0; j <= n; j++) {
      const score = matrix[i][j];
      const isPath = pathSet.has(`${i},${j}`);
      const isInit = (i === 0 || j === 0);
      const dirs   = traceback[i][j] || [];

      let cls = 'data-cell';
      if (isPath) cls += ' path-cell';
      else if (isInit) cls += ' init-cell';

      const arrowHtml = buildArrows(dirs);

      html += `<td class="${cls}"
                   data-i="${i}" data-j="${j}"
                   data-score="${score}"
                   data-dirs="${dirs.join(',')}"
                   onmouseenter="onCellHover(this)"
                   onmouseleave="clearHover()">
                 <div class="cell-inner">
                   <span class="score-val">${score}</span>
                   <span class="arrows">${arrowHtml}</span>
                 </div>
               </td>`;
    }
    html += '</tr>';
  }
  html += '</table>';

  document.getElementById('matrix-container').innerHTML = html;
}

function buildArrows(dirs) {
  let parts = [];
  if (dirs.includes('D')) parts.push('<span class="arrow-D">↖</span>');
  if (dirs.includes('U')) parts.push('<span class="arrow-U">↑</span>');
  if (dirs.includes('L')) parts.push('<span class="arrow-L">←</span>');
  return parts.join('');
}

function onCellHover(td) {
  const i = parseInt(td.dataset.i);
  const j = parseInt(td.dataset.j);
  const score = td.dataset.score;
  const dirs  = td.dataset.dirs ? td.dataset.dirs.split(',').filter(Boolean) : [];

  // Highlight parents
  clearHover();
  const parents = [];
  if (dirs.includes('D') && i > 0 && j > 0) parents.push([i-1, j-1]);
  if (dirs.includes('U') && i > 0)           parents.push([i-1, j]);
  if (dirs.includes('L') && j > 0)           parents.push([i, j-1]);

  parents.forEach(([pi, pj]) => {
    const el = document.querySelector(`td[data-i="${pi}"][data-j="${pj}"]`);
    if (el) el.classList.add('highlight-parent');
  });

  // Info bar
  let info = `Cell (${i},${j}) — score: <em>${score}</em>`;
  if (i > 0 || j > 0) {
    const dirNames = dirs.map(d => d === 'D' ? '↖ diagonal' : d === 'U' ? '↑ up' : '← left');
    info += ` — arrived via: <em>${dirNames.join(', ') || '—'}</em>`;
  }
  const el = document.getElementById('cell-info');
  el.innerHTML = info;
  el.classList.add('visible');
}

function clearHover() {
  document.querySelectorAll('.highlight-parent').forEach(el => el.classList.remove('highlight-parent'));
  const el = document.getElementById('cell-info');
  el.classList.remove('visible');
}

function renderAlignment(data) {
  const { aligned1, aligned2, seq1, seq2, match_score, mismatch_penalty } = data;
  const len = aligned1.length;

  let charHtml1 = '', charHtml2 = '', barHtml = '';

  for (let k = 0; k < len; k++) {
    const c1 = aligned1[k], c2 = aligned2[k];
    let cls1 = '', cls2 = '', symCls = '', sym = '·';

    if (c1 === '-' || c2 === '-') {
      cls1 = cls2 = 'gap'; sym = ' '; symCls = 'mismatch-sym';
    } else if (c1 === c2) {
      cls1 = cls2 = 'match'; sym = '|';
    } else {
      cls1 = cls2 = 'mismatch'; sym = '✕'; symCls = 'mismatch-sym';
    }

    charHtml1 += `<span class="seq-char ${cls1}">${c1}</span>`;
    charHtml2 += `<span class="seq-char ${cls2}">${c2}</span>`;
    barHtml   += `<span class="match-sym ${symCls}">${sym}</span>`;
  }

  document.getElementById('alignment-display').innerHTML = `
    <div class="seq-row">
      <span class="seq-label-tag">S1</span>
      <div class="seq-chars">${charHtml1}</div>
    </div>
    <div class="seq-row" style="padding-left:calc(28px + 0.7rem)">
      <div class="match-bar">${barHtml}</div>
    </div>
    <div class="seq-row">
      <span class="seq-label-tag">S2</span>
      <div class="seq-chars">${charHtml2}</div>
    </div>
  `;

  document.getElementById('alignment-panel').style.display = 'block';
  document.getElementById('alignment-panel').classList.add('fade-in');
}

// Run on page load with default values
window.addEventListener('DOMContentLoaded', () => runAlignment());
</script>
</body>
</html>
"""


# ──────────────────────────────────────────────
#  ROUTES
# ──────────────────────────────────────────────

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)


@app.route('/align', methods=['POST'])
def align():
    body     = request.get_json(force=True)
    seq1     = body.get('seq1', '').upper().strip()
    seq2     = body.get('seq2', '').upper().strip()
    match    = int(body.get('match', 1))
    mismatch = int(body.get('mismatch', -1))
    gap      = int(body.get('gap', -2))

    matrix, traceback, path, (aligned1, aligned2) = needleman_wunsch(
        seq1, seq2, match, mismatch, gap
    )

    return jsonify({
        'seq1':            seq1,
        'seq2':            seq2,
        'matrix':          matrix,
        'traceback':       traceback,
        'path':            path,
        'aligned1':        aligned1,
        'aligned2':        aligned2,
        'optimal_score':   matrix[len(seq1)][len(seq2)],
        'match_score':     match,
        'mismatch_penalty':mismatch,
    })


# ──────────────────────────────────────────────
#  ENTRY POINT
# ──────────────────────────────────────────────

if __name__ == '__main__':
    import webbrowser, threading, time

    def open_browser():
        time.sleep(0.9)
        webbrowser.open('http://localhost:5000')

    threading.Thread(target=open_browser, daemon=True).start()
    print("\n  Needleman–Wunsch Visualizer")
    print("  ─────────────────────────────")
    print("  Open: http://localhost:5000")
    print("  Press Ctrl+C to quit\n")
    app.run(host='0.0.0.0', port=5000, debug=False)