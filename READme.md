# Needleman_Wusch_Visualisations
Interactive visualization of the Needleman–Wunsch global sequence alignment algorithm. Step-by-step scoring matrix, traceback path highlighting, and directional arrows — built with Python + Flask.

# 🧬 Needleman–Wunsch Visualizer

An interactive, browser-based visualization tool for the **Needleman–Wunsch global sequence alignment algorithm** — built for bioinformatics students and educators who want to *see* the dynamic programming matrix come to life.

---

## ✨ Features

- **Interactive scoring matrix** — fully computed and rendered as an annotated HTML table
- **Traceback path highlighting** — the optimal alignment path is visually distinct at a glance
- **Directional arrows per cell** — ↖ diagonal, ↑ up, ← left, color-coded by direction type
- **Parent-cell hover** — mouse over any cell to illuminate its three predecessor cells
- **Configurable parameters** — set your own match score, mismatch penalty, and gap penalty
- **Aligned sequence output** — final alignment displayed as color-coded character tiles (match / mismatch / gap)
- **Zero frontend dependencies** — pure HTML/CSS/JS; no Node, no npm, no bundler

---

## 🖥️ Demo

| Sequences | Gap | Match | Mismatch |
|-----------|-----|-------|----------|
| `CTTAACT` vs `CGGATCAT` | −3 | +1 | −1 |

The initialization row produces `0, −3, −6, −9, −12, −15, −18, −21, −24` — consistent with standard Needleman–Wunsch textbook examples.

---

## 🚀 Quick Start

### Prerequisites

- Python 3.10+
- Flask

```bash
pip install flask
```

### Run

```bash
python needleman_wunsch.py
```

The server starts on `http://localhost:5000` and opens your browser automatically.

---

## 🧪 Usage

1. Enter **Sequence 1** (vertical axis) and **Sequence 2** (horizontal axis) in the input fields.
2. Set your scoring parameters: **Match**, **Mismatch**, and **Gap** penalty.
3. Click **▶ Run Alignment**.
4. Explore the matrix:
   - **Purple-tinted cells** = optimal traceback path
   - **Blue-tinted cells** = initialization row/column
   - **Hover** any cell to highlight its parent cells in amber
5. Read the final aligned sequences at the bottom of the page.

> Sequences are capped at **30 characters** for readability. Longer sequences will warn before rendering.

---

## 📁 Project Structure

```
needleman_wunsch.py   ← Single-file Python app (server + embedded HTML/CSS/JS)
README.md
```

The entire front-end (HTML template, CSS design tokens, JavaScript rendering logic) is embedded as a string inside the Python file. No separate template directory needed.

---

## ⚙️ Algorithm Details

The implementation follows the standard three-step Needleman–Wunsch procedure:

**1. Initialization**
```
F(i, 0) = i × gap_penalty
F(0, j) = j × gap_penalty
```

**2. Recurrence (fill)**
```
F(i, j) = max(
    F(i−1, j−1) + s(a_i, b_j),   ← diagonal (match/mismatch)
    F(i−1, j)   + gap,             ← up (gap in sequence 2)
    F(i,   j−1) + gap              ← left (gap in sequence 1)
)
```

**3. Traceback**

Starting from `F(m, n)`, the path follows the highest-scoring predecessor back to `F(0, 0)`. When multiple predecessors share the maximum score, **all** directional arrows are rendered — but the traceback follows a single canonical path (diagonal preferred).

---

## 🎨 Tech Stack

| Layer | Technology |
|-------|-----------|
| Algorithm | Pure Python 3.10 |
| Server | Flask (single route + JSON API) |
| Frontend | Vanilla HTML5 / CSS3 / ES6 JS |
| Fonts | Space Mono · DM Sans (Google Fonts) |

---

## 🔭 Potential Extensions

- [ ] **Smith–Waterman** local alignment variant
- [ ] **Affine gap penalties** (gap open vs. gap extend)
- [ ] **Multiple traceback paths** rendered simultaneously
- [ ] **FASTA file input** support
- [ ] **Exportable matrix** as PNG or CSV
- [ ] **Step-by-step animation** mode (fill one cell at a time)

---

## 📄 License

MIT — free to use, modify, and redistribute.

---

## 🙏 Acknowledgements

Algorithm based on:
> Needleman, S.B. & Wunsch, C.D. (1970). *A general method applicable to the search for similarities in the amino acid sequence of two proteins.* Journal of Molecular Biology, 48(3), 443–453.