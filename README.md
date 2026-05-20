# DNA Mutation Analyser

A locally-run Flask web application that takes two DNA sequences (reference and mutant) along with a protein name, detects mutations, translates them to protein level, fetches the real 3D structure from PDB, and displays everything on an interactive results page.

---

## Features

- **Sequence Alignment** — Needleman-Wunsch global alignment via BioPython
- **Mutation Detection** — substitutions, insertions, and deletions with exact positions
- **Protein Translation** — DNA → Protein using the genetic code (BioPython Seq)
- **3D Structure Viewer** — Interactive py3Dmol viewer with mutation site highlighted in red
- **UniProt Info Panel** — Live protein function, organism, and length from UniProt REST API
- **Mutation Chart** — Plotly scatter chart showing mutation positions across the sequence

---

## Project Structure

```
dna_mutation_analyser/
├── app.py                  Flask app — routes and controller logic
├── aligner.py              SequenceAligner class (Needleman-Wunsch + mutation detection)
├── translator.py           ProteinTranslator class (DNA → Protein + AA comparison)
├── pdb_fetcher.py          PDBFetcher class (RCSB PDB search + structure download)
├── uniprot_fetcher.py      UniProtFetcher class (UniProt annotation retrieval)
├── templates/
│   ├── index.html          Input form page
│   └── results.html        Results page (3D viewer, mutation table, charts)
├── requirements.txt        Python dependencies
├── .gitignore
└── README.md
```

---

## Setup and Run

### 1. Clone the repository
```bash
git clone https://github.com/YOUR_USERNAME/dna-mutation-analyser.git
cd dna-mutation-analyser
```

### 2. Create a virtual environment
```bash
python -m venv venv
venv\Scripts\activate        # Windows
source venv/bin/activate     # macOS / Linux
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

### 4. Run the app
```bash
python app.py
```

Open your browser at **http://localhost:5000**

---

## Example Input

| Field | Value |
|-------|-------|
| Reference DNA | `ATGCGTAAAGGTCCAATCGTAGCTGAAGATCGTTTACGC` |
| Mutant DNA | `ATGCATAAAGGTCCAATCGTAGCTGAAGATCGTTTACGC` |
| Protein name | `HLA-C` |

This detects a G→A substitution at position 6, resulting in an Arg→His amino acid change.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| Flask | Web framework |
| BioPython | Sequence alignment and translation |
| Requests | HTTP calls to PDB and UniProt APIs |
| Plotly | Mutation position chart |
| Pandas | Data handling |

---

## Notes

- PDB and UniProt fetching require an internet connection. All other processing (alignment, translation) works offline.
- Runs on a single machine at `localhost:5000`. Not intended for multi-user deployment.
- No database — results exist only for the session duration unless saved manually.
