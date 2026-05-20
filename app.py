"""
DNA Mutation Analyser — Flask Application
==========================================
Main controller that wires together all modules:
  - SequenceAligner  : Needleman-Wunsch alignment + mutation detection
  - ProteinTranslator: DNA → Protein translation + AA comparison
  - PDBFetcher       : RCSB PDB 3D structure retrieval
  - UniProtFetcher   : UniProt protein annotation retrieval

Run:
    python app.py
Then open http://localhost:5000 in your browser.
"""

from flask import Flask, render_template, request
import json

from aligner import SequenceAligner
from translator import ProteinTranslator
from pdb_fetcher import PDBFetcher
from uniprot_fetcher import UniProtFetcher

app = Flask(__name__)
app.secret_key = "dna_analyser_secret_key"


# ---------------------------------------------------------------------------
# Input Validation
# ---------------------------------------------------------------------------

def validate_inputs(ref_dna: str, mut_dna: str, protein_name: str) -> list:
    """
    Validates the three form inputs.
    Returns a list of error strings (empty list = all valid).
    """
    errors = []
    valid_bases = set("ATGCN-")

    if not ref_dna:
        errors.append("Reference DNA sequence is required.")
    elif not set(ref_dna.upper()).issubset(valid_bases):
        errors.append("Reference DNA contains invalid characters. Use only A, T, G, C.")

    if not mut_dna:
        errors.append("Mutant DNA sequence is required.")
    elif not set(mut_dna.upper()).issubset(valid_bases):
        errors.append("Mutant DNA contains invalid characters. Use only A, T, G, C.")

    if not protein_name:
        errors.append("Protein name is required.")

    if ref_dna and mut_dna and ref_dna.upper() == mut_dna.upper():
        errors.append("Reference and mutant sequences are identical — no mutation to detect.")

    return errors


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    """Renders the input form (Page 1)."""
    return render_template("index.html")


@app.route("/submit", methods=["POST"])
def submit():
    """
    Handles form submission.
    Runs the full analysis pipeline and renders the results page (Page 2).
    Returns to the input form with errors if validation fails at any step.
    """
    ref_dna      = request.form.get("ref_dna", "").strip()
    mut_dna      = request.form.get("mut_dna", "").strip()
    protein_name = request.form.get("protein_name", "").strip()

    # --- Step 0: Validate ---
    errors = validate_inputs(ref_dna, mut_dna, protein_name)
    if errors:
        return render_template("index.html", errors=errors,
                               ref_dna=ref_dna, mut_dna=mut_dna,
                               protein_name=protein_name)

    # --- Step 1 & 2: Align sequences and detect mutations ---
    aligner = SequenceAligner(ref_dna, mut_dna)
    try:
        aligner.align()
        mutations = aligner.detect_mutations()
    except Exception as e:
        return render_template("index.html",
                               errors=[f"Alignment error: {str(e)}"],
                               ref_dna=ref_dna, mut_dna=mut_dna,
                               protein_name=protein_name)

    if not mutations:
        return render_template("index.html",
                               errors=["No mutations detected between the two sequences."],
                               ref_dna=ref_dna, mut_dna=mut_dna,
                               protein_name=protein_name)

    # --- Step 3: Translate DNA → Protein and compare AA sequences ---
    translator = ProteinTranslator()
    try:
        aa_changes = translator.compare_proteins(ref_dna, mut_dna)
    except Exception:
        aa_changes = []

    # Attach the first matching AA change label to each mutation object
    for i, mut in enumerate(mutations):
        if i < len(aa_changes):
            mut.aa_change = aa_changes[i]["change"]

    # --- Step 4: Fetch 3D structure from RCSB PDB ---
    pdb_id        = None
    pdb_structure = None
    pdb_error     = None
    try:
        pdb = PDBFetcher(protein_name)
        pdb_id = pdb.search_pdb()
        if pdb_id:
            pdb_structure = pdb.fetch_structure(pdb_id)
    except ConnectionError as e:
        pdb_error = str(e)

    # --- Step 5: Fetch protein annotation from UniProt ---
    uniprot_info  = {}
    uniprot_error = None
    try:
        uni = UniProtFetcher(protein_name)
        uniprot_info = uni.get_info()
    except ConnectionError as e:
        uniprot_error = str(e)

    # --- Build Plotly chart data ---
    seq_length = len(ref_dna)
    chart_data = {
        "positions":  [m.position for m in mutations],
        "types":      [m.type     for m in mutations],
        "labels":     [f"{m.ref_base}→{m.mut_base}" for m in mutations],
        "seq_length": seq_length
    }

    return render_template(
        "results.html",
        mutations=mutations,
        aa_changes=aa_changes,
        pdb_id=pdb_id,
        pdb_structure=pdb_structure,
        pdb_error=pdb_error,
        uniprot_info=uniprot_info,
        uniprot_error=uniprot_error,
        protein_name=protein_name,
        ref_protein=translator.ref_protein,
        mut_protein=translator.mut_protein,
        chart_data=json.dumps(chart_data),
        seq_length=seq_length
    )


@app.route("/reset")
def reset():
    """Clears the form — navigates back to a fresh input page."""
    return render_template("index.html")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
   import webbrowser
webbrowser.open("http://127.0.0.1:5000")
app.run(debug=True, host="0.0.0.0", port=5000)
