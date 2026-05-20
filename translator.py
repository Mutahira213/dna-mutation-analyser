from Bio.Seq import Seq


class ProteinTranslator:
    """
    Translates DNA sequences to protein using BioPython's Seq.translate(),
    then compares the resulting protein sequences to identify amino acid changes.
    Pipeline: DNA → RNA (implicit) → Protein via genetic code.
    """

    def __init__(self):
        self.ref_protein = None
        self.mut_protein = None

    def translate(self, dna_seq: str) -> str:
        """
        Translates a DNA string to an amino acid sequence.
        - Strips gaps and whitespace
        - Trims to nearest multiple of 3 (avoids partial codon errors)
        - Stops at first stop codon (to_stop=True)
        """
        dna = Seq(dna_seq.upper().strip().replace('-', ''))
        trimmed = dna[:len(dna) - len(dna) % 3]
        protein = trimmed.translate(to_stop=True)
        return str(protein)

    def compare_proteins(self, ref_dna: str, mut_dna: str) -> list:
        """
        Translates both DNA sequences and returns a list of amino acid changes.
        Each change is a dict with: position, ref_aa, mut_aa, change (formatted string).
        Also flags length differences caused by frameshifts.
        """
        self.ref_protein = self.translate(ref_dna)
        self.mut_protein = self.translate(mut_dna)

        changes = []
        min_len = min(len(self.ref_protein), len(self.mut_protein))

        for i in range(min_len):
            if self.ref_protein[i] != self.mut_protein[i]:
                changes.append({
                    "position": i + 1,
                    "ref_aa": self.ref_protein[i],
                    "mut_aa": self.mut_protein[i],
                    "change": f"{self.ref_protein[i]} → {self.mut_protein[i]}"
                })

        # Handle length difference from frameshift / premature stop codon
        if len(self.ref_protein) != len(self.mut_protein):
            changes.append({
                "position": min_len + 1,
                "ref_aa": self.ref_protein[min_len] if min_len < len(self.ref_protein) else "-",
                "mut_aa": self.mut_protein[min_len] if min_len < len(self.mut_protein) else "-",
                "change": "Length difference (possible frameshift)"
            })

        return changes
