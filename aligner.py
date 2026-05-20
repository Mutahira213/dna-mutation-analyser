from Bio.Align import PairwiseAligner
from dataclasses import dataclass


@dataclass
class MutationResult:
    """Data container for a single detected mutation."""
    position: int
    type: str        # substitution | insertion | deletion
    ref_base: str
    mut_base: str
    aa_change: str = ""


class SequenceAligner:
    """
    Aligns two DNA sequences using the Needleman-Wunsch global alignment
    algorithm (via BioPython PairwiseAligner) and detects mutations.
    """

    def __init__(self, ref_seq: str, mut_seq: str):
        self.ref_seq = ref_seq.upper().strip()
        self.mut_seq = mut_seq.upper().strip()
        self._aligned_ref = None
        self._aligned_mut = None

    def align(self):
        """
        Runs global pairwise alignment and stores the best aligned strings.
        Returns: (aligned_ref, aligned_mut) as strings with gap characters.
        """
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        alignments = list(aligner.align(self.ref_seq, self.mut_seq))
        if not alignments:
            raise ValueError("Alignment failed — sequences may be incompatible.")
        best = alignments[0]

        # Reconstruct aligned strings from coordinate blocks
        ref_aligned = []
        mut_aligned = []
        coords = best.coordinates
        for i in range(coords.shape[1] - 1):
            ref_start, ref_end = coords[0, i], coords[0, i + 1]
            mut_start, mut_end = coords[1, i], coords[1, i + 1]
            ref_len = ref_end - ref_start
            mut_len = mut_end - mut_start
            max_len = max(ref_len, mut_len)
            ref_seg = self.ref_seq[ref_start:ref_end] + '-' * (max_len - ref_len)
            mut_seg = self.mut_seq[mut_start:mut_end] + '-' * (max_len - mut_len)
            ref_aligned.append(ref_seg)
            mut_aligned.append(mut_seg)

        self._aligned_ref = ''.join(ref_aligned)
        self._aligned_mut = ''.join(mut_aligned)
        return self._aligned_ref, self._aligned_mut

    def detect_mutations(self) -> list:
        """
        Compares aligned sequences position by position to detect:
        - substitutions (base change)
        - insertions (gap in reference)
        - deletions (gap in mutant)
        Returns a list of MutationResult objects.
        """
        if self._aligned_ref is None:
            self.align()

        mutations = []
        ref_pos = 0

        for r, m in zip(self._aligned_ref, self._aligned_mut):
            if r == m:
                if r != '-':
                    ref_pos += 1
                continue

            if r == '-':
                # Gap in reference = insertion in mutant
                mutations.append(MutationResult(
                    position=ref_pos,
                    type="insertion",
                    ref_base="-",
                    mut_base=m
                ))
            elif m == '-':
                # Gap in mutant = deletion from reference
                ref_pos += 1
                mutations.append(MutationResult(
                    position=ref_pos,
                    type="deletion",
                    ref_base=r,
                    mut_base="-"
                ))
            else:
                # Different bases = substitution
                ref_pos += 1
                mutations.append(MutationResult(
                    position=ref_pos,
                    type="substitution",
                    ref_base=r,
                    mut_base=m
                ))

        return mutations
