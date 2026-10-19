"""Rate at which the L1 endonuclease motif is found by chance with the criterion used in annotate_known_and_motif.py:
any 6-mer in a 9 bp window within one mismatch of 5'-TTAAAA-3' or its reverse complement (5'-TTTTAA-3').
Random 9 bp sequences are drawn with the base composition of the human genome (59% A+T)."""
import random

def motif_ok(seq):
    for target in ("TTAAAA", "TTTTAA"):
        for i in range(len(seq) - 5):
            if sum(a != b for a, b in zip(seq[i:i + 6], target)) <= 1:
                return True
    return False

rng = random.Random(1)
weights = [0.295, 0.205, 0.205, 0.295]  # A, C, G, T
n = 200_000
hits = sum(motif_ok("".join(rng.choices("ACGT", weights, k=9))) for _ in range(n))
print(f"motif present by chance in {hits / n:.1%} of {n:,} random 9 bp windows")
