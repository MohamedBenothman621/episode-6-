# -----------------------------------------------------------------------------
#  Ep 6 - Hamming Distance & Approximate Pattern Matching
# Find pattern positions allowing up to 'd' mismatches in a DNA sequence
# -----------------------------------------------------------------------------

def HammingDistance(seq1, seq2):
    """
    Calculate the Hamming Distance between two DNA sequences of equal length.

    The Hamming Distance counts the number of positions where the two
    sequences differ — useful for measuring similarity between DNA strands.

    Args:
        seq1 (str): The first DNA sequence.
        seq2 (str): The second DNA sequence.

    Returns:
        int: The number of differing positions between the two sequences.

    Example:
        >>> HammingDistance("ATGC", "AAGC")
        1
    """
    counter = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            counter += 1
    return counter


def ApproximatePatternMatching(Text, Pattern, d):
    """
    Find all positions where a Pattern matches the Text with at most 'd' mismatches.

    Unlike exact pattern matching, this function allows up to 'd' mismatches
    between the pattern and the substring — making it robust to mutations
    or sequencing errors in real genomic data.

    Args:
        Text    (str): The parent DNA sequence to search within.
        Pattern (str): The DNA motif/pattern to search for.
        d       (int): Maximum number of allowed mismatches (Hamming Distance).

    Returns:
        list: A list of starting indices where Pattern approximately matches Text.

    Example:
        >>> ApproximatePatternMatching("ATGCATATG", "ATAG", 1)
        [4]
    """
    positions = []
    n = len(Text)
    m = len(Pattern)

    for i in range(n - m + 1):
        if HammingDistance(Text[i:i+m], Pattern) <= d:
            positions.append(i)

    return positions


# Example usage - Ep 6
Text    = "ATGCATATGACTACTAGATACTGATACTGATACATA"
Pattern = "ATAG"
d       = 1

print("=" * 50)
print(" Ep 6 - Approximate Pattern Matching")
print(f"   Sequence      : {Text}")
print(f"   Pattern       : {Pattern}")
print(f"   Max Mismatches: {d}")
print(f"   Positions     : {ApproximatePatternMatching(Text, Pattern, d)}")
print()
print("=" * 50)
print("✅ Episode 6 executed successfully.")
print("=" * 50)
