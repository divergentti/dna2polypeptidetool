"""
Utility functions for the DNA to Polypeptide Encoder.
"""

def reverse_complement(seq):
    """Generate the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))