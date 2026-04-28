# Custom imports:
from biolib.classes.files.fasta import FastaFile
from biolib.util import regex_patterns

# Standard library:
from pathlib import Path
import sys

"""
This script takes a fasta file as the first positional argument and returns
which sequences match the canonical GxSxG patterns representing the catalytic serine
in esterases.
"""

def main():
    path_to_fasta: Path = Path(sys.argv[1])

    myfasta: FastaFile = FastaFile(path_to_fasta)
    
    print(f'Detected {myfasta.count} sequences.')
    # Look for sequences matching the catalytic serine pattern in this fasta file:
    matches: list = myfasta.matchPattern(regex_patterns.aa_esterase_catalytic_serine)

    unique_matches: set = {match[0] for match in matches}
    print(f'Found {len(matches)} matches across {len(unique_matches)} sequences.')
    print(*matches, sep='\n')

if __name__ == "__main__":
    main()
