#!/usr/bin/env python3

from pathlib import Path
from biolib.files.fasta import FastaFile
import sys

def main():

    in_file: Path = Path(sys.argv[1])
    out_dir: Path = Path(sys.argv[2])

    # Convert in_file to FastaFile object:
    myfasta: FastaFile = FastaFile(in_file)

    # Split in multiple files:
    myfasta.splitToSeparateFiles(out_dir)
    
    return 0

if __name__ == "__main__":
    main()
