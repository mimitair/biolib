#! /usr/bin/env python3

"""
Script to find three neigbouring residues based on distance thresholds in a cif file.
"""

from pathlib import Path
from biolib.files.pdbcif import PdbCifFile
import sys

RES1_NAME: str = "SER"
RES2_NAME: str = "HIS"
RES3_NAME: str = "ASP"
ATOM1_NAME: str = "CA"
ATOM2_NAME: str = "CA"
ATOM3_NAME: str = "CA"
MAX_DIST1_2: float = 8.5
MAX_DIST2_3: float = 5.1

def main():
    # Read the arguments:
    in_file: Path = Path(sys.argv[1])
#    out_file: Path = Path(sys.argv[2])

    # Convert input file to Cif object:
    cif: PdbCifFile = PdbCifFile(in_file)
    
    # Find the defined triad of residues:
    triads: list = cif.findTriads(RES1_NAME,
                                  RES2_NAME,
                                  RES3_NAME,
                                  MAX_DIST1_2,
                                  MAX_DIST2_3)

    print("Found the following triads:")
    print([(atom_array.res_id,atom_array.chain_id) for atom_array in triads])
    return 0
    
if __name__ == "__main__":
    main()
