from pathlib import Path
import sys
from biolib.files.pdbcif import PdbCifFileCollection

def main():
    in_dir: Path = Path(sys.argv[1])
    out_file: Path = Path(sys.argv[2])

    # Convert in_dir to collection object:
    mycoll: PdbCifFileCollection = PdbCifFileCollection(in_dir)
    return 0

if __name__ == "__main__":
    main()
