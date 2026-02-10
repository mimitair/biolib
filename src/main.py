from files.pdbcif import PDBCIFFile
from pathlib import Path

def main():
    my_file = PDBCIFFile(Path("/home/green/projects/bio/data/260131_pdb_ec3.1.1.74/1XZL_full.cif"))
    
    print(my_file.getHeteroAtoms())
    
if __name__ == "__main__":
    main()
