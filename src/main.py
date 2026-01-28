from files.pdbx_mmcif import PDBCIFFile
from pathlib import Path

def main():
    my_file = PDBCIFFile(Path("5XJH.cif"))
    
    df_atoms = my_file.loopCategoryToDf("atom_site")

    ligands = my_file.getHeteroAtoms(True)

    print(ligands)
    
if __name__ == "__main__":
    main()
