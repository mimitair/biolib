from files.pdbx_mmcif_file import PDBCIFFile
from pathlib import Path

def main():
    my_file = PDBCIFFile(Path("5XJH.cif"))
    
    df_audit = my_file.loopCategoryToDf("site")

    print(df_audit)


if __name__ == "__main__":
    main()
