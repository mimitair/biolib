from biolib.classes.files.pdbcif import PDBCIFFile, PDBCIFFileCollection
from biolib.classes.files.fasta import FastaFile
from pathlib import Path

def test

def main():
    #my_file = PDBCIFFile(Path("test_data/3jti_full.cif"))

    #print(my_file.matchCategory('entity_poly'))
    #print(my_file.getHeteroAtoms())
    #print(my_file.categoryToDict('entity_poly')['pdbx_seq_one_letter_code'])
    #my_collection = PDBCIFFileCollection(Path('test_data'))
    #my_collection.writeSequencesToFasta()

    #print(my_file.matchLoopCategory('entity_poly'))
   
    my_fasta = FastaFile(Path('test_data/esterases.fasta'))

    print(my_fasta.getAmountOfEntries())
if __name__ == "__main__":
    main()
