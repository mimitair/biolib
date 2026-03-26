import unittest
from pathlib import Path

from biolib.classes.files.pdbcif import PdbCifFile, PdbCifFileCollection

class TestPdbCifFile(unittest.TestCase):

    def setUp(self):
        cif_files_path: Path = Path('cif_files')

    def test_getAmountOfDataBlocks(self):
        pass
    
    def test_init(self):
        print(self.mycif.file_path)

    def test_getAminoAcidSequences(self):
        print(self.mycif.getAminoAcidSequences())

    
if __name__ == "__main__":
    unittest.main()
