# Standard library imports:
from pathlib import Path

# External libraries:
import unittest
import pandas as pd

# Self-made libraries:
from biolib.classes.files.pdbcif import PdbCifFile, PdbCifFileCollection

class TestPdbCifFile(unittest.TestCase):

    def setUp(self):
        self.cif_simple: PdbCifFile = PdbCifFile(Path('cif_files/simple.cif'))
        self.cif_with_quotes: PdbCifFile = PdbCifFile(Path('cif_files/with_quotes.cif'))
        self.cif_multiple_polymer_entities = PdbCifFile(Path('cif_files/multiple_polymer_entities.cif'))
        
    def test_notPathObject(self):
        with self.assertRaises(TypeError):
            PdbCifFile('not/a/path/object')
            PdbCifFile(2)

    def test_nonExistingFile(self):
        with self.assertRaises(FileNotFoundError):
            PdbCifFile(Path('cif_files/non_existing.cif'))

    def test_notAFile(self):
        with self.assertRaises(FileNotFoundError):
            PdbCifFile(Path('cif_files'))
        
    def test_wrongSuffix(self):
        with self.assertRaises(ValueError):
            PdbCifFile(Path('cif_files/wrong_suffix.ciff'))

    def test_countDataBlocks(self):
        self.assertEqual(4, self.cif_simple.countDataBlocks())

    def test_countLoopBlocks(self):
        self.assertEqual(1, self.cif_simple.countLoopBlocks())

    def test_categoryExists(self):
        self.assertTrue(self.cif_simple.categoryExists('_entity_poly'))  # non-loop block
        self.assertFalse(self.cif_simple.categoryExists('_atom_site'))  # existing loop block, should return false
        self.assertFalse(self.cif_simple.categoryExists('_symmetry'))  # not present

    def test_loopCategoryExists(self):
        self.assertTrue(self.cif_simple.loopCategoryExists('_atom_site'))
        self.assertFalse(self.cif_simple.loopCategoryExists('_entity_poly')) # non-loop block
        self.assertFalse(self.cif_simple.loopCategoryExists('_struct_sheet')) # not present

    def test_categoryToDf(self):
        cif_simple_atom_site: pd.DataFrame = pd.DataFrame(
            data=[['ATOM', '1', 'ALA', 'CA', '0.5', '0.5', '0.5'],['HETATM', '2', 'LIG', 'C1', '0.0', '1.0', '2.0']],
            columns = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.label_comp_id', '_atom_site.label_atom_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

        #print(cif_simple_atom_site)
        #print(self.cif_simple.loopCategoryToDf('_atom_site'))
        
        pd.testing.assert_frame_equal(cif_simple_atom_site, self.cif_simple.categoryToDf('_atom_site'))

    def test_getAminoAcidSequences(self):
        self.assertEqual(['ACDE'], self.cif_simple.getAminoAcidSequences())
        self.assertEqual(['ACDE', 'AGYFKR'], self.cif_multiple_polymer_entities.getAminoAcidSequences())

    def test_countEntities(self):
        self.assertEqual(1, self.cif_simple.countEntities())
        self.assertEqual(2, self.cif_multiple_polymer_entities.countEntities())
        
    def test_countPolymerEntities(self):
        self.assertEqual(1, self.cif_simple.countPolymerEntities())
        self.assertEqual(2, self.cif_multiple_polymer_entities.countPolymerEntities())

        
class TestPdbCifFileCollection(unittest.TestCase):

    def setUp(self):
        self.cif_collection = PdbCifFileCollection(Path('cif_files'))

    def test_writeSequencesToFasta(self):
        pass

if __name__ == "__main__":
    unittest.main()
