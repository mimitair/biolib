# Standard library imports:
from pathlib import Path

# External libraries:
import unittest
import pandas as pd
import numpy as np

# Self-made libraries:
from biolib.files.pdbcif import PdbCifFile, PdbCifFileCollection

class TestPdbCifFile(unittest.TestCase):

    def setUp(self):
        self.cif_simple: PdbCifFile = PdbCifFile(Path('cif_files/simple.cif'))
        self.cif_with_quotes: PdbCifFile = PdbCifFile(Path('cif_files/with_quotes.cif'))
        self.cif_multiple_polymer_entities = PdbCifFile(Path('cif_files/multiple_polymer_entities.cif'))
        self.cif_multiline = PdbCifFile(Path('cif_files/multiline.cif'))
        
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

    def test_categoryToDf_simple(self):
        # Define what the dataframe should look like:
        cif_simple_atom_site: pd.DataFrame = pd.DataFrame(
            data=[['ATOM', '1', 'ALA', 'CA', '0.5', '0.5', '0.5'],['HETATM', '2', 'LIG', 'C1', '0.0', '1.0', '2.0']],
            columns = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.label_comp_id', '_atom_site.label_atom_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

        # Do the test:
        pd.testing.assert_frame_equal(cif_simple_atom_site, self.cif_simple.categoryToDf('_atom_site'))

    def test_categoryToDf_with_quotes1(self):
        # Define what the dataframe should look like:
        cif_with_quotes_atom_site: pd.DataFrame = pd.DataFrame(
            data=[['ATOM', 'GLY', 'CA', '0.0', '0.0', '0.0'],
                  ['ATOM', 'LYS', 'NZ', '1.1', '1.2', '1.3'],
                  ['HETATM', 'LIG', "C1'", '1.0', '2.0', '3.0']],
            columns=['_atom_site.group_PDB', '_atom_site.label_comp_id', '_atom_site.label_atom_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']
        )
        
        # Do the test:
        pd.testing.assert_frame_equal(cif_with_quotes_atom_site, self.cif_with_quotes.categoryToDf('_atom_site'))
        
    def test_categoryToDf_with_quotes2(self):
        # Define what the dataframe should look like:
        cif_with_quotes_chem_comp: pd.DataFrame = pd.DataFrame(
            data=[["'C3 H7 N O2'", '69.420', 'ALA', "'l-peptide linking'", '?'],
                  ["'C6 H15 N4 O2 1'", '175.2', 'ARG', "'l-peptide linking'", '?']],
            columns=['_chem_comp.formula', '_chem_comp.formula_weight', '_chem_comp.id', '_chem_comp.name', '_chem_comp.pdbx_synonyms']
        )

        # Do the test:
        pd.testing.assert_frame_equal(cif_with_quotes_chem_comp, self.cif_with_quotes.categoryToDf('_chem_comp'))
    
    def test_categoryToDf_multiline(self):
        # Define what the dictionary should look like:
        cif_multiline_entity_poly: dict = {'_entity_poly.entity_id': '1', '_entity_poly.pdbx_seq_one_letter_code': 'ACDEFGHIKLMN'}

        # Do the test
        self.assertEqual(cif_multiline_entity_poly, self.cif_multiline.categoryToDf('_entity_poly'))


    def test_getHeteroAtoms(self):
        pass
    
    def test_getAminoAcidSequences_simple(self):
        self.assertEqual(['ACDE'], self.cif_simple.getAminoAcidSequences())

    def test_getAminoAcidSequences_multiple_polymer_entities(self):
        self.assertEqual(['ACDE', 'AGYFKR'], self.cif_multiple_polymer_entities.getAminoAcidSequences())

    def test_countEntities(self):
        self.assertEqual(1, self.cif_simple.countEntities())
        self.assertEqual(3, self.cif_multiple_polymer_entities.countEntities())
        
    def test_countPolymerEntities(self):
        self.assertEqual(1, self.cif_simple.countPolymerEntities())
        self.assertEqual(3, self.cif_multiple_polymer_entities.countPolymerEntities())

    def test_countPolypeptideEntities(self):
        self.assertEqual(1, self.cif_simple.countPolypeptideEntities())
        self.assertEqual(2, self.cif_multiple_polymer_entities.countPolypeptideEntities())

    def test_atomSiteToDf_simple(self):
        # Define what the dataframe should look like (with correct typing):
        cif_simple_atom_site: pd.DataFrame = pd.DataFrame(
            data=[['ATOM', np.uint32(1), 'ALA', 'CA', np.float16(0.5), np.float16(0.5), np.float16(0.5)],['HETATM', np.uint32(2), 'LIG', 'C1', np.float16(0.0), np.float16(1.0), np.float16(2.0)]],
            columns = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.label_comp_id', '_atom_site.label_atom_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

        # Do the test:
        pd.testing.assert_frame_equal(cif_simple_atom_site, self.cif_simple.atomSiteToDf(), check_dtype=True, check_column_type=True)
        
    def test_filterAtomSite_simple_on_alanine_residue(self):
        # Define what the result should look like
        cif_simple_atom_site_filtered: pd.DataFrame = pd.DataFrame(
            data=[['ATOM', np.uint32(1), 'ALA', 'CA', np.float16(0.5), np.float16(0.5), np.float16(0.5)]],
            columns = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.label_comp_id', '_atom_site.label_atom_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

        # Do the test:
        pd.testing.assert_frame_equal(cif_simple_atom_site_filtered, self.cif_simple.filterAtomSite(residue_names=['ALA'], atom_names=['CA']), check_dtype=True, check_column_type=True)
    

    
class TestPdbCifFileCollection(unittest.TestCase):

    def setUp(self):
        self.cif_collection = PdbCifFileCollection(Path('cif_collection'))

    def test_writeSequencesToFasta(self):
        pass

if __name__ == "__main__":
    unittest.main()
