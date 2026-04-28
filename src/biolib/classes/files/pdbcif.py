# Default library imports:
import sys
import re
from pathlib import Path
import logging
import itertools

# External libraries:
import pandas as pd
from parsnip import CifFile
import numpy as np

# Logger config:
logger = logging.getLogger(__name__)

class PdbCifFile:
    """
    This class represents a PDBx/mmCIF file as described in https://mmcif.wwpdb.org/, and contains methods to parse and manipulate them.
    """
    def __init__(self, path_to_pdbcif: Path) -> None:
        """
        Initializes the PdbCifFFile object.
        
        Input:
            - path_to_pdbcif: Path: Path object to the cif file.

        Returns:
            - None
        """
        logger.info("--- Initializing PdbCifFile object ---")
        
        ### DEFENSIVE CHECKS ###
        if not isinstance(path_to_pdbcif, Path):
            raise TypeError("Path to PDB/mmCIF file must be a Path object. Convert it to a Path object before initiating the class.")
        
        if not path_to_pdbcif.exists():
            raise FileNotFoundError(f"{path_to_pdbcif} does not exist. Check if you have provided the correct file path.")
            
        if not path_to_pdbcif.is_file():
            raise FileNotFoundError(f"{path_to_pdbcif} is not a file. Check if you have provided a file path instead of a directory.")

        if not path_to_pdbcif.suffix == '.cif':
            raise ValueError(f"{path_to_pdbcif} is not a .cif file. This class only accepts .cif files.")
        
        logger.info("Defensive checks OK.")

        ### INIT ###
        self.full_path: Path = path_to_pdbcif.resolve()  # Resolved path to the CIF file
        self.name: str = path_to_pdbcif.stem  # Name of the file without suffix and prefix
        #parsnip: CifFile = CifFile(path_to_pdbcif)  # The given cif file as a parsnip CifFile object

        logger.info(f"Full path to CIF file: {self.full_path}")        

    def countDataBlocks(self) -> int:
        """
        Returns the total amount of data blocks in this cif file.
        """
        return self.full_path.read_text().count('#') -1

    def countLoopBlocks(self) -> int:
        """
        Returns the amount of loop data blocks in this cif file.
        """
        return self.full_path.read_text().count('loop_')
        
    def categoryExists(self, category: str) -> bool:
        """
        Returns true if the given category exists in the cif file based on regex matches. Works only for non-loop data blocks.

        Input:
            - category: str: the desired category to extract.

        Returns:
            - bool: True if the given category exists, False if not.
        """
        # Extract file contents:
        file_contents: str = self.full_path.read_text()  # Content of the file as a string.
        
        # Define the pattern:
        pattern = rf"#\s*\n{category}[^#]+#"

        # Search in the file contents:
        m = re.search(pattern, file_contents)

        return True if m is not None else False

    def loopCategoryExists(self, category: str) -> bool:
        """
        Returns True if the given loop data block exists in the cif file based on regex matches.

        Input:
            - category: str: The loop category to parse

        Returns
            - bool: True if the given category exists, False if not.
        """
        file_contents: str = self.full_path.read_text()  # Content of the file as a string.
        
        pattern = rf"#\s*\nloop_\s*\n{category}[^#]+#"
        
        # Search for the pattern:
        m = re.search(pattern, file_contents)

        return True if m is not None else False

    def categoryToDf(self, category: str) -> pd.DataFrame | dict: 
        """
        Convert any data block, given its category name, to a pandas dataframe or series.
        This function uses the parnsip external library.

        Input:
            - category: str: The category to extract (f.e.: "atom_site")

        Returns:
            - pd.DataFrame: A loop data block as a dataframe. Beware that no further processing is done. Every value is essentially a string.
            - dict: A non-loop data block as a dict.

        Raises:
            - ValueError: If the given category does not exist.
        """
        if self.categoryExists(category):
            return {key:value for key,value in self.toParsnip().pairs.items() if category + '.' in key}
        
        elif self.loopCategoryExists(category):
            return pd.DataFrame([arr for arr in self.toParsnip().loops if category + '.' in ''.join(arr.dtype.names)][0].reshape(-1))

        else:
            raise ValueError(f'{category} does not exist in {self.name}')
        
    def getHeteroAtoms(self) -> set:
        """
        Returns a set of hetero atoms (labeled as 'HETATM') in the _atom_site category.

        Returns:
            - set: A set of HETATM names in this PDBx/mmCIF file. Empty set if nothing is found

        """
        # Convert the atom_site category to a dataframe:
        df_atoms = self.loopCategoryToDf('atom_site')

        # Return the 'label_comp_id' column for each row where 'group_PDB' == 'HETATM' as a set
        return set(df_atoms[df_atoms['group_PDB'] == 'HETATM']['label_comp_id'].to_list())

    def toParsnip(self):
        return CifFile(self.full_path)
    
    def getAminoAcidSequences(self) -> list:
        """
        Returns the amino acid sequences of each polymer entity in this file as a list of strings.
        If the cif file has only one polymer entity, a list of length 1 is returned.

        Returns:
            - list: A list of the amino acid sequences for each polymer entity in this PDBx/mmCIF file.

        Raises:
            - ValueError: If there are no polymer entities.
        """
        # Extract the data block where we can find the amino acid sequences:
        df: pd.DataFrame | pd.Series = self.categoryToDf('_entity_poly')

        # TODO checking type every time is too slow?
        if isinstance(df, pd.DataFrame):
            # Extract only the 
            df = df[(df['_entity_poly.type'] == 'polypeptide(L)') | (df['_entity_poly.type'] == 'polypetide(D)')]
            return [seq.replace('\n', '').replace(';', '') for seq in df['_entity_poly.pdbx_seq_one_letter_code'].tolist()]

        elif isinstance(df, dict):
            if (df['_entity_poly.type'] == 'polypeptide(L)') or (df['_entity_poly.type'] == 'polypeptide(D)'):
                return [df['_entity_poly.pdbx_seq_one_letter_code'].replace('\n', '').replace(';', '')]
            else:
                return []

    def countEntities(self) -> int:
        """
        Returns the amount of entities in this cif file.
        """
        return len(self.toParsnip()['_entity.id'])

    def countPolymerEntities(self) -> int:
        """
        Counts the amount of entities labeled as 'polymer' in this cif file.
        """
        return len([entity for entity in np.nditer(self.toParsnip()['_entity.type']) if entity == 'polymer'])

    def countPolypeptideEntities(self) -> int:
        pass
        

        
    ##################################################
    ##### EVERYTHING UNDERNEATH IS NOT FUNCTIONAL ####
    ##################################################
    
    def residueNumberToResidueName(self, residue_number: int) -> str:
        """
        Input:
            - residue_number: Number of an amino acid residue in the PDB/mmCIF file.
            
        Returns:
            - str: The name of the residue as a string.
        """
        return self.df_atoms.loc[self.df_atoms['label_seq_id'] == str(residue_number)]['label_comp_id'].iloc[0]
    
    
    def residueAtomNamesToPoints(self, residue_name: str, atom_name: str) -> dict:
        """
        Input:
            - residue_name: Name of an amino acid residue (fe 'GLY').
            - atom_name: Name of the atom as it is represented in the PDB/mmCIF file (fe 'CA').
        
        Returns:
            - dict: A dictionary with residue numbers as keys and the atom, represented as a Point object, as values.
        """
        # Empty dict to store results:
        result: dict = {}
        
        # Extract the rows where atom and residue equal that of what is given:
        df = self.df_atoms.loc[(self.df_atoms['label_comp_id'] == residue_name) & (self.df_atoms['label_atom_id'] == atom_name)]

        # Convert the xyz coordinates of every row to a Point object, and store in a dictionary with the residue number as key:
        # Iterate over rows as named tuples:
        # We can access the column name through 'row.<column_name>'
        # row.label_seq_id is the residue number in string format, so we typecast to int.
        for row in df.itertuples():
            result[int(row.label_seq_id)] = self.atomNumberToPoint(int(row.Index))        
        
        return result
    
    
    def atomNumberToSeries(self, atom_number: int) -> pd.Series:
        """
        Input:
            - atom_number: The number of the atom in the PDB/mmCIF file.
            
        Returns:
            - pd.Series: The atom as a pandas Series.
        """
        return self.df_atoms.loc[[str(atom_number)]]
    
    
    def atomNumberToPoint(self, atom_number: int) -> 'Point':
        """
        Input:
            - atom_number: The number of the atom in the PDB/mmCIF file ('id' column).
            
        Returns:
            - Point: The atom as a Point object with x,y,z coordinates.
        """
        atom = self.df_atoms.loc[[str(atom_number)]]
        x = float(atom['Cartn_x'].iloc[0])
        y = float(atom['Cartn_y'].iloc[0])
        z = float(atom['Cartn_z'].iloc[0])
        
        return point.Point(x,y,z)
    
    
    def atomNumberToResidueNumber(self, atom_number: int) -> int:
        """
        Input: 
            - atom_number: Number of the atom in the PDB/mmCIF file ('id' column).
        
        Returns:
            - int: The residue number of which this atom is part.
        
        """
        return int(self.df_atoms.loc[[str(atom_number)]]['label_seq_id'].iloc[0])
        
    
    def findTriads(self, 
                   res1: str, 
                   res2: str, 
                   res3: str, 
                   atom1: str = 'CA',
                   atom2: str = 'CA',
                   atom3: str = 'CA',
                   min_dist_1_2: float = 1, 
                   min_dist_2_3: float = 1,
                   min_dist_3_1: float = 1,
                   max_dist_1_2: float = 4,
                   max_dist_2_3: float = 4,
                   max_dist_3_1: float = 4):
        """
        Input:
            - res1: Name of the first residue in the triad.
            - res2: Name of the second residue in the triad.
            - res3: Name of the third residue in the triad.
            - atom1: Name of the atom to consider in res1. [Default = 'CA', which is the alpha carbon]
            - atom2: Name of the atom to consider in res2. [Default = 'CA', which is the alpha carbon]
            - atom3: Name of the atom to consider in res3. [Default = 'CA', which is the alpha carbon]
            - min_dist_1_2: Minimum distance between atom1 and atom2, in Angstroms. [Default = 1]
            - min_dist_2_3: Minimum distance between atom2 and atom3, in Angstroms. [Default = 1]
            - min_dist_3_1: Minimum distance between atom3 and atom1, in Angstroms. [Default = 1]
            - max_dist_1_2: Maximum distance between atom1 and atom2, in Angstroms. [Default = 4]
            - max_dist_2_3: Maximum distance between atom2 and atom3, in Angstroms. [Default = 4]
            - max_dist_3_1: Maximum distance between atom3 and atom1, in Angstroms. [Default = 4]
            
        Returns:
            ?
        """
        # Import the Point class and printTriad function::
        from .util import Point, printTriad

        # Print a pretty ASCII art triad:
        print("Looking for the following triad(s):")
        printTriad(res1, res2, res3, atom1, atom2, atom3, min_dist_1_2, min_dist_2_3, min_dist_3_1, max_dist_1_2, max_dist_2_3, max_dist_3_1)
        
        # Define the atom loop block as a dataframe
        df_atoms = self.loopCategoryToDf("atom_site")

        # Initiate empty lists to store distances between pairs of Points (atoms):
        distances_1_2: list = []
        distances_2_3: list = []
        distances_3_1: list = []
        
        # Extract only the rows where residue name and atom name match what is given.
        df_atoms = self.df_atoms.loc[(self.df_atoms['label_comp_id'] == residue_name) & (self.df_atoms['label_atom_id'] == atom_name)]

        # Store each given atom in a {residue_number : Point} dictionary
         ## do this directly with pandas, needs further implementation, use dict comprehension?
        res1_points: dict
        res1_points: dict = self.residueAtomNamesToPoints(res1, atom1)
        res2_points: dict = self.residueAtomNamesToPoints(res2, atom2)
        res3_points: dict = self.residueAtomNamesToPoints(res3, atom3)
        
        # Now we have to calculate the distances between each pair of Point objects.
        # This requires some nested for loops
        # Iterate over the given atoms of residue 1:
        for key1, point1 in res1_points.items():  # key1=residue number, point1=Point object representing an atom of the aa residue
            # Iterate over the given atoms of residue 2:
            for key2, point2 in res2_points.items():
                distance_1_2 = point1.distance(point2)  # Distance between p1 and p2 (see Point class in util)
                # If it falls within the given thresholds, append to the list:
                if (distance_1_2 >= min_dist_1_2) & (distance_1_2 <= max_dist_1_2):
                    distances_1_2.append((key1, key2, distance_1_2))  # Append to list
            # Now the same for 3 and 1, since we are iterating over 1 anyway:
            for key3, point3 in res3_points.items():
                distance_3_1 = point1.distance(point3)  # d(p1-p3)
                # Append to lists:
                if (distance_3_1 >= min_dist_3_1) & (distance_3_1 <= max_dist_3_1):
                    distances_3_1.append((key3, key1, distance_3_1))
            
        # Lastly for the 2-3 pairs:
        for key2, point2 in res2_points.items():
            for key3, point3 in res3_points.items():
                distance_2_3 = point2.distance(point3)  # d(p2-p3)
                # Append to lists:
                if (distance_2_3 >= min_dist_2_3) & (distance_2_3 <= max_dist_2_3):
                    distances_2_3.append((key2, key3, distance_2_3))
        
        # Now we have, for each possible pair in the triad, a list of tuples with each tuple containing
        # (res_numberX, res_numberY, d(atomX-atomY))
        # Only the distances within the given thresholds are stored.
        
        # The next step is to extract only the residue numbers that can actually form a triangle of points.
        # i.e. a RES1-RES2 pair in which RES1 does not occur in RES1-RES3 pairs, 
        # or RES2 does not occur in RES2-RES3 pair, needs to be filtered out
        # Fot it to be a valid triad, the residue number of a given RES MUST occur in its two possible pairings
        
        # Define res1 in res1-res2 pairs and in res3-res1 pairs
        res1a: set = {res1 for res1,res2,_ in distances_1_2}    
        res1b: set = {res1 for res3,res1,_ in distances_3_1}
        # Take the intersection of both sets (i.e. filter out the ones that do not occur in both pairs):
        res1_set: set = res1a.intersection(res1b)
    
        # Repeat for res2 and res3:
        res2a: set = {res2 for res2,res3,_ in distances_2_3}
        res2b: set = {res2 for res1,res2,_ in distances_1_2}  
        res2_set: set = res2a.intersection(res2b)
        
        res3a: set = {res3 for res3,res1,_ in distances_3_1}
        res3b: set = {res3 for res2,res3,_ in distances_2_3}
        res3_set: set = res3a.intersection(res3b)
        
        # Now we should have the residue numbers of each given atom capable of forming a triangle within the given distance thresholds.
        # We can store these points in a tuple of Triad objects and return it:
        print(f"{res1}: {res1_set}")
        print(f"{res2}: {res2_set}")
        print(f"{res3}: {res3_set}")
                    
    def alignTo(self, other: 'PDBCIFFile'):
        pass
    
    def plot2DStructure(self):
        pass

    def get2DStructure(self):
        pass

    def openPyMol(self):
        pass

    def openChimera(self):
        pass
    
    def getStructSite(self): 
        pass

# Test triad Ser-His-Asp:
#my_file.findTriads('SER', 'HIS', 'ASP', 'OG', 'ND1', 'OD1', 1, 1, 6, 5, 4, 7)
   
class PdbCifFileCollection():
    """
    Class that represents a collection of PDBx/mmCIF files and methods to manipulate them.
    """

    def __init__(self, path_to_cif_collection: Path):

        ### DEFENSIVE CHECKS: ###
        if not isinstance(path_to_cif_collection, Path):
          pass

        if not path_to_cif_collection.is_dir():
            pass

        if not path_to_cif_collection.exists():
            pass

        ### INIT: ###
        self.full_path: Path = path_to_cif_collection.resolve()  # Full path to the cif collection.
        
        try:
            self.pdbcif_files: list = [PdbCifFile(child) for child in path_to_cif_collection.iterdir()]

        except ValueError:
            print('Some files in this collection are not .cif files. These will be excluded from the object.')
            self.pdbcif_files: tuple = (PdbCifFile(child) for child in path_to_cif_collection.iterdir() if child.suffix == '.cif')

        return None

    @property
    def size(self):
        return len(self.pdbcif_files)
    
    def writeSequencesToFasta(self, out_file: Path) -> None:
        """
        Writes all the amino acid sequences in this PDBxCIF file collection to a fasta file.
        Header of each entry is the name of the file.
        If multiple amino acid sequences are in the cif file, they will be written as:
            > <name>_<number>
        """
        result: dict = {}
        for pdbcif_file in self.pdbcif_files:
            try:
                aa_sequences: list = pdbcif_file.getAminoAcidSequences()
                if len(aa_sequences) > 1:
                    for i in range(len(aa_sequences)):
                        result[pdbcif_file.name + '_' + str((i+1))] = aa_sequences[i]
                elif len(aa_sequences) == 1:
                    result[pdbcif_file.name] = aa_sequences[0]
            except TypeError:
                print('Encounted TypeError, likely when parsing cif file using parsnip.')
                continue
    
        # Write reusults to file:
        with out_file.open('w') as f:
            for header, sequence in result.items():
                f.write('>' + header + '\n' + sequence + '\n')
        
        
            
    def summarize(self):
        pass

    def getNumberOfMutants(self):
        pass

    def getNumberOfLigandBound(self):
        pass

    def getLigands(self):
        pass

    def alignPairwise(self):
        pass

    def alignAll(self):
        pass

    def plotSimilarityNetwork(self):
        pass

    


