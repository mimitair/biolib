# Default library imports:
import sys
import re
from pathlib import Path
import logging
import itertools
import warnings

# External libraries:
import pandas as pd
from parsnip import CifFile
import numpy as np
import matplotlib.pyplot as plt
import biotite.structure.io.pdbx as pdbxio
import biotite.structure.io as strucio
import biotite.structure as struc

# Imports from this project:
from biolib.util import util

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
        
        ### DEFENSIVE CHECKS ###
        if not isinstance(path_to_pdbcif, Path):
            raise TypeError("Path to PDB/mmCIF file must be a Path object. Convert it to a Path object before initiating the class.")
        
        if not path_to_pdbcif.exists():
            raise FileNotFoundError(f"{path_to_pdbcif} does not exist. Check if you have provided the correct file path.")
            
        if not path_to_pdbcif.is_file():
            raise FileNotFoundError(f"{path_to_pdbcif} is not a file. Check if you have provided a file path instead of a directory.")

        if not path_to_pdbcif.suffix == '.cif':
            raise ValueError(f"{path_to_pdbcif} is not a .cif file. This class only accepts .cif files.")
        
        ### INIT ###
        self.full_path: Path = path_to_pdbcif.resolve()  # Resolved path to the CIF file
        self.name: str = path_to_pdbcif.stem  # Name of the file without suffix and prefix

    
    def getAccession(self) -> str:
        """Returns the PDB accession of this cif file based on the _entry.id column
        """
        return self.toBiotiteCifFile().block['entry']['id'].as_item()

    def countDataBlocks(self) -> int:
        """ Returns the total amount of data blocks in this cif file.
        """
        return self.full_path.read_text().count('#') -1

    
    def countLoopBlocks(self) -> int:
        """Returns the amount of loop data blocks in this cif file.
        """
        return self.full_path.read_text().count('loop_')

    
    def categoryExists(self, category: str) -> bool:
        """Returns true if the given category exists in the cif file based on regex matches. Works only for non-loop data blocks.

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

    def toBiotiteCifFile(self) -> pdbxio.CIFFile:
        return pdbxio.CIFFile.read(self.full_path)

    def toBiotiteAtomArray(self):
        return strucio.load_structure(self.full_path)

    @classmethod
    def biotiteAtomArrayToDf(atom_array: struc.AtomArray) -> pd.DataFrame:
        pass
    
    @classmethod
    def filterAtomArrayResidues(atom_array: struc.AtomArray, res_names: list) -> struc.AtomArray:
        """Filter a given atom array based on residue names.
        """
        mask = np.isin(atom_array.res_name, res_names)
        return atom_array[mask]

    @classmethod
    def filterAtomArrayAtoms(atom_array: struc.AtomArray, atom_names: list) -> struc.AtomArray:
        """ Filter a given atom array based on atom names.
        """
        mask = np.isin(atom_array.atom_name, atom_names)
        return atom_array[mask]
    
    def getChainCount(self) -> int:
        return struc.get_chain_count(self.toBiotiteAtomArray())

    def getModelCount(self) -> int:
        return pdbxio.get_model_count(self.toBiotiteCifFile())
        
    def categoryToDf(self, category: str) -> pd.DataFrame | dict: 
        """
        Convert any data block, given its category name, to a pandas dataframe or dictionary.
        This function uses the parnsip external library.

        Input:
            - category: str: The category to extract (f.e.: "atom_site")

        Returns:
            - pd.DataFrame: A loop data block as a dataframe. Beware that no further processing is done. Every value is essentially a string.
            - dict: A non-loop data block as a dict.

        Raises:
            - ValueError: If the given category cannot be found in the cif file.
        """
        if self.categoryExists(category):
            return {key:value for key,value in self.toParsnip().pairs.items() if category + '.' in key}
        
        elif self.loopCategoryExists(category):
            return pd.DataFrame([arr for arr in self.toParsnip().loops if category + '.' in arr.dtype.names[0]][0].reshape(-1))

        else:
            raise ValueError(f'{category} does not exist in {self.name}')

    def categoryToDf2(self, category: str) -> pd.DataFrame():
        """Convert a given data block to a pandas dataframe.
        NEEDS TESTING
        """
        # Convert to biotite cif file object:
        cif_biotite = self.toBiotiteCifFile()
        # First we have to get all the column names of this category
        keys: list = [key for key in cif_biotite.block[category].keys()]
        values: list = [cif_biotite.block[category][key].as_array() for key in keys]
        data: dict = {key: value for (key,value) in zip(keys, values)}
        
        return pd.DataFrame(data=data)

    def findAtomPairsWithinDistance(self,
                                    res_name1: str,
                                    res_name2: str,
                                    atom_name1: str,
                                    atom_name2: str,
                                    max_dist: float) -> list:
        """
        Find all atom pairs of two different residues based on a maximum distance threshold.
        For example; I want to find all serine-histidine pairs where the CA atoms are within 8.5A of each other.
        Note this only works when the two residues are non-identical. This will not work where res1_name == res2_name
        NEEDS TESTING
        
        Input:
           -
        Returns:
           list: A list of biotite AtomArray objects. Each AtomArray consists of two Atom objects comprising a pair within the given distance threshold.
        """
        # Convert to atom array:
        atom_array: struc.AtomArray = self.toBiotiteAtomArray()

        # Filter based on res names and atom names, hetero atoms are discared by default:
        atom_array = atom_array[((atom_array.res_name==res_name1) & (atom_array.atom_name==atom_name1)) | ((atom_array.res_name==res_name2) & (atom_array.atom_name==atom_name2)) & atom_array.hetero==False]

        # Now get the index of each atom pair within the given distance threshold via cell list object for efficient distance calculations:
        idx_pairs: np.ndarray = struc.CellList(atom_array, cell_size=max_dist).get_atoms(atom_array.coord, radius=max_dist, result_format=struc.CellList.Result.PAIRS)

        # Omit the index pairs along the diagonal of the 'distance matrix':
        # Also filter out redundant pairs such as [[a,b][b,a]] by sorting and uniqueing along the appropriate axes
        idx_pairs_filtered: np.ndarray = np.unique(np.sort(idx_pairs[idx_pairs[:,0]!=idx_pairs[:,1]], axis=1), axis=0)

        # Return as a list of atom arrays if the residue names of the pairs are not identical:
        return [struc.array([atom_array[i], atom_array[j]]) for i,j in idx_pairs if atom_array[i].res_name != atom_array[j].res_name]

    def findTriads3(self,
                    res_name1: str,
                    res_name2: str,
                    res_name3: str,
                    max_dist1_2: float,
                    max_dist2_3: float,
                    max_dist3_1: float,
                    atom_name1: str = "CA",
                    atom_name2: str = "CA",
                    atom_name3: str = "CA") -> list:

        # Get the residue IDs for all pairs first
        pairs1_2: list = [atom_array.res_id for atom_array in self.getAtomPairsWithinDistance(res_name1, res_name2, atom_name1, atom_name2, max_dist1_2)]
        pairs2_3: list = [atom_array.res_id for atom_array in self.getAtomPairsWithinDistance(res_name2, res_name3, atom_name2, atom_name3, max_dist2_3)]
        pairs3_1: list = [atom_array.res_id for atom_array in self.getAtomPairsWithinDistance(res_name3, res_name1, atom_name3, atom_name1, max_dist3_1)]

        # convert to numpy array:
        a = np.array([pairs1_2, pairs2_3, pairs3_1])
        
        result: list = []
        
        # Now filter those that actually make a triad
        for atom_array in pairs_1_2:
            
            

        
    def findTriads(self,
                   res1_name: str,
                   res2_name: str,
                   res3_name: str,
                   max_distance1_2: float,
                   max_distance2_3: float,
                   atom1_name: str = "CA",
                   atom2_name: str = "CA",
                   atom3_name: str = "CA") -> list:
        """Detect a triad of residues based on distance threshold between their atoms.
        Needs testing
        """
        # First we need the current cif file as an AtomArray object:
        atom_array: AtomArray = self.toBiotiteAtomArray()
        
        # Then we extract the atoms for each residue:
        res1_atoms: AtomArray = atom_array[(atom_array.res_name==res1_name) & (atom_array.atom_name==atom1_name)]
        res2_atoms: AtomArray = atom_array[(atom_array.res_name==res2_name) & (atom_array.atom_name==atom2_name)]
        res3_atoms: AtomArray = atom_array[(atom_array.res_name==res3_name) & (atom_array.atom_name==atom3_name)]

        result: list = []
        # Handle each chain independently:
        for chain in np.unique(atom_array.chain_id):
            chain_res1_atoms = res1_atoms[res1_atoms.chain_id==chain]
            chain_res2_atoms = res2_atoms[res2_atoms.chain_id==chain]
            chain_res3_atoms = res3_atoms[res3_atoms.chain_id==chain]

            # Then we define two lists that keep track of the atoms that pass the distance threshold:
            distance1_2_pass: list = []
            distance2_3_pass: list = []

            # Now we loop over the atoms of res1 and res2:
            for atom1 in chain_res1_atoms:
                for atom2 in chain_res2_atoms:
                    # Set distance threshold and omit neighboring residues.
                    if struc.distance(atom1, atom2) < max_distance1_2 and abs(atom1.res_id - atom2.res_id) != 1:
                        distance1_2_pass.append((atom1, atom2))

            # Now we loop over the atoms of res2 that passed the distance threhsold with res1:
            for _,atom2 in distance1_2_pass:
                for atom3 in chain_res3_atoms:
                    if struc.distance(atom2, atom3) < max_distance2_3 and abs(atom2.res_id - atom3.res_id) != 1:
                        distance2_3_pass.append((atom2, atom3))

            # Now we complete the threesome. If atom2 passes both distance thresholds, it completes the triad:
            for atom1,atom2a in distance1_2_pass:
                for atom2b,atom3 in distance2_3_pass:
                    if atom2a.res_id == atom2b.res_id:
                        result.append(struc.array([atom1, atom2a, atom3]))

        # Return the result as a list of biotite atom arrays:
        return result
                                
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
        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                return CifFile(self.full_path)
        
            except Warning as w:
                warnings.filterwarnings("ignore")
                #print(f"Warning in {self.name}")
                return CifFile(self.full_path)

    
    def getAminoAcidSequences(self) -> list:
        """
        Returns the amino acid sequences of each polymer entity in this file as a list of strings.
        If the cif file has only one polymer entity, a list of length 1 is returned.

        Returns:
            - list: A list of the amino acid sequences for each polymer entity in this PDBx/mmCIF file. Empty list if no polypeptide sequences were found
        """
        # Extract the data block where we can find the amino acid sequences:
        df: pd.DataFrame | dict = self.categoryToDf('_entity_poly')

        # TODO checking type every time is too slow?
        if isinstance(df, pd.DataFrame):
            # Extract only the rows where _entity.type == 'polypeptide' 
            df = df.loc[(df['_entity_poly.type'] == "'polypeptide(L)'") | (df['_entity_poly.type'] == "polypeptide(L)") | (df['_entity_poly.type'] == "'polypeptide(D)'") | (df['_entity_poly.type'] == "polypeptide(D)")]
            return [seq.replace('\n', '').replace(';', '') for seq in df['_entity_poly.pdbx_seq_one_letter_code'].tolist()]

        elif isinstance(df, dict):
            if (df['_entity_poly.type'] == "'polypeptide(L)'") or (df['_entity_poly.type'] == "polypeptide(L)") or (df['_entity_poly.type'] == "polypeptide(D)") or (df['_entity_poly.type'] == "'polypeptide(D)'"):
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
        """ Counts the amount of polymer entities labeled as 'polypeptide(L)' or 'polypeptide(D)'

        Returns:
           - int: Amount of polypetide entities in this cif file.
        """
        
        if self.categoryExists('_entity_poly'):
            df: dict = self.categoryToDf('_entity_poly')
            assert isinstance(df, dict)
            if df['_entity_poly.type'] == "'polypeptide(D)'" or df['_entity_poly.type'] == "'polypeptide(L)'" or df['_entity_poly.type'] == "polypeptide(D)" or df['_entity_poly.type'] == "polypeptide(L)":
                return 1
            else:
                return 0

        elif self.loopCategoryExists('_entity_poly'):
            df: pd.DataFrame = self.categoryToDf('_entity_poly')
            assert isinstance(df, pd.DataFrame)
            df = df.loc[(df['_entity_poly.type'] == "polypeptide(D)") | (df['_entity_poly.type'] == "polypeptide(L)") | (df['_entity_poly.type'] == "'polypeptide(D)'") | (df['_entity_poly.type'] == "'polypeptide(L)'")]
            return len(df)

        else:
            print('_entity_poly does not exist')
            return 0
                        
    def atomSiteToDf(self, filter: dict=None) -> pd.DataFrame:
        """ Returns the _atom_site loop block as a dataframe.
        This method is preferred over categoryToDf(), since the columns are converted to appropriate dtypes for memory efficiency.
        Additionally, a filter can be applied to select only the desired rows.
        TODO: check if all necessary columns are present

        Input:
           -filter: dict: Keys are column names, values are lists of values allowed for that columns (e.g., {'_atom_site.group_PDB': 'ATOM'})

        Returns:
           -pd.DataFrame: The _atom_site loop block as a dataframe. Possibly with some rows omitted due to the filter argument
        """

        # Get the atom site dataframe:
        df = self.categoryToDf('_atom_site')

        # Dictionary that maps column names to desired dtypes:
        dtype_mapping: dict = {"_atom_site.id": "uint32",
                               "_atom_site.Cartn_x": "float16",
                               "_atom_site.Cartn_y": "float16",
                               "_atom_site.Cartn_z": "float16",
                               "_atom_site.pdbx_PDB_model_num": "uint16",
                               "_atom_site.auth_seq_id": "uint16",
                               "_atom_site.label_seq_id": "uint16",
                               "_atom_site.label_entity_id": "uint16"}

            
        # Return with adapted typing:
        return util.safeCastColumns(df, dtype_mapping)

        
    ##################################################
    ##### EVERYTHING UNDERNEATH IS NOT FUNCTIONAL ####
    ##################################################

    def filterAtomSite(self, residue_names: list, atom_names: list) -> pd.DataFrame:
        """
        Filter the atom site loop block based on residue names and atom names.
        Only those passed to the function will be retained and returned as a dataframe.
        TODO: handle empty lists (retain everything)
        """
        df = self.atomSiteToDf()

        return df.loc[(df['_atom_site.label_comp_id'].isin(residue_names)) & (df['_atom_site.label_atom_id'].isin(atom_names))]

    
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
        
                    
    def alignTo(self, other: 'PDBCIFFile'):
        pass
    
    def plot2DStructure(self):
        pass

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
    def count(self):
        """Returns the amount of cif files in this collection.
        """
        return len(self.pdbcif_files)
    
    def writeSequencesToFasta(self, out_file: Path) -> dict:
        """
        Writes all the amino acid sequences in this PDBxCIF file collection to a fasta file.
        Header of each entry is the name of the file.
        If multiple amino acid sequences are in the cif file, they will be written as:
            > <name>_<number>

        Input:
           - out_file: Path: The file path to which the amino acid sequences will be written

        Returns:
           -dict: Dictionary containing {'name_of_cif_file': [<list of sequences>]}
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
            except TypeError as e:
                print(f'Encountered TypeError when reading from {pdbcif_file.name}, likely when parsing cif file using parsnip.')
                print(e)
                continue
            except ValueError as e:
                print(f'Encountered ValueError when reading from {pdbcif_file.name}')
                print(e)
                continue
    
        # Write reusults to file:
        with out_file.open('w') as f:
            for header, sequence in result.items():
                f.write('>' + header + '\n' + sequence + '\n')

        # Return the dictionary for testing purposes:
        return result

    def getModelCounts(self) -> tuple:
        result = [pdbcif_file.getModelCount() for pdbcif_file in self.pdbcif_files]
        return tuple(result)
    
    def getChainCounts(self) -> tuple:
        result: list = []

        for pdbcif_file in self.pdbcif_files:
            result.append(pdbcif_file.getChainCount())

        return tuple(result)
    
    def getPolypeptideCounts(self) -> tuple:
        result: list = []
        
        for pdbcif_file in self.pdbcif_files:
            result.append(pdbcif_file.countPolypeptideEntities())

        return tuple(result)
    
    def plotPolypeptideCount(self) -> None:
        """ Plots the distirbution of polypetide counts for each CIF file in this collection as a histogram.
        """
        data: tuple = self.getPolypeptideCounts()
     
        #TODO fig, ax?
        fig, ax = plt.subplots()
        
        ax.hist(data)
        plt.show()

        return None

    def getPDBIDs(self) -> list:
        result: list = []
        for pdbcif_file in self.pdbcif_files:
            result.append(pdbcif_file.getPDBID())

        return result
    
    def writePDBIDsToFile(self, out_file: Path) -> None:
        """Write all the PDB accessions of the collection to a .txt file.
        Each PDB accession is placed on a new line

        Input:
           - out_file: Path: File path to write the output to
        """
        result: list = self.getPDBIDs()
        
        with out_file.open("w") as f:
            f.write('\n'.join(result))

        return None

        
    def toDf(self):
        data: dict = {
            "file_name": [cif_file.name for cif_file in self.pdbcif_files],
            "pdb_id": self.getPDBIDs(),
            "polypeptide_count": self.getPolypeptideCounts(),
            "chain_count": self.getChainCounts(),
            "model_count": self.getModelCounts()
        }
        df: pd.DataFrame = pd.DataFrame.from_dict(data)
        return df

    def getAuthors(self):
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

    


