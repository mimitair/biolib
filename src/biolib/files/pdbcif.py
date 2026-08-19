# Default library imports:
import re
from pathlib import Path
import itertools
import warnings

# External libraries:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import biotite.structure.io.pdbx as pdbxio
import biotite.structure.io as strucio
import biotite.structure as struc
from biotite import DeserializationError
import networkx as nx

# Imports from this project:
from biolib.util import util

class PdbCifFile:
    """
    This class represents a PDBx/mmCIF file as described in https://mmcif.wwpdb.org/, and contains methods to parse and manipulate them.
    """
    def __init__(self, path_to_pdbcif: str) -> None:
        """
        Initializes the PdbCifFFile object.
        
        Input:
           - path_to_pdbcif: Path: path to the cif file.

        Returns:
           - None
        """
        ### DEFENSIVE CHECKS ###
        path: Path = Path(path_to_pdbcif)
        
        if not path.exists():
            raise FileNotFoundError(f"{path} does not exist. Check if you have provided the correct file path.")
            
        if not path.is_file():
            raise FileNotFoundError(f"{path} is not a file. Check if you have provided a file path instead of a directory.")

        ### INIT ###
        self.full_path: Path = path.resolve()  # Resolved path to the CIF file
        self.name: str = path.stem  # Name of the file without suffix and prefix

    
    def getPDBAccession(self) -> str:
        """Returns the PDB accession of this cif file based on the _entry.id column
        TODO: What if the _entry.id column is not there? What does biotite do? (e.g., predicted structures)
        """
        return self.toBiotiteCifFile().block['entry']['id'].as_item()

    def getUniprotAccession(self) -> str:
        """
        NOT YET IMPLEMENTED
        Returns the Uniprot ID that maps to the PDB ID found in the file based on Uniprot mapping service.
        TODO: what if no PDB ID found?
        """
        pass
    
    def countDataBlocks(self) -> int:
        """ Returns the total amount of data blocks in this cif file.
        """
        return self.full_path.read_text().count('#') -1

    def countLoopBlocks(self) -> int:
        """Returns the amount of loop data blocks in this cif file.
        """
        return self.full_path.read_text().count('loop_')

    
    def categoryExists(self, category: str) -> bool:
        """Returns true if the given category exists in the cif file based on regex matches.
        Works only for non-loop data blocks.

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
        """Converts this cif file to a biotite CifFile object

        Returns:
           - pdbxio.CIFFile: this cif file as a biotite CIFFile object
        """
        return pdbxio.CIFFile.read(self.full_path)

    def toBiotiteAtomArray(self) -> struc.AtomArray:
        """Convert the _atom_site data block to a biotite AtomArray object.
        The first model is chosen by default.
        The atom_id and charge annotoation categories are added by default.

        Returns:
           - struc.AtomArray: atom_site data as a biotite AtomArray object
        """
        atom_array: struc.AtomArray = strucio.load_structure(self.full_path, model=1, extra_fields=['atom_id', 'charge'])
        return atom_array

    def getChainCount(self) -> int:
        """Get the amount of chains in this cif file.
        """
        return struc.get_chain_count(self.toBiotiteAtomArray())

    def getModelCount(self) -> int:
        """Get the amount of models in this cif file.
        """
        return pdbxio.get_model_count(self.toBiotiteCifFile())
        
    def categoryToDf(self, category: str) -> pd.DataFrame():
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

    @staticmethod
    def findAtomPairs(atom_array: struc.AtomArray,
                      pair: tuple,
                      max_dist: float) -> struc.AtomArrayStack:
        """Find all atom pairs of two distinct residues based on a maximum distance threshold.
        For example; I want to find all serine-histidine pairs where the CA atoms are within 8.5A of each other.
        Note this only works when the two residues are non-identical. This will not work where res1_name == res2_name
        Neigbouring residues are omitted.
        NEEDS TESTING
        
        Input:
           - pair: tuple: Pair of atoms to search in the following format: ('res_name1@atom_name1', 'res_name2@atom_name2')
           - max_dist: float: The maximum distcance (in angstroms) between the first and second atom.
           
        Returns:
           list[struc.AtomArray]: A list of biotite AtomArray objects. Each AtomArray consists of two Atom objects comprising a pair within the given distance threshold.
        """
        # Inititate result list:
        result: list = []

        # Extract residue names and atom names from pair input:
        res_name1: str = pair[0].split('@')[0]
        atom_name1: str = pair[0].split('@')[1]
        res_name2: str = pair[1].split('@')[0]
        atom_name2: str = pair[1].split('@')[1]

        # Only consider the amino acids (no hetero atoms):
        atom_array = atom_array[struc.filter_amino_acids(atom_array)]
        
        # Consider each chain idnividually:
        for chain in struc.chain_iter(atom_array):
            # Filter based on res names and atom names, hetero atoms are discared by default:
            chain = chain[((chain.res_name==res_name1) & (chain.atom_name==atom_name1)) | ((chain.res_name==res_name2) & (chain.atom_name==atom_name2))]

            # Now get the index of each atom pair within the given distance threshold via cell list object for efficient distance calculations:
            idx_pairs: np.ndarray = struc.CellList(chain, cell_size=max_dist).get_atoms(chain.coord, radius=max_dist, result_format=struc.CellList.Result.PAIRS)

            # Omit the index pairs along the diagonal of the 'distance matrix':
            # Also filter out redundant pairs such as [[a,b][b,a]] by sorting and uniqueing along the appropriate axes
            idx_pairs_filtered: np.ndarray = np.unique(np.sort(idx_pairs[idx_pairs[:,0]!=idx_pairs[:,1]], axis=1), axis=0)

            # append atom array to result list if the residue names of the pairs are not identical and, not next to each other and in the same chain
            result.extend([struc.array([chain[i], chain[j]]) for i,j in idx_pairs_filtered if chain[i].res_name != chain[j].res_name and abs(chain[i].res_id - chain[j].res_id) != 1])

        # Return list of atom arrays:
        return result
    
    @staticmethod
    def findTriads(atom_array: struc.AtomArray,
                   pair1: tuple,
                   pair2: tuple,
                   pair3: tuple,
                   max_dist_pair1: float,
                   max_dist_pair2: float,
                   max_dist_pair3: float) -> list[struc.AtomArray]:
        """Find atoms of three distinct residues within a given distance of each other.
      
        Input:
           -

        Returns:
           - list[struc.AtomArray]: list of biotite AtomArray objects. Each AtomArray contains Three Atom objects that fulfill the distance thresholds.
        """
        # Initiate result list:
        result: list = []

        # Extract the atom names to filter later on:
        atom_names: list = [pair.split('@')[1] for pair in pair1+pair2+pair3]
        
        # Consider each chain individually:
        for chain in struc.chain_iter(atom_array):
            # Initiate graph:
            g: nx.Graph = nx.Graph()
            
            # Add the residue ID of each atom pair as an edge. F.e. 126-167 167-189 
            g.add_edges_from([chain.res_id for chain in PdbCifFile.findAtomPairs(chain, pair1, max_dist_pair1) + PdbCifFile.findAtomPairs(chain, pair2, max_dist_pair2) + PdbCifFile.findAtomPairs(chain, pair3,  max_dist_pair3)])

            # Now we have to couple each residue back to its Atom object based on the res_id.
            # Then use networkx all_triangles() function to find edges that form a triangle. This returns an iterator where each element contains the residue IDs of a triad.
            # We're using res_id instead of atom_id to allow two atoms in the same residue to complete the triangle (f.e.: two nitrogens in HIS)
            # Also only return the atoms that were specified in the input (otherwise all atoms of the residues are returned)
            result.extend([chain[(np.isin(chain.res_id, triad)) & (np.isin(chain.atom_name, atom_names))] for triad in nx.all_triangles(g)])

            # Clear the graph for the next chain:
            g.clear()

        # Return result as a list of atom arrays:
        return result

    @staticmethod
    def applySasaToBiotiteAtomArray(atom_array: struc.AtomArray) -> struc.AtomArray:
        atom_array.set_annotation('sasa', struc.sasa(atom_array))
        return atom_array
    
    def isMutant(self) -> bool:
        """Check if this cif file has mutated residues based on the _entity.pdbx_mutation column.

        Returns:
           - bool: True if any value besides '?' or '.' is encountered in the pdbx_mutation column 
        """
        try:
            column: pdbx.CifColumn = self.toBiotiteCifFile().block['entity']['pdbx_mutation']
            for value in self.toBiotiteCifFile().block['entity']['pdbx_mutation'].as_array():
                if value != '?' and value != '.':
                    return True
            return False
     
        except DeserializationError as e:
            print(self.name, e)
            return False
        
    def getHetero(self, omit_water: bool=True) -> np.ndarray:
        """Returns the residue names of all hetero atoms
        """
        atom_array: struc.AtomArray = self.toBiotiteAtomArray()
        hetero: struc.AtomArray = atom_array[atom_array.hetero==True]
        if omit_water:
            return np.unique(hetero[hetero.res_name!='HOH'].res_name)
        else:
            return np.unique(hetero.res_name)


    @staticmethod
    def atomArrayToDf(atom_array: struc.AtomArray) -> pd.DataFrame:
        result: dict = {}
        columns: list = atom_array.get_annotation_categories()
        for annotation in columns:
            result[annotation] = atom_array.get_annotation(annotation)
        return pd.DataFrame.from_dict(result)

    
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
        result: dict = {}

        for pdbcif_file in self.pdbcif_files:
            result[pdbcif_file.name] = pdbcif_file.getChainCount()

        return result
    
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

    def findTriads(self,
                   pair1: tuple,
                   pair2: tuple,
                   pair3: tuple,
                   max_dist_pair1,
                   max_dist_pair2,
                   max_dist_pair3) -> pd.DataFrame:
        """Detect triads in a cif file collection.
        Returns a dataframe where each row represents a detected triad.
        """
        # Initiate empty dataframe to store results
        result: pd.DataFrame = pd.DataFrame()
        
        # Start triad count at 1
        count = 1
        
        # Loop over each file 
        for cif_file in self.pdbcif_files:
            # Loop over every triad found in the file
            for triad in PdbCifFile.findTriads(cif_file.toBiotiteAtomArray(), pair1, pair2, pair3, max_dist_pair1, max_dist_pair2, max_dist_pair3):
                df: pd.DataFrame = PdbCifFile.atomArrayToDf(triad)
                # Add additional columns to track file name and triad id:
                df['file_name'] = cif_file.name
                df['triad_id'] = count
                count += 1 # increment
                # Concatenate to df:
                result = pd.concat([df, result], ignore_index=True)

        return result
             
    
    def getMutants(self):
        """Returns the file names of cif files with mutated residues.
        """
        return [cif_file.name for cif_file in self.pdbcif_files if cif_file.isMutant()]
        
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

    


