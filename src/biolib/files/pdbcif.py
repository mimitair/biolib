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

class CifFile:
    """
    This class represents a PDBx/mmCIF file as described in https://mmcif.wwpdb.org/, and contains methods to parse and manipulate them.
    """
    def __init__(self, path_to_pdbcif: str) -> None:
        """
        Initializes the CifFFile object.
        
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
        self.name: str = path.name  # Name of the file without prefix
        self.stem: str = path.stem  # Name without prefix and suffix

    
    def getID(self) -> str:
        """Returns the value in the _entry.id column
        If it does not exist, return the the part before the first '_' character of stem of the file name (e.g.; 'abcd' for 'abcd_full_A.cif')
        """
        try:
            return self.toBiotiteCifFile().block['entry']['id'].as_item().strip().upper()
        
        except KeyError as e: # When the 'id' column cannot be found
            return self.stem.split('_')[0].strip().upper()
    
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
        try:
            return strucio.load_structure(self.full_path, include_bonds=True, model=1, extra_fields=['atom_id', 'charge', 'b_factor', 'occupancy'])
        
        except KeyError as e:
            mycif: pdbxio.CIFFile = self.toBiotiteCifFile()
            mycif.block['atom_site']['pdbx_PDB_model_num'] = np.ones(len(mycif.block['atom_site']['id'].as_array()), dtype=np.int32)
            return pdbxio.get_structure(mycif, include_bonds=True, model=1, extra_fields=['atom_id', 'charge', 'b_factor', 'occupancy'])

        
    def getChainCount(self) -> int:
        """Get the amount of chains in this cif file.
        """
        return struc.get_chain_count(self.toBiotiteAtomArray())

    
    def getUniqueChainCount(self) -> int:
        """Returns the amount of uniqe chain identifiers in this cif file
        """
        return np.unique(struc.get_chains(self.toBiotiteAtomArray())).size

    
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
                      max_dist: float) -> list[struc.AtomArray]:
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
            
            # We could get an empty chain here (none of the residues are present)
            if len(chain) != 0:
                # Now get the index of each atom pair within the given distance threshold via cell list object for efficient distance calculations:
                idx_pairs: np.ndarray = struc.CellList(chain, cell_size=max_dist).get_atoms(chain.coord, radius=max_dist, result_format=struc.CellList.Result.PAIRS)

                # Omit the index pairs along the diagonal of the 'distance matrix':
                # Also filter out redundant pairs such as [[a,b][b,a]] by sorting and uniqueing along the appropriate axes
                idx_pairs_filtered: np.ndarray = np.unique(np.sort(idx_pairs[idx_pairs[:,0]!=idx_pairs[:,1]], axis=1), axis=0)

                # append atom array to result list if the residue names of the pairs are not identical and, not next to each other and in the same chain
                result.extend([struc.array([chain[i], chain[j]]) for i,j in idx_pairs_filtered if chain[i].res_name != chain[j].res_name and abs(chain[i].res_id - chain[j].res_id) != 1])
            else:
                continue

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
            g.add_edges_from([chain.res_id for chain in CifFile.findAtomPairs(chain, pair1, max_dist_pair1) + CifFile.findAtomPairs(chain, pair2, max_dist_pair2) + CifFile.findAtomPairs(chain, pair3,  max_dist_pair3)])

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
        """Returns the set of residue names of all hetero atoms.
        Omits 'HOH' (water molecules) by default.

        Input:
           - omit_water:bool: Whether to omit 'HOH' molecules from the output [True]

        Returns:
           - np.ndarray: numpy array containing hetero atom names as strings
        """
        atom_array: struc.AtomArray = self.toBiotiteAtomArray()
        hetero: struc.AtomArray = atom_array[atom_array.hetero==True]
        if omit_water:
            return np.unique(hetero[hetero.res_name!='HOH'].res_name)
        else:
            return np.unique(hetero.res_name)


    @staticmethod
    def atomArrayToDf(atom_array: struc.AtomArray) -> pd.DataFrame:
        """ Convert an atom array to a pandas dataframe.
        """
        result: dict = {}
        columns: list = atom_array.get_annotation_categories()
        for annotation in columns:
            result[annotation] = atom_array.get_annotation(annotation)

        result['x_coord'] = atom_array.coord[:,0]
        result['y_coord'] = atom_array.coord[:,1]
        result['z_coord'] = atom_array.coord[:,2]
            
        return pd.DataFrame.from_dict(result)


    def getEnzymmCatalyticSite(self, path_to_enzymm_out: str) ->  tuple | None:
        """After running enzymm on a collection of cif files, find the catalytic site associated with this entry ID.
        Note that it is assumed that the enzymm output has been filtered by RMSD beforehand (see EnzymmOut.filterByRMSD())
        Input:
           - path_to_enzymm_out: str: path to a .tsv file containing enzymm output (filtered by RMSD such that only one row per query_id is present)

        Returns:
           - tuple(struc.AtomArray, list[tuple(res_name,chain_id,res_id)])
           - None: if this query ID cannot be found in the enzymm output (e.g., no catalytic motif was found)
        """
        def parseMatchedResidues(s: str) -> list:
            # Helper function to parse the 'matched_residues' column of enzymm output
            # Returns [(res_name,chain_id,res_id),(res_name,chain_id,res_id),...]
            return [tuple(residue.split('_')) for residue in s.split(',')]

        # Read the enzymm output:
        df: pd.DataFrame = pd.read_csv(path_to_enzymm_out, sep='\t', comment='#')
        try:
            # Locate the row where this ID matches the query_id column in enzymm:
            row = df.loc[df['query_id']==self.getID().upper()].squeeze() # squeeze to force it into 1D (for some reason it returns a dataframe instead of a series)
            # Apparently row can be empty here (normally it should then raise a KeyError but pandas has its own ways :))
            if not row.empty:
                # Parse the matched_residues column using helper function:
                matched_residues: list = parseMatchedResidues(row['matched_residues'])
                # Now extract the residues as an atom array:
                atom_array: struc.AtomArray = self.toBiotiteAtomArray()
                atom_array = atom_array[(np.isin(atom_array.res_id, [int(res[2]) for res in matched_residues])) & (atom_array.chain_id==matched_residues[0][1])] # assuming all residues are in the same chain for chain_id mask
                return (atom_array, matched_residues)
            else:
                return None

        except KeyError as e: # when pandas does not find the row
            return None
        
    
    def getNearestPocket(self, coord: np.ndarray, path_to_fpocket_out: str) -> tuple:
        """Find the pocket closest to the given coord based on fpocket output

        Input:
           - coord: np.ndarray: [x y z] coordinates as used by biotite
           - path_to_fpocket_out: str: A directory containing the fpocket output for this cif file.

        Returns:
           - tuple(pocket_id: int, distance_to_pocket: float)
        """
        # Import here to avoid circular imports:
        from biolib.files.fpocket_out import FpocketOut
        
        # Obtain the centroids of all pockets from the fpocket output:
        pocket_centroids: dict = FpocketOut(path_to_fpocket_out).getPocketCentroids()

        # Calculate the distance from coord to each centroid:
        distances: dict = {pocket_id: struc.distance(coord, pocket_coord) for pocket_id, pocket_coord in pocket_centroids.items() }

        # This should return the key (pocket_id) and value (distance) where the value is lowest as a tuple (id,distance):
        return min(distances.items(), key=lambda x: x[1]) # wtf does lambda x: x[1] do??

    
    def runFpocket(self, args: str, out_dir: str|Path) -> None:
        """Run fpocket on this cif file.
        """
        command: list = ["fpocket", "-f", f"{self.full_path}"]
        if args: # If non-empty arguments are given
            command.extend(args.split(' '))
        
        from subprocess import run
        out = run(command, check=True, capture_output=True)
        print(out)

        # This should be the path that fpocket creates upon execution. In the same directory as where the input file was
        fpocket_out_path: Path = self.full_path.parent / f"{self.full_path.stem}_out"
        # Check if it exists, then move it to the desired output folder
        if fpocket_out_path.exists():
            out_dir: Path = Path(out_dir)
            out_dir.mkdir(exist_ok=True)
            fpocket_out_path.move_into(out_dir)
            return None

        else:
            print(f"could not find fpocket out for {self.name}")
            return None


    def chainToCifFile(self, chain_id: str) -> pdbxio.CIFFile:
        """
        Returns a CIFFile object only including the specified chain (most of the dictionaries are lost, only atom_site and chem_comp...)
        """
        atom_array: struc.AtomArray = self.toBiotiteAtomArray()
        atom_array = atom_array[(struc.filter_amino_acids(atom_array)) & (atom_array.chain_id == chain_id)]

        # Initiate new CIFFile object with the same ID:
        categories: pdbxio.CIFCategory = pdbxio.CIFCategory({'id':self.getID()})
        block: pdbxio.CIFBlock = pdbxio.CIFBlock({'entry':categories})
        cif_file: pdbxio.CIFFile = pdbxio.CIFFile({self.getID():block})

        # Write the atom array to the cif file and return:
        pdbxio.set_structure(cif_file, atom_array)
        return cif_file


    def dockToLigand(self,
                     ligand: struc.AtomArray,
                     docking_coord: np.ndarray,
                     search_space: list,
                     path_to_vina_bin: str|Path):
        """
        """
        import biotite.application.autodock as autodock
        # Prepare the receptor:
        receptor: struc.AtomArray = self.toBiotiteAtomArray()
        receptor = receptor[struc.filter_amino_acids(receptor)]
        receptor.charge = struc.partial_charges(receptor)  # Adds Gasteiger charges

        # Start the autodock vina app
        app = autodock.VinaApp(ligand, receptor, docking_coord, search_space, bin_path=path_to_vina_bin)
        
        # Initialize some parameters:
        app.set_seed(0)
        app.set_cpu(1)
        app.set_max_number_of_models(100)
        app.set_energy_range(100.0)

        # Start docking run
        app.start()
        app.join()

        # Get docking coordinates for each binding mode
        docked_coord: np.ndarray = app.get_ligand_coord()
        # Create an AtomArrayStack for all docked binding modes
        docked_ligand: struc.AtomArrayStack = struc.from_template(ligand, docked_coord)
        # As Vina discards all nonpolar hydrogen atoms, their respective coordinates are NaN -> remove these atoms
        docked_ligand = docked_ligand[..., ~np.isnan(docked_ligand.coord[0]).any(axis=-1)]
        
        # Get energies for each binding pose and add to the atom array stack
        energies: np.ndarray = app.get_energies()

        return docked_ligand, energies

    
class CifFileCollection():
    """
    Class that represents a collection of PDBx/mmCIF files and methods to manipulate them.
    """

    def __init__(self, path_to_cif_collection: str|Path):

        path: Path = Path(path_to_cif_collection)
        ### DEFENSIVE CHECKS: ###
        if not isinstance(path, Path):
            pass

        if not path.is_dir():
            pass

        if not path.exists():
            pass

        ### INIT: ###
        self.full_path: Path = path.resolve()  # Full path to the cif collection.
        self.name: str = path.name
        
        return None

    
    @property
    def size(self):
        """Returns the amount of cif files in this collection.
        """
        return len([f for f in self.iterFiles()])

    
    def merge(self, other):
        pass

    
    def iterFiles(self):
        """
        Returns an iterator over every cif file in this collection.
        TO BE IMPLEMENTED
        """
        for child in self.full_path.iterdir():
            if child.suffix=='.cif':
                yield CifFile(child)
        
    
    def writeSequencesToFasta(self, out_file: Path) -> dict:
        """
        DEPRECTATED, NEEDS NEW IMPLEMENTAION
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
        for pdbcif_file in self.iterFiles():
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
        result = [pdbcif_file.getModelCount() for pdbcif_file in self.iterFiles()]
        return tuple(result)
    

    def getChainCounts(self) -> tuple:
        result: dict = {}

        for pdbcif_file in self.iterFiles():
            result[pdbcif_file.name] = pdbcif_file.getChainCount()

        return result


    def getUniqueChainCounts(self) -> dict:
        """Returns dict of {cif_file_name : unique_chain_counts}
        """
        return {cif.name:cif.getUniqueChainCount() for cif in self.iterFiles()}


    def plotUniqueChainCounts(self) -> None:
        data: dict = self.getUniqueChainCounts()
        counts: np.ndarray = np.array(data.values())

        fig,ax = plt.subplots()

        ax.hist(data)

        plt.show()

        return None
    

    def getPolypeptideCounts(self) -> tuple:
        result: list = []
        
        for pdbcif_file in self.iterFiles():
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


    def getIDs(self) -> list:
        return [cif.getID() for cif in self.iterFiles()]

    
    def writeIDsToFile(self, out_file: str) -> None:
        """Write all the IDs of the collection to a .txt file.
        Each ID is placed on a new line

        Input:
           - out_file: Path: File path to write the output to
        """
        result: list = self.getIDs()
        
        with Path(out_file).open("w") as f:
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
        for cif_file in self.iterFiles():
            # Loop over every triad found in the file
            for triad in CifFile.findTriads(cif_file.toBiotiteAtomArray(), pair1, pair2, pair3, max_dist_pair1, max_dist_pair2, max_dist_pair3):
                df: pd.DataFrame = CifFile.atomArrayToDf(triad)
                # Add additional columns to track file name and triad id:
                df['file_name'] = cif_file.name
                df['triad_id'] = count
                count += 1 # increment
                # Concatenate to df:
                result = pd.concat([df, result], ignore_index=True)

        # Sometimes the same triad is found across different chains. We remove these dupliactes here:
        filtered_indices: pd.Index = result.drop_duplicates(subset=['res_id', 'atom_name', 'file_name', 'res_name', 'element']).index
        
        return result.iloc[filtered_indices]
             

    def getMutants(self):
        """Returns the file names of cif files with mutated residues.
        """
        return [cif_file.name for cif_file in self.iterFiles() if cif_file.isMutant()]


    def filterMutants(self, out_dir: str) -> list:
        """
        Returns the file names of mutants.
        Creates a directory (out_dir)  with symlinks to the non-mutants.
        """
        # Make the output dir if it does not exist already:
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)

        mutants: list = []
        for cif_file in self.iterFiles():
            # If it's a mutant, add to mutants list
            if cif_file.isMutant():
                mutants.append(cif_file.name)
            # Otherwise make a symbolic link in the new directory:
            else:
                path: Path = out_dir / cif_file.name
                path.symlink_to(cif_file.full_path)

        return mutants

    
    def filterByChainCount(self, max_chain_count: int, out_dir: str = None) -> list:
        """Returns the file names of cif files exceeding the chain count
        Creates a directory (out_dir) containing with symlinks to only the cif files that have less chains than max_chain_count
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)
        
        outliers: list = []
        for cif in self.iterFiles():
            if cif.getUniqueChainCount() > max_chain_count:
                outliers.append(cif.name)
            else:
                if out_dir is not None:
                    path: Path = out_dir / cif.name
                    path.symlink_to(cif.full_path)
                else:
                    continue

        return outliers


    def getEnzymmCatalyticSites(self, path_to_enzymm_out: str) -> dict:
        """Returns dictionary {cif_file_name: (catalytic_site_atom_array, parsed_matched_residues)}
        """
        return {cif.name:cif.getEnzymmCatalyticSite(path_to_enzymm_out) for cif in self.iterFiles()}


    def runFpocket(self, args: str, out_dir: str|Path) -> None:
        """Run fpocket for each file in this collection.
        """
        for cif in self.iterFiles():
            cif.runFpocket(args, out_dir)

        return None
    
   
    def getCatalyticSitePockets(self, path_to_enzymm_out: str, path_to_fpocket_out_collection: str) -> pd.DataFrame:

        from biolib.files.fpocket_out import FpocketOutCollection, FpocketOut
        fpocket_out_coll: FpocketOutCollection = FpocketOutCollection(path_to_fpocket_out_collection)

        df_dict: dict = {}
        for cif in self.iterFiles():
            catalytic_site: np.ndarray = cif.getEnzymmCatalyticSite(path_to_enzymm_out)
            if catalytic_site is not None:
                # Calculate the centroid of the identified catalytic site
                catalytic_site_coord: np.ndarray = struc.centroid(catalytic_site[0])
                # Find the fpocket_out corresponding to this id:
                try:
                    fpocket_out: FpocketOut = FpocketOut(fpocket_out_coll.getFpocketOutFromID(cif.getID()))
                except IndexError as e: # the ID does not match:
                    fpocket_out: FpocketOut = FpocketOut(fpocket_out_coll.getFpocketOutFromID(cif.stem))
                # Get the nearest pocket:
                nearest_pocket: tuple = cif.getNearestPocket(catalytic_site_coord, fpocket_out.full_path)
                pocket_info: dict = fpocket_out.getPocketInfoAt(nearest_pocket[0])
                # Add the distance as well:
                pocket_info['distance_to_catalytic_site'] = nearest_pocket[1]
                pocket_info['pocket_id'] = nearest_pocket[0]
                # Now add to the master dictionary:
                df_dict[cif.getID()]=pocket_info

            else:
                continue

        return pd.DataFrame.from_dict(df_dict, orient='index')

    
    def getLigands(self):
        ligands: list = []
        for cif_file in self.iterFiles():
            ligands.extend(list(cif_file.getHetero()))
            
        return ligands


    def extractFoldseekClusters(self, path_to_foldseek_cluster: str, out_dir: str) -> None:
        """Makes a new directory containing only the foldseek cluster representatives in the adjacency matrix from foldseek (which might also have been length filtered separately).
        This function will only retain the appropriate chains that foldseek mentions (use foldseek with chain-name-mode 1 !!)
        
        Input:
           - path_to_foldseek_cluster: str: Usually a .tsv file
           - out_dir: str: Directory to which the 
        """
        # Prepare the output dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)

        # read cluster file as pandas dataframe
        df_adjacency: pd.DataFrame = pd.read_csv(path_to_foldseek_cluster, names=['representative', 'node'])
        # Extract the list of representatives, their IDs (first part of the string) and chain IDs (last part of the string)
        rep_list: dict = {s.split('_')[0]:s.split('_')[-1] for s in list(df_adjacency['representative'].unique())}
        
        for cif in self.iterFiles():
            file_id: str = cif.stem.split('_')[0] # This is the first part of the file used to match it to the cluster list from foldseek
            if file_id in list(rep_list.keys()):
                chain_id: str = rep_list[file_id]
                cif_file: pdbxio.CIFFile = cif.chainToCifFile(chain_id)
                file_name: str = cif.stem + f'_{chain_id}.cif'
                out_path: Path = out_dir / file_name
                cif_file.write(out_path)
            
        return None

    
    def removeIDs(self, path_to_ids: str, out_dir: str = None) -> list:
        """Removes entries with ids (PDB accessions/AF accessions) in the given list.
        Makes new directory (out_dir) with symlinks
        Returns the IDs of cif files that were removed
        """
        # Retrieve the list of accessions to be removed:
        accessions: list = []
        with Path(path_to_ids).open('r') as f:
            for line in f:
                accessions.append(line.strip())

        # Make the output directory:
        out_dir: Path = Path(out_dir)
        out_dir.mkdir(exist_ok=True)

        # Initiate list that will hold the IDs of the files that were removed
        removed: list = []
        # Loop over the cif files:
        for cif in self.iterFiles():
            # Remove if it is in the list:
            if (cif.getID() in accessions) or (cif.stem in accessions):
                removed.append(cif.getID())
            # Create symlink to og file in out_dir:
            else:
                if out_dir is not None:
                    path: Path = out_dir / cif.name
                    path.symlink_to(cif.full_path)

        return removed
    
    def alignPairwise(self):
        pass

    def alignAll(self):
        pass

    def plotSimilarityNetwork(self):
        pass

    


