# Default library imports:
import pandas as pd
import sys
import re
from pathlib import Path

class PDBCIFFile:
    """
    This class represents a PDBx/mmCIF file as described in https://mmcif.wwpdb.org/, and contains methods to parse and manipulate them.
    """
    def __init__(self, path_to_pdbcif: Path):
        ### DEFENSIVE CHECKS ###
        if not isinstance(path_to_pdbcif, Path):
            raise TypeError("Path to PDB/mmCIF file must be a Path object.")
        
        if not path_to_pdbcif.exists():
            print(f"{path_to_pdbcif} does not exist.")
            sys.exit()
            
        if not path_to_pdbcif.is_file():
            print(f"{path_to_pdbcif} is not a file.")
            sys.exit()
        
        ### INIT ###
        self.file_path: Path = path_to_pdbcif  # Path to the PDBx/mmCIF
        self.file_name: str = path_to_pdbcif.name  # Name of the file
        self.file_contents: str = path_to_pdbcif.read_text()  # Content of the file as a string

    def matchCategory(self, category: str) -> str | None:
        """
        Matches any non-loop category and returns the data block as a string

        Input:
            - category: str: the desired category to extract

        Returns:
            - str: the data block as a string
            - None: if no match is found
        """
        # Define the pattern:
        pattern = rf"#\s*\n_{category}[^#]+#"

        # Search in the file contents:
        m = re.search(pattern, self.file_contents)

        # return the full match if one was found:
        if m is not None:
            return m.group()
        else:
            return None
        
    def categoryToDict(self, category: str) -> dict | None:
        """
        Extract all key/value pairs belonging to a given CIF category
        (non-loop blocks only).

        Input:
            - category: str: Name of the category to extract
        
        Returns:
            - dict: {field_name : value_string}.
            - None: If the category could not be matched.
        """
        # Pattern for CIF semicolon-delimited multiline values
        semicolon_value = r"""
            ;                       # opening semicolon
            ([\s\S]*?)              # non-greedy capture of all characters
            ;                       # closing semicolon on its own line
        """

        # Pattern for normal (single-line) values
        single_value = r"(\S+)"

        # Full regex:
        # 1. Match keys like: _entity_poly.xxx
        # 2. Match either:
        #       - a semicolon block
        #       - or a single token
        pattern = re.compile(
            rf"""
            _{category}\.(\S+)              # group(1): the field name after the dot
            \s*                             # whitespace
            (?:                             # start value group
                {semicolon_value}           # group(2): multiline value OR
                |                           # OR
                {single_value}              # group(3): single-line token
            )
            #
            """,
            re.VERBOSE
        )



        # First match the data block of the category as a string:
        data_block: str = self.matchCategory(category)

        # If the category was found:
        if data_block is not None:
            # Inititate result dictionary:
            result = {}
            
            # Iterate over the lines that matched the pattern:
            for match in pattern.finditer(data_block, re.MULTILINE):
                # Extract the field name:
                field = match.group(1)
                # Extract the value, which is either a semicolon multiline value, or a single value
                if match.group(2) is not None:
                    # Multiline value: strip trailing newline and newlines inside the match
                    value = match.group(2).strip().replace("\n", "")
                else:
                    value = match.group(3).strip()

                # Add the result in the dictionary:
                result[field] = value

            # Return the result dictionary:
            return result

        # If no match was found, return None:
        else:
            return None


    def matchLoopCategory(self, category: str) -> str | None:
        """
        Function to match any loop block category using regular expressions.
        Can also be used to check if the given category is present in the file.
        Returns None if it isn't.

        The following pattern is matched:

            #
            loop_
            _<category>
            <Anything that is not a hashtag>
            #

        Input:
            - category: str: The loop category to parse

        Returns
            - str: The loop block from its starting "#" till its ending "#" as a string
            - None: If no match is found
        """

        pattern = rf"#\s*\nloop_\s*\n_{category}[^#]+#"
        
        # Search for the pattern:
        match = re.search(pattern, self.file_contents)
        
        # If there is no match, we need to inform the user:
        if match is None:
            return None
        
        else:
            # We can now return the entire matched string by calling match.group():
            return match.group()


    def loopCategoryToDf(self, category: str) -> pd.DataFrame | None: 
        """
        Convert any loop block, given its category name, to a pandas dataframe.

        Input:
            - category: str: The category to extract (f.e.: "atom_site")

        Returns:
            - pd.DataFrame: The loop block in a dataframe. Beware that no further processing is done. Every value is essentially a string.
            - None: If the category cannot be found.
        """

        # Initiate list for column names and data:
        columns: list = []
        data: list = []

        # First we extract the loop block as a string using the matchLoopCategory() function:
        loop_block: str = self.matchLoopCategory(category)

        # If the category cannot be matched, return None
        if loop_block is None:
            return None

        # Split the string by newlines:
        loop_block: list[str] = loop_block.split("\n")

        # Loop over the loop block, which is now a list of strings where every string is a line of text in the loop block:
        for line in loop_block:
            if line.startswith("#") or line.startswith("loop_"):
                # These are the 'paddings' and can be omitted
                continue
            
            elif line.startswith(f"_{category}."):
                # These contain the column names. We need the part after the dot:
                columns.append(line.split(".")[-1].strip())
                continue
            
            else:
                # What remains are the data lines. We should split them by whitespace, strip them, and add them to the data list
                # PROBLEM: Some values are enclosed within ' ' and contain spaces within them. Hence, split does not work!
                # FIX: Use regex
                # PROBLEM: Some carbon atoms are labeled as "C1'", so this regex pattern does not work because it will start matching everything between that apostroph and the next.
                # FIX: Simply put a space before the first apostophe of the second pattern. This prevents C1' <multiple values> C1' from being matched, but it will match <whitespace>'<some value>' 
                pattern: str = r"[^\s']+| '[^']*'"  # AI generated, matches anything that is not a whitespace or ' OR something enclosed within '' and preceded by a space.
                
                # findall will return all matches as list:
                data_line: list = re.findall(pattern, line)
                
                # Append to the data list
                data.append(data_line)

        # Now we can convert the columns and data list into a dataframe and return it:
        return pd.DataFrame(columns=columns, data=data)

    def getHeteroAtoms(self) -> set:
        """
        Returns a set of hetero atoms (labeled as HETATM) in the atom_site category

        Returns:
            - set: A set of HETATM names in this PDBx/mmCIF file. Empty set if nothing is found

        """
        # Convert the atom_site category to a dataframe:
        df_atoms = self.loopCategoryToDf('atom_site')

        # Return the 'label_comp_id' column for each row where 'group_PDB' == 'HETATM' as a set
        return set(df_atoms[df_atoms['group_PDB'] == 'HETATM']['label_comp_id'].to_list())


    def getAminoAcidSequences(self) -> list:
        """
        Returns the amino acid sequences of each polymer entity in this file
        """

        # Extract the entity_poly category containing the sequence information:
        sequence_data_block: dict = self.categoryToDict('entity_poly')

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
    
    
    def atomNumberToPoint(self, atom_number: int) -> Point:
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
   
class PDBCIFFileCollection():
    """
    Class that represents a collection of PDBx/mmCIF files and methods to manipulate them.
    """

    def __init__(self, path_to_folder: Path):
        ### DEFENSIVE CHECKS:
        if not isinstance(path_to_folder, Path):
          pass

        if not path_to_folder.is_dir():
            pass

        if not path_to_folder.exists():
            pass

        self.pdbcif_files: list = [PDBCIFFile(child) for child in path_to_folder.iterdir()]


        return None

    def writeSequencesToFasta(self, out_file: Path) -> None:
        """
        Writes all the amino acid sequences to a fasta file.
        Header of the faste file is the name of the file.
        If multiple amino acid sequences are in the cif file, they will be written as:
            > <file_name><number>
        """
        
        for pdbcif_file in self.pdbcif_files:
            sequence_data_block: dict = pdbcif_file.categoryToDict('entity_poly')
            if sequence_data_block is not None:
                pass
                
            else:
                print(pdbcif_file.file_name)
        
        
            
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

    


