from pathlib import Path
import pandas as pd
import sys
from pathlib import Path
import re

from ..util import misc, point

class PDBCIFFile:
    """
    This class represents a PDBx/mmCIF file as described in https://mmcif.wwpdb.org/
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
        self.path_to_pdbcif: Path = path_to_pdbcif  # Path to the file
        self.df_atoms: pd.DataFrame = None  # Atom information in a dataframe. Empty by default.
        self.df_site: pd.DataFrame = None
    
    
    def matchLoopCategory(self, category: str) -> str:
        """
        Function to match any loop block using regular expressions.

        Input:
            - category: The loop category to parse

        Returns
            - str: the loop block from its starting "#" till its ending "#" as a string
    
        """
       
        # First we construct a regular expression pattern to capture a loop block
        # Let's break it down:
        # Group 1 (start of the loop block): "(^#\s*\n^loop_\s*\n)"
        #       "^" matches the start of the string
        #       "#" is the hashtag that ends and starts a new data block
        #       "\s*\n" means any whitespace character zero or more times, and ending with an explicit newline.
        #       "^loop_\s*\n" matches the start of the loop block followed by zero or more whitespaces and ending with a newline
        #
        # Group 2 (column names): "(^_" + re.escape(category) + r"\.[a-zA-Z_]+\s*\n)+"
        #       "^_" matches the start of a category name
        #       "+ re.escape(category) here we briefly interrupt the string to insert the given category name. we use re.escape to avoid whitespace character being interpreted literally in the string.
        #       "\." matches a literal dot
        #       "[a-zA-Z_]+" Matches any alphabetical character and underscores one or more times. So far no attribute names contain numbers or other symbols to my knowledge
        #       "\s*\n)+" Again any whitespace or none can follow, and ending with a newline. The group is closed here and can be matched one or more times.
        #
        # Group 3 (the actual data lines): "(^[^#]+\s*\n)+"
        #       "^[^#]" the first ^ matches the start of the string. the part in square brackets: anything that is NOT a hashtag. First I used [^\n]+ to match anything that is not a newline, but then the expression matched the final hashtag as well, and the match ran till the end of the file
        #       "\s*\n" again any whitespace character zero or more times and ending with a newline
        #       ")+" End the group and allow it to match one or more  times
        # Finally "^#" defines the end of the loop block
        
        pattern = r"(^#\s*\n^loop_\s*\n)(^_" + re.escape(category) + r"\.[a-zA-Z_]+\s*\n)+(^[^#]+\s*\n)+^#"
 
        # read the contents of the file as a string.
        content = self.path_to_pdbcif.read_text()
 
        # Search for the pattern, MULTILINE for correct interpretation of the '^' character:
        match = re.search(pattern, content, re.MULTILINE)
        
        # If there is no match, we need to inform the user:
        if match is None:
            return f"Could not match {category}"

        # We can now return the entire loop block of choice:
        return match.group()


    def loopCategoryToDf(self, category: str) -> pd.DataFrame: 
        # Initiate list for column names and data:
        columns = []
        data = []

        # First we extract the loop block as a string using the matchLoopCategory() function:
        loop_block = self.matchLoopCategory(category)

        # Split the string by newlines:
        loop_block = loop_block.split("\n")

        # Loop over the loop block, which is now a list of strings where every string is a line of text in the loop block:
        for line in loop_block:
            if line.startswith("#") or line.startswith("loop_"):
                # These are the 'paddings' and can be omitted
                continue
            
            elif line.startswith(f"_{category}."):
                # These contain the column names:
                columns.append(line.split(".")[-1].strip())
                continue
            
            else:
                # What remains are the data lines. We should split them by whitespace, strip them, and add them to the data list
                #TODO Some values are within '' and contain spaces within them. Hence, split does not work!
                print(line.strip().split())
                data.append(line.strip().split())

        # Now we can convert the columns and data list into a dataframe and return it:
        return pd.DataFrame(columns=columns, data=data)

    
    def atomsToDf(self): 
        """
        Reads this PDBx/mmCIF file and extracts the atom information into a pandas dataframe.
        Identifies the '_atom_site' category and stores every attribute as a column name.
        Then identifies every line that starts with 'ATOM' and fills the dataframe with it.
        The atoms dataframe is then stored internally.
        """
        
        print(f"Extracting atoms from {self.path_to_pdbcif.name}.")
        
        atom_site_attributes: list = []  # list of strings
        atoms: list = []  # list of lists
        
        # Open the PDB file:
        with self.path_to_pdbcif.open(mode='r') as pdb_file:
            # Iterate over each line
            for line in pdb_file:
                # Extract the category attributes:
                if line.startswith("_atom_site."):  # The dot is necessary to avoid _atom_sites category being parsed
                    # Divide the line at the dot symbol
                    attribute = line.strip().split(".")[-1] # Part that comes after the dot
                    # Append to the attributes list
                    atom_site_attributes.append(attribute)
                    
                # If the line starts with "ATOM"
                elif line.startswith("ATOM"):
                        # Strip new lines and split the string at every whitespace:
                        line = line.strip().split()
                        # Now we append this list as a row to the atoms list:
                        atoms.append(line)
                
            # Now we pour this into a dataframe:
            df = pd.DataFrame(data=atoms, columns=atom_site_attributes)
            
            # Drop the redundant 'group_PDB' column and set the index to the id attribute:
            df = df.drop(columns=['group_PDB'])
            df = df.set_index('id')
            
            print(f"Extracted {df.shape[0]} atoms and {df.shape[1]} attributes.")
            
        # Store internally:
        self.df_atoms = df
            
        return None
    
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
        # Print a pretty ASCII art triad:
        print("Looking for the following triad(s):")
        misc.printTriad(res1, res2, res3, atom1, atom2, atom3, min_dist_1_2, min_dist_2_3, min_dist_3_1, max_dist_1_2, max_dist_2_3, max_dist_3_1)

        # Initiate empty lists to store distances between pairs of Points (atoms associated to residue numbers):
        distances_1_2: list = []
        distances_2_3: list = []
        distances_3_1: list = []
        
        # Store each given atom in a {residue_number : Point} dictionary
        res1_points: dict = self.residueAtomNamesToPoints(res1, atom1)
        res2_points: dict = self.residueAtomNamesToPoints(res2, atom2)
        res3_points: dict = self.residueAtomNamesToPoints(res3, atom3)
        
        # Now we have to calculate the distances between each pair of Point objects.
        # This requires some nested for loops
        # Iterate over the given atoms of residue 1:
        for key1, point1 in res1_points.items():  # key1=residue number, point1=Point object representing an atom of the aa residue
            # Iterate over the given atoms of residue 2:
            for key2, point2 in res2_points.items():
                distance_1_2 = point1.distance(point2)  # Distance between p1 and p2 (see Point class)
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
                    


def main():
    path_to_pdb = Path('/home/green/projects/bio/biolib/5XJH.cif')
        
    my_file = PDBCIFFile(path_to_pdb)

    print(my_file.matchLoopCategory('pdbx_audit_revision_group'))
    print(my_file.loopCategoryToDf('pdbx_audit_revision_group'))  

    #my_file.atomsToDf()
    
    #my_file.findTriads('SER', 'HIS',                        'ASP',                       'OG',                        'ND1',                        'OD1',                        1,                        1,                        6,                        5,                        4,                        7)
    
    
if __name__ == "__main__":
    main()

   
