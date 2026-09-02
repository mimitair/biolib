from pathlib import Path
import sys
import pandas as pd
import re
import matplotlib.pyplot as plt

class FastaFile:
    """
    This class represents a FASTA file, and contains methods to parse and manipulate them.
    """

    def __init__(self, path_to_fasta: Path):
        
        def isFasta(path_to_fasta: Path) -> bool:
            """
            Helper function. Checks if the suffix is one of the accepeted fasta suffices.
            Returns True if it's all good.
            """
            # define set of accepted suffices:
            accepted_suffices: set = {".fa", ".fasta", ".fas", ".fna", ".faa"}
            
            if path_to_fasta.suffix in accepted_suffices:
                return True
            
            return False

        ### DEFENSIVE CHECKS ###
        if not path_to_fasta.exists():
            raise FileNotFoundError(f"{path_to_fasta} does not exist.")
        
        if not path_to_fasta.is_file():
            raise FileNotFoundError(f"{path_to_fasta} is not a file.")

        if not isFasta(path_to_fasta):
            raise ValueError(f"The given file is not a fasta file")

        ### INIT ###
        self.path_to_fasta: Path = path_to_fasta

    @property
    def count(self) -> int:
        """
        Returns the amount of sequences in this fasta file.
        """
        return len(self.toDict())

    def toDict(self) -> dict:
        """
        Read this FASTA file and return as a dictionary of {header:sequence}
        """
        # inititate empty dictionary to store the result
        result: dict = {}
        
        # open the file
        with self.path_to_fasta.open('r') as f:
            # loop over the lines
            for line in f:
                line = line.strip()  # remove leading and trailing whitespace
                # if line starts with '>', add it as a key to result dictionary
                if line.startswith('>'):
                    sequence_id = line[1:].strip() # omit the '>' character and strip again
                    # initiate empty string to store the coming sequence
                    result[sequence_id] = ""
                    continue
                # if not, append the string to the current sequence ID:
                result[sequence_id] += line
            
        return result

    
    def getHeaders(self) -> list:
        return list(self.toDict().keys())

    
    def getSequences(self) -> list:
        return list(self.toDict().values())

    
    def toDf(self) -> pd.DataFrame:
        """
        Read this FASTA file and return as a dataframe with column 'header' and 'sequence'
        """
        # convert to dictionary:
        fasta_dict: dict = self.fastaToDict()
        
        # convert dictionary to pandas dataframe with keys as index for rows:
        df_fasta: pd.DataFrame = pd.DataFrame.from_dict(fasta_dict, orient='index', columns=['Sequence'])

        # name index col:
        df_fasta.index.name = "Header"
        
        return df_fasta


    def toCsv(self, out_file: str) -> None:
        """
        Converts this FASTA file to csv file.
        """
        # Convert output file to Path object:
        out_file: Path = Path(out_file)
       
       # Make the directories if they do not exist already:
        out_file.parent.mkdir(parents=True, exist_ok=True)

        # Create dataframe from fasta file:
        df_fasta = self.fastaToDf()
       
       # Write to csv:
        df_fasta.to_csv(path_to_out)
        
        return None

    @staticmethod
    def csvToFasta(in_file: Path, out_file: Path) -> int:
        """
        Converts a csv file with 'id' and 'aa_seq' columns to a fasta file

        Input:
           - in_file: Path: Input CSV file
           - out_file: Path: Output FASTA file
        """
        df: pd.DataFrame = pd.read_csv(in_file)

        d: dict = df.set_index('id')['aa_seq'].to_dict()
        
        with out_file.open("w") as f:
            for key, value in d.items():
                f.write(">" + str(key) + "\n" + str(value) + "\n")

        return 0

    def splitToSeparateFiles(self, out_dir: Path) -> int:
        # Make the output dir if it does not exist
        out_dir.mkdir(parents=True, exist_ok=True)

        # Convert fasta to dict:
        d: dict = self.toDict()

        # Loop over dict items:
        for header,sequence in d.items():
            header = str(header)
            file_name: str = header + '.fasta'
            out_file: Path = out_dir / file_name
            with out_file.open("w") as f:
                f.write(">" + header + "\n" + str(sequence))

        return 0
        
    def filterByLength(self, min_length: int, max_length: int, out_file: Path = None):
        """
        Filters out all sequences below min_length or above max_length.
        Writes cleaned fasta file to out_file
        """        
        result: dict = {}

        count: int = 0
        for header,sequence in self.toDict().items():
            if len(sequence) < min_length or len(sequence) > max_length:
                count += 1
            else:
                result[header] = sequence

        print(f"Removed {count} sequences according to length thresholds. Writing filtered fasta file to {out_file}")

        if out_file is not None:
            with out_file.open("w") as f:
                for header, seq in result.items():
                    f.write(">" + str(header) + "\n" + str(seq) + "\n")
        
        return result
                
    def matchPattern(self, pattern: str) -> list[(str, int, int, str)]:
        """
        Returns a list of all entries in this fasta file that match the given pattern
        Each element of the list is a tuple formatted as: (entry header, start position, end position, matched pattern)
        All non-overlapping matches are returned, meaning that a sequence can have multiple matches.
        Returns an empty list if no match is found
        TODO Detect overlapping matches as well?
        """
        result = []
        pattern = re.compile(pattern)
        
        for header, sequence in self.toDict().items():
            matches = pattern.finditer(sequence)
            for m in matches:
                result.append((header, m.start(), m.end(), m.group()))

        return result

    def getLengths(self) -> list:
        return [len(sequence) for sequence in self.toDict().values()]

    
    def plotLengthHistogram(self, xlabel: str, ylabel: str, title: str) -> None:
       data: list = self.getLengths()

       fig, ax = plt.subplots()

       ax.hist(data, bins=50)

       plt.show()

       return None

    def plotLengthBox(self):
        data: tuple = self.getLengths()

        fig, ax = plt.subplots()

        ax.boxplot(data)

        plt.show()

        return None

    
class FastaFileCollection:
    def __init__(self):
        pass
