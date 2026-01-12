! /home/u0173836/miniconda3/envs/basic_data_science/bin/python3

from pathlib import Path
import pandas as pd
import sys

class FastaFile:
    """
    This class represents a FASTA file and contains method to read them.
    """

    def __init__(path_to_fasta: Path):
        
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

        ### DEFENSIVE CHECKS:
        if not path_to_fasta.exists():
            raise FileNotFoundError(f"{path_to_fasta} does not exist.")
        
        if not path_to_fasta.is_file():
            raise FileNotFoundError(f"{path_to_fasta} is not a file.")

        if not isFasta(path_to_fasta):
            raise ValueError(f"The given file is not a fasta file")


        ### INIT:
        self.path_to_fasta = path_to_fasta


    def toDict(self) -> dict:
        """
        Read this FASTA file and return as a dictionary of header:sequence
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


    def toCsv(self, path_to_out: str) -> None:
        """
        Converts this FASTA file to csv file.
        """
        # Convert output file to Path object:
        path_to_out = Path(path_to_out)
       
       # Make the directories if they do not exist already:
        path_to_out.parent.mkdir(parents=True, exist_ok=True)

        # Create dataframe from fasta file:
        df_fasta = self.fastaToDf()
       
       # Write to csv:
        df_fasta.to_csv(path_to_out)
        
        return None
    
