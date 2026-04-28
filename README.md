## What is this repo?
A Python library for a wide variety of bioinformatics tasks written in an object-oriented manner.
The idea is to represent commonly encountered bioinformatics objects as classes and write functions within these classes to manipulate the objects and extract information from them. Currently the following entities are represented in classes:
- Files
  - FASTA
  - PDBx/mmCIF
- Databases
  - PDB
  - InterPro

Different file formats and database APIs can easily be added due to the modular architecture.
I'm also considering taking a more biology-oriented approach by for example representing protein structures as a class (as another layer of abstraction upon the PDBx/mmCIF format). Where possible, each class and its functions are tested using the python unittest framework. The test suite can be found under the test/ folder.

## Examples
### Downloading structures from PDB given an EC number
TODO

### Finding sequences that match a given pattern in a .FASTA file.
In this example I'm using the canonical catalytic serine found in a wide variety of enzymes (proteases, cutinases, carboxylic esterases...).
Usually the catalytic serine is found in close vicinity of a histidine and aspartic acid residue at the structural level.
The serine alcohol, kept in a negatively charged state due to the basic properties of the histidine residue, acts as a nucleophile capable of attacking electron-poor carbons (e.g., in an ester bond). At the amino acid sequence level, this catalytic serine is strongly conserved and can be identified by the pattern 'GxSxG', where x can be any amino acid.
Now let's say we have a FASTA file with some amino acid sequences, and we want to know which sequences match this pattern, and at which location in the string.
Let's use biolib to  write a small python script for this:

```
# Let's import the necessary FastaFile class first:
from biolib.classes.files.fasta import FastaFile
# Standard library import to work with file paths:
from pathlib import Path  

# Now we can load FASTA files by providing the path and call some methods on it:
path_to_fasta: Path = Path('path/to/fasta')  # Specifying the path to our fasta file.
my_fasta_file: FastaFile = FastaFile(path_to_fasta)  # Initiating the FastaFile object.

# We can now use the matchPattern() function to find our catalytic serine.
pattern: str = r'G.S.G'  # Here I'm defining the pattern in regular expression language
# Now we call the relevant function on our FastaFile object:
my_fasta_file.matchPattern(pattern)  
```

And voila, that's it! If any matches are found, the output will be a list of tuples formatted as follows:
```
[('sequence_header', start_position, end_position, 'matched_pattern')]
```
Where 'sequence_header' is what comes after the '>' in your fasta file for that specific sequence, and 'matched_pattern' is the actual string that was matched.
You can put this code in a main() function and call it from your terminal (don't forget `if __name__ == "__main__: main()"`!).
Alternatively, under the scripts/ folder you can find a functional script 'find_catalytic_serine_in_fasta.py' that you can run directly.
You can simply provide the path to your fasta file and it will print the output to your terminal.
Note that for this specific example you would still have to confirm that this is indeed the desired catalytic serine by loading the structure in Pymol/ChimeraX... and selecting the correct amino acid position.
Then you can check if any His or Asp residues are close by.
On a test run containing sequences of multicopper oxidases (presumed not to have catalytic serines), I still found 143/538 sequences matching the pattern.
This method also does not find overlapping matches, something I might implement in the future.
Lastly, notice that one can essentially define any pattern, given some basic regex knowledge, and use this function to find which sequences contain the pattern.

### Extracting all polypeptide sequences from a collection of .CIF files.
TODO