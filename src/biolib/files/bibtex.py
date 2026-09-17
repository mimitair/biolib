from pathlib import Path

import bibtexparser as btp
import networkx as nx
import nameparser as np

class Bibtex:
    def __init__(self, path_to_bibtex: str|Path):
        path: Path = Path(path_to_bibtex).resolve()

        self.full_path: Path = path
        self.library: btp.Library = btp.parse_file(path)

        return None

    
    @staticmethod
    def standardizeName(name: str) -> str:
        parsed_name: ParsedName = np.parse(name.strip())
        if parsed_name.given:
            return parsed_name.family.lower() + ' ' + parsed_name.given.lower()[0]
        elif parsed_name.suffix:
            return parsed_name.family.lower() + ' ' + parsed_name.suffix.lower()[0]
        elif parsed_name.title:
            return parsed_name.family.lower() + ' ' + parsed_name.title.lower()[0]
        else:
            return parsed_name.family.lower()

    def standardizeNames(self, names: list) -> list:
        return [self.standardizeName(name) for name in names]

    
    def getAmountOfEntries(self) -> int:
        pass

    
    def getAmountOfAuthorsPerEntry(self) -> dict:
        pass

    
    def getAuthors(self) -> list:
        """
        """
        result: list = []
        
        for entry in self.library.entries:
            author_names: list = self.standardizeNames(entry.get('author').value.split(' and '))
            result.extend(author_names)
            
        return result

    
    def getAuthorPairs(self) -> set:
        """
        """
        from itertools import combinations

        result: list = []
        for entry in self.library.entries:
            author_names: list = self.standardizeNames(entry.get('author').value.split(' and '))
            author_combinations: list[tuple] = list(combinations(author_names, 2))
            result.extend(author_combinations)

        return set([tuple(sorted(m)) for m in result])  # sort each pair alphabetically and remove duplicates by converting to set


    def countAuthorPairOccurences():
        pass


    def getCoAuthorNetwork(self):
        g = nx.Graph()
        g.add_nodes_from(self.getAuthors())
        g.add_edges_from(self.getAuthorPairs())

        return g
        
