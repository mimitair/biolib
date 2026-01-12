from pathlib import Path

class PDB:
    """
    """
    def __init__(self):
        pass
    
    @staticmethod
    def getAccessionsFromEC(ec_number: str):
        """
        Returns a list of PDB accessions annotated with the given EC number
        """
        # Import the rcsbapi module:
        from rcsbapi.search import AttributeQuery
        
        print(f"Querying PDB for EC number {ec_number}")

        # The attribute for EC numbers is the following (see  https://search.rcsb.org/structure-search-attributes.html)
        attribute = "rcsb_polymer_entity.rcsb_ec_lineage.id" 
        
        # Building the query:
        query = AttributeQuery(
                attribute = attribute,
                operator = "exact_match",
                value = ec_number
                )

        # Run the query as a function to obtain the results:
        results = list(query())
    
        print(f"Got {len(results)} PDB accessions as a result.")

        return results

    @staticmethod
    def downloadCIFFromEC(ec_number: str, out_dir: Path):
        """
        Downloads all PDB entries associatex with the given EC number and dumps them in out_dir
        """
        # Get the list of accession IDs for this EC number:
        accessions_list: list = PDB.getAccessionsFromEC(ec_number)

        # Import the ModelQuery class from the rcsbapi module:
        from rcsbapi.model import ModelQuery
        
        # Initiate a query:
        query: ModelQuery = ModelQuery()
        
        print(f"Downloading CIF files from PDB for EC {ec_number}. files will be stored in {out_dir}.")
        # Execute the query for multiple structures
        # We could also use the get_multiple_structures() function, but this returns a dictionary and does not work for writing to disk for some reason.
        # Hence, we just loop over the list and use the get_full_structure() function:
    
        query.get_multiple_structures(entry_ids = ["1AGY", "1A2A"],
                                      query_type = "full",
                                      encoding = "cif",
                                      download = True,
                                      compress_gzip = False,
                                      file_directory = out_dir)
        
        return None


def main():
    out_dir = Path("/home/u0173836/playground/downloads")

    PDB.downloadCifFromEC("3.1.1.1", out_dir)
    


if __name__ == "__main__":
    main()

