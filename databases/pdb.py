from pathlib import Path

class PDB:
    """
    Class to interact with the PDB database.
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
        
        print(f"Querying PDB for EC number {ec_number}.")

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
    
        print(f"Got {len(results)} PDB accessions.")

        return results


    @staticmethod
    def downloadCIFFromAccessions(accessions: list, out_dir: Path) -> None:
        """
        Downloads the given list of PDB accessions in out_dir
        """
        # Import the ModelQuery class from the rcsbapi module:
        from rcsbapi.model import ModelQuery
        
        # Create out-dir if it does not exist yet:
        out_dir.mkdir(exist_ok=True, parents=True)

        # Initiate a query:
        query: ModelQuery = ModelQuery()
        
        print(f"Downloading {len(accessions)} CIF files.")
        print(f"Files will be stored in {out_dir}.")
        # Execute the query for multiple structures
        query.get_multiple_structures(entry_ids = accessions,
                                      query_type = "full",
                                      encoding = "cif",
                                      download = True,
                                      compress_gzip = False,
                                      file_directory = out_dir)

        print("Downloading finished.")

        return None        
 
def main():
    out_dir = Path("/home/u0173836/playground/downloads")

    PDB.downloadCifFromEC("3.1.1.1", out_dir)
    


if __name__ == "__main__":
    main()

