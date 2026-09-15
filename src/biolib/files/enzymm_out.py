import pandas as pd
from pathlib import Path

class EnzymmOut:
    def __init__(self, path_to_enzymm_out: str|Path):
        ### DEFENSIVE CHECKS ###
        path: Path = Path(path_to_enzymm_out)

        if not path.exists():
            raise FileNotFoundError

        if not path.is_file():
            raise FileNotFoundError
        
        ### INIT ###
        if path.suffix == '.parquet':
            self.df = pd.read_parquet(path_to_enzymm_out)
            
        elif path.suffix == '.tsv':
            self.df = pd.read_csv(path, sep='\t', comment='#')
        
    def filterByRMSD(self, out_file: str|Path) -> pd.DataFrame:
        """Filters by lowest rmsd for each query id and retains only relevant columns
        """
        df = pd.DataFrame()
        for name, group in self.df.groupby('query_id'):
            to_retain: pd.Series = group.nsmallest(1, columns='rmsd')
            df = pd.concat([df, to_retain])

        df.to_csv(out_file, columns=['query_id', 'template_pdb_id', 'template_mcsa_id', 'template_ec', 'matched_residues', 'rmsd', 'query_residue_count'], sep='\t', index=False)

        return df[['query_id', 'template_pdb_id', 'template_mcsa_id', 'template_ec', 'matched_residues', 'rmsd', 'query_residue_count']]

    def getQueryChain(self, query_id: str) -> str:
        return self.df.loc[self.df['query_id']==query_id]['matched_residues'].to_string().split(',')[0].split('_')[1]
