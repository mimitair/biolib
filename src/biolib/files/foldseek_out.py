from pathlib import Path
import pandas as pd

class EasyClustOut():
    def __init__(self, path_to_easyclust_out: str):
        path: Path = Path(path_to_easyclust_out).resolve()

        ### TODO DEFENSIVE CHECKS

        self.full_path: Path = path
        self.adjacency_path: Path = list(path.glob('*_cluster.tsv'))[0]
        self.all_seqs_path: Path = list(path.glob('*_all_seqs.fasta'))[0]
        self.rep_seqs_path: Path = list(path.glob('*_rep_seq.fasta'))[0]

        self.df = pd.read_csv(self.adjacency_path, sep='\t', names=['representative', 'node'])

    def getClusterCount(self) -> int:
        return self.df['representative'].nunique()

    
    def getAverageClusterSize(self) -> int:
        pass

    
    def getMaxClusterSize(self) -> int:
        pass

    
    def getMinClusterSize(self) -> int:
        pass

    
    def getSingletonCount(self) -> int:
        df_grouped = self.df.groupby('representative').size()
        return sum(df_grouped==1)

    
    def getRepresentatives(self) -> list:
        """Returns the representatives of each cluster as a list
        """
        return self.df['representative'].unique()

    def summarize():
        pass
    
    def filterRepresentativesByLength(self, min_length: int, max_length: int, out_dir: str = None) -> pd.DataFrame:
        """
        Mutates
           - self.df
        """
        from biolib.files.fasta import FastaFile

        rep_fasta: FastaFile = FastaFile(self.rep_seqs_path)
        filtered_headers = list(rep_fasta.filterByLength(min_length, max_length).keys())

        self.df = self.df.loc[self.df['representative'].isin(filtered_headers)]

        if out_dir is not None:
            out_file: Path = Path(self.adjacency_path.stem + f'_filtered_length{min_length}-{max_length}.tsv')
            path_to_filtered_cluster_out: Path = Path(out_dir) / out_file
            self.df.to_csv(path_to_filtered_cluster_out, header=False, index=False)
            
        return self.df
        

    
