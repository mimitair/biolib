import biotite.structure as struc
import biotite.structure.io.pdbx as pdbxio
from pathlib import Path
import numpy as np
import pandas as pd
import re
from biolib.files.pdbcif import CifFile

class FpocketOut():
    
    def __init__(self, path_to_fpocket_out: str | Path):
        ### DEFENSIVE CHECKS ###
        path: Path = Path(path_to_fpocket_out).resolve()
        if not path.exists():
            raise FileNotFoundError(f"{path} does not exist")
        if not path.is_dir():
            raise FileNotFoundError(f"{path} is not a directory")

        ### INIT ###
        self.full_path: Path = path
        self.info_path: Path = list(path.glob('*_info.txt'))[0]
        self.pdbcif_path: Path = list(path.glob('*_out.cif'))[0]


    def getID(self) -> str:
        return self.pdbcif_path.split('_')[0].upper()

    
    def countPockets(self) -> int:
        with self.info_path.open('r') as f:
            pass
    

    def getPocketCentroids(self) -> dict:
        """Returns the centroid coordinates for each pocket as a dictionary {pocket_id : centroid_coord}
        """
        atom_array: struc.AtomArray = CifFile(self.pdbcif_path).toBiotiteAtomArray()
        pockets: struc.AtomArray = atom_array[atom_array.res_name=='STP']
        centroids: dict = {}
        for pocket_id in np.unique(pockets.res_id):
            centroids[pocket_id] = struc.centroid(pockets[pockets.res_id==pocket_id])
        
        return centroids


    def pocketInfoToDf(self) -> pd.DataFrame:
        result: dict = {}
        with self.info_path.open('r') as f:
            pocket_id: int = 0
            data: list = []
            for line in f:
                if line.startswith('Pocket'):
                    pocket_id += 1
                    continue
                elif bool(line.strip())==False: # empty line, represents the end of a pocket
                    result[pocket_id] = data # add data to result dictionary
                    data = [] # clear data for the next pocket
                else:
                    data.append(line.split(':')[-1].strip()) # append value to data for the current pocket

        return pd.DataFrame.from_dict(result, orient='index', columns=['score', 'druggability_score', 'alpha_sphere_count', 'total_sasa', 'polar_sasa', 'apolar_sasa',
                                                                       'volume', 'mean_local_hydrophobic_density', 'mean_alpha_sphere_radius', 'mean_alpha_sphere_solvent_access',
                                                                       'apolar_alpha_sphere_proportion', 'hydrophobicity_score', 'volume_score', 'polarity_score', 'charge_score',
                                                                       'proportion_polar_atoms', 'alpha_sphere_density', 'center_of_mass_alpha_sphere_max_dist', 'flexibility'])
                    
                    
    def getPocketInfoAt(self, pocket_id: int) -> dict:
        """Returns the pocket characteristics for a given pocket id
        """
        # This is probably not the most efficient:
        # Get all the pockets as a df first:
        df = self.pocketInfoToDf()
        # Return as dict (no iloc because it is zero-based, indexing starts at 1 as the pockets in the original output file)
        return df.loc[pocket_id].to_dict()
    
    def getAminoAcidsNearPocket(self, pocket_id: int, zone: float) -> tuple:
        pass

        
class FpocketOutCollection():
    def __init__(self, path_to_fpocket_out_collection: str):
        self.full_path: Path = Path(path_to_fpocket_out_collection).resolve()

        return None

    def iterDirs(self):
        """Iterates over every fpocket output folder in this collection
        Yields FpocketOut objects
        """
        for child in self.full_path.iterdir():
            yield FpocketOut(child)

            
    def pocketInfoToDf(self, out_file: str):
        """Transform all the pocket info in each fpocket_out directory to one big csv file.
        IDs of the structures are added as a column
        """
        # Initiate empty dataframe
        df: pd.DataFrame = pd.DataFrame()
        # Iterate over all the subdirectories
        for fpocket_out in self.iterDirs():
            # Extract pocket info from the current fpocket_out
            df_info: pd.DataFrame = fpocket_out.pocketInfoToDf()
            # Add the id
            df_info['id'] = fpocket_out.getID()
            # Concatenate to master df
            df = pd.concat([df,df_info])

        # Write the result to the specified output file
        df.to_csv(out_file)


    def getFpocketOutFromID(self, id: str) -> Path:
        return list( self.full_path.glob(id+'*', case_sensitive=False))[0]


    
    
