import biotite.structure as struc
import biotite.structure.io.pdbx as pdbxio
from pathlib import Path
import numpy as np
import pandas as pd
import re

class FpocketOut():
    
    def __init__(self, path_to_fpocket_out: str):
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
        
    def countPockets(self) -> int:
        with self.info_path.open('r') as f:
            pass
    
    def toBiotiteAtomArray(self) -> struc.AtomArray:
        """Convert to biotite atom array object
        """
        mycif: pdbxio.CIFFile = pdbxio.CIFFile.read(self.pdbcif_path)
        mycif.block['atom_site']['pdbx_PDB_model_num'] = np.ones(len(mycif.block['atom_site']['id'].as_array()), dtype=np.int32)

        atom_array: struc.AtomArray = pdbxio.get_structure(mycif, include_bonds=True, extra_fields=['atom_id', 'charge'])

        return atom_array

        
    def getPocketCentroids(self) -> dict:
        """Returns the centroid coordinates for each pocket in a dictionary {pocket_id : centroid}
        """
        atom_array: struc.AtomArray = self.toBiotiteAtomArray()
        pockets: struc.AtomArray = atom_array[atom_array.hetero==True]
        centroids: dict = {}
        for pocket_id in np.unique(pockets.res_id):
            centroids[pocket_id] = struc.centroid(pockets[pockets.res_id==pocket_id])
        
        return centroids

    def getNearestPocket(self, atom: struc.Atom) -> tuple:
        return 0

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
                    
                    
    def findTriads(self, *args, **kwargs) -> struc.AtomArray:
        return PdbCifFile(self.pdbcif_path).findTriads(args, kwargs)


class FpocketOutCollection():
    pass
    
