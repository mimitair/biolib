from rdkit import Chem
import biotite.structure as struc
import biotite.structure.io.mol as molio

class SDFFile():
    def __init__(self, path_to_sdf_file: str|Path):
        pass

    @staticmethod
    def createFromSmiles(smiles: str, out_file: str|Path) -> Self:
        mol = Chem.MolFromSmiles(smiles)
        with Chem.SDWriter(out_file) as writer:
            writer.write(mol)

        return None

    def prepareForDocking(self) -> struc.AtomArray:
        pass
