from time import time
from utils import *


# A class GraphMolecule representing a Molecule in a form of graph


class GraphMolecule:
    def __init__(self, mol):
        self.mol = AutoMol(mol) # convert auto into rdkit mol




def main():
    iptacopan_smiles = "CCOC1CCN(C(C1)C2=CC=C(C=C2)C(=O)O)CC3=C(C=C(C4=C3C=CN4)C)OC"
    ipta = Chem.MolFromSmiles(iptacopan_smiles)

    Graphmol = GraphMolecule(iptacopan_smiles)
    print(Graphmol)

if __name__ == '__main__':
    main()