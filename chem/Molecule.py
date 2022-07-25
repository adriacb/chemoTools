from chem.utils import *

class Molecule:
    def __init__(self, molecule: Chem.rdchem.Mol):
        self.molecule = molecule
        self.name = molecule.GetProp("_Name")
        self.properties = self.calcProp()
    

    def calcProp(self):
        """
        Calculate properties of the molecule
        """
        self.logP = Descriptors.MolLogP(self.molecule)
        self.MWt = Descriptors.MolWt(self.molecule)
        self.nHBD = Descriptors.NumHDonors(self.molecule)
        self.NumAromaticRings = Chem.Lipinski.NumAromaticRings(self.molecule)
        self.PFI = self.logP + self.NumAromaticRings
        self.properties = np.asarray([self.logP, self.MWt, self.nHBD, self.NumAromaticRings, self.PFI])
        
        return self.properties
    
    
        