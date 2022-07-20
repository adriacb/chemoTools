
from chem.utils import *


class Molecules:
    def __init__(self, molecules, setname=None):
        self.setName = str(uuid.uuid1()) if setname is None else setname
        self.molecules = molecules
        self.num_molecules = len(molecules)
        self.properties = None
        
        if type(molecules) == Chem.rdmolfiles.SDMolSupplier or all(isinstance(n, Chem.rdchem.Mol) for n in molecules):
            self.properties = self._calculate_properties()
        else:
            self.molecules = [f"Molecule_{i+1}" for i in range(self.num_molecules)]
            self.properties = molecules
        
        # It will store the results of the Dimensionality Reduction
        self.coordinates = None

    
    def _calculate_properties(self):
        
        properties = list()

        for mol in self.molecules:
            logp = Descriptors.MolLogP(mol)
            mwt = Descriptors.MolWt(mol)
            nhbd = Descriptors.NumHDonors(mol)
            nar = Chem.Lipinski.NumAromaticRings(mol)
            pfi = logp + nar
            ps = np.asarray([logp, mwt, nhbd, nar, pfi])

            mol.SetProp("MWt", str(mwt))
            mol.SetProp("logP", str(logp))
            mol.SetProp("nHBD", str(nhbd))
            mol.SetProp("NumAromaticRings", str(nar))
            mol.SetProp("PFI", str(pfi))
            
            if properties is None:
                properties = ps
                print(properties)
            else:
                properties.append(ps)

        return np.array([x for x in properties])
    
    def to_df(self):
        """
        Convert the Molecules object into a Pandas DataFrame
        ----------
        Returns:
            Pandas DataFrame
        """
        columns = ["ID", "Molecule", "logP", "MWt", "nHBD", "NumAromaticRings", "PFI"]
        # Add the coordinates if they are available
        if self.coordinates is not None:
            for i in range(self.coordinates.shape[1]):
                columns.append("PC" + str(i))

        df = pd.DataFrame(columns=columns)
        # Add the properties
        df["ID"] = [str(i) for i in range(self.num_molecules)]
        df["Molecule"] = self.molecules
        df["logP"] = self.properties[:, 0]
        df["MWt"] = self.properties[:, 1]
        df["nHBD"] = self.properties[:, 2]
        df["NumAromaticRings"] = self.properties[:, 3]
        df["PFI"] = self.properties[:, 4]

        # Add the coordinates if they are available
        if self.coordinates is not None:
            for i in range(self.coordinates.shape[1]):
                df["PC" + str(i)] = self.coordinates[:, i]

        return df
    

    def to_sdf(self, path=None):
        """
        Convert the Molecules to sdf file.
        ----------
        path: str, optional
            The path to save the sdf file. If not provided, the file will be saved in the current directory.
        """
        if path is None:
            path = os.getcwd()
        
        Chem.PandasTools.WriteSDF(self.to_df(), path + "/" + self.setName + ".sdf")
        print("SDF file saved in {}/{}.sdf".format(path, self.setName))
        
        


    def __add__(self):
        """
        TO DO:
        - Add two or more sets of molecules and transform them into a ChemSpace object
        """
        pass

    def __str__(self):
        return "<Molecules> object with {} molecules.\nMolecules set name {}\nProperties dimension: {}".format(len(self.molecules), self.setName, self.properties.shape)


