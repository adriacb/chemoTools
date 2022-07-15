import pickle
import uuid
import time

import pandas as pd
import numpy as np
np.set_printoptions(suppress=True)
pd.set_option('display.float_format', lambda x: '%.5f' % x)

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, Descriptors

from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler

from scipy.spatial import distance

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


class Molecules:
    def __init__(self, molecules, setname=None):
        self.setName = str(uuid.uuid1()) if setname is None else setname
        self.molecules = molecules
        self.num_molecules = len(molecules)
        self.properties = None
        
        if type(molecules) == Chem.rdmolfiles.SDMolSupplier or type(molecules[0]) == Chem.rdchem.Mol:
            self.properties = self._calculate_properties()
        else:
            self.molecules = [f"Molecule_{i+1}" for i in range(len(molecules))]
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
            pfi = Descriptors.MolLogP(mol) + Chem.Lipinski.NumAromaticRings(mol)
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

    def __add__(self):
        """
        TO DO:
        - Add two or more sets of molecules and transform them into a ChemSpace object
        """
        pass

    def __str__(self):
        return "<Molecules> object with {} molecules.\nMolecules set name {}\nProperties dimension: {}".format(len(self.molecules), self.setName, self.properties.shape)






class ChemSpace:
    """
    Parameters
        chemspace - list of <Molecules> objects
    ===========================
    """

    def __init__(self, chemspace: list[Molecules]):

        assert type(chemspace) == list, "Chemspace must be a list of Molecules objects."
        assert len(chemspace) > 0, "Chemspace must contain at least one molecule."
        assert type(chemspace[0]) == Molecules, "Chemspace must be a list of Molecules objects."


        self.hash_space = {}
        self.index_space = []
        self.NumOfProperties = list()

        # Create a dictionary of Molecules objects, with the molecule name as the key.
        # Save the length of the molecules list as the value.
        for molecules in chemspace:
            self.hash_space[molecules.setName] = molecules
            self.index_space.append(molecules.properties.shape[0])
            self.NumOfProperties.append(molecules.properties.shape[1])

        assert len(set(self.NumOfProperties)) == 1, "All the sets of molecules must have the same number of properties."

        # Dimensionality reduction configuration
        self.nComponents = None

        self.models = {'PCA': PCA(n_components=self.nComponents),
                       'SVD': TruncatedSVD(n_components=self.nComponents),
                       'LDA': LinearDiscriminantAnalysis(n_components=self.nComponents)
                       }

        self.num_molecules = len(self.hash_space)
    
        self.num_properties = set([x.properties.shape[1] for x in self.hash_space.values()])
        
        self.properties = np.vstack([x.properties for x in self.hash_space.values()])


        # Store the PCA coordinates in the Molecules objects
        self.coordinates = None
        self.model = None

        self.dataframe = None


    def _store_pca_coordinates(self):
        """
        For each molecule in the chemspace, store the PCA coordinates in the molecule object 
        by indexing the self.coordinates with the molecule index (self.index_space).
        """

        current_index = 0
        for i, molecules in enumerate(self.hash_space.values()):
            length = molecules.num_molecules
            molecules.coordinates = self.coordinates[np.r_[current_index:length+current_index, :]]
            current_index += length


    def apply_model(self, model='PCA', scale=True, nComponents=None):
        self.nComponents = nComponents

        assert model in ('PCA', 'SVD', 'LDA') or type(model) in (PCA, TruncatedSVD, LinearDiscriminantAnalysis), "Choose a model (PCA, SVM, LDA) or enter a pretrainded model."
        if model in ('PCA', 'SVD', 'LDA'):
            # Standardize data
            scaler = StandardScaler()
            descriptors_scaled = scaler.fit_transform(self.properties)

            # Apply selected model to data (PCA, TruncatedSVD, LinearDiscriminantAnalysis).
            m = self.models[model] # get model
            print("Applying {} model...".format(model))
            self.coordinates = m.fit_transform(descriptors_scaled)[:,:nComponents]
            print(f"Number of components: {self.coordinates.shape[1]}")

        else:
            m = model
            assert type(model) in (PCA, TruncatedSVD, LinearDiscriminantAnalysis), "Selected model must be one of PCA, SVD, LDA from sklearn.\n e.g. sklearn.decomposition.PCA,\nsklearn.decomposition.TruncatedSVD,\nsklearn.discriminant_analysis.LinearDiscriminantAnalysis"
            #Standardize data
            scaler = StandardScaler()
            descriptors_scaled = scaler.fit_transform(self.properties)

            # Apply selected model to data (PCA, TruncatedSVD, LinearDiscriminantAnalysis).
            self._coordinates = m.fit_transform(descriptors_scaled)[:,:nComponents]

        self.model = m

        # Store the PCA coordinates in the Molecules objects
        self._store_pca_coordinates()

        return self.coordinates, self.model, scaler

    def __str__(self):
        
        return "<ChemSpace> object with {} sets of molecules.\nChemSpace dimensions: {}".format(self.num_molecules, self.properties.shape)
    

    def to_df(self):
        """
        Convert the ChemSpace object to a pandas.DataFrame.
        """
        # Dataframe to store PCA coordinates
        df = pd.DataFrame(self.coordinates, columns=['PC'+str(i+1) for i in range(self.nComponents)])
        # Names of the sets
        names = [x.setName for x in self.hash_space.values()]
        # Number of molecules in each set
        num_mol_by_set = [x.properties.shape[0] for x in self.hash_space.values()]
        # Add the set names and number of molecules to the dataframe
        listnames = []

        for i, name in enumerate(names):
            for j in range(num_mol_by_set[i]):
                listnames.append(name)
        
        df['setName'] = listnames
        

        # Add the properties to the dataframe
        for i in range(self.nComponents):
            # each column is a Pcomponent
            df['PC'+str(i+1)] = self.coordinates[:,i] # Column i is the Pcomponent

        self.dataframe = df
        
        return self.dataframe
        

def compute_distance(space_ref: np.array, space_target: np.array) -> float:
    """
    Compute distance between all molecules in space_ref and space_comp.

    Parameters
    ----------
    space_ref : np.array
        numpy array containing dimensional coordinates of molecules in space_ref.
    space_target : np.array
        numpy array containing dimensional coordinates of molecules in space_target.

    Returns
    -------
    float
        minimum distance between space_ref and space_target
    """

    distance_matrix = distance.cdist(space_ref, space_target, 'euclidean')

    return distance_matrix.min(axis=1)
        
def plot_distances(df: pd.DataFrame, ref_region: str=None, filename: str=None):
    """
    Plot the distances between the all molecules in the chemspace
    and the "target" region of the chemspace.
    """
    region = df[df.setName == ref_region]
    allsets = df[df.setName != ref_region]
    
    fig, ax = plt.subplots(1,1, figsize=(10, 6))
    plt.grid()

    ax1 = ax.scatter(data=region, x='PC1', y='PC2', marker='o', s=125, color="whitesmoke", alpha=1)
    ax2 = ax.scatter(data=allsets, x='PC1', y='PC2', marker='+', c=allsets['EuclDist'], cmap='plasma')
    ax3 = plt.scatter(data=allsets[allsets.setName == "IPTACOPAN"],
                    x='PC1',
                    y='PC2',
                    marker='X',
                    s=223)

    cbar = fig.colorbar(ax2, ax=ax)
    cbar.set_label('Euclidean Distance')

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    


    # LEGEND

    patch1 = mpatches.Patch(color='whitesmoke', label='Target Chemical Space')


    ipt = mlines.Line2D([], [], color='indigo', marker='X', linestyle='None',
                            markersize=10, label='Iptacopan', markeredgewidth=4)




    plt.legend(handles=[patch1, ipt], frameon=False, loc='upper left')

    plt.tight_layout()
    if filename is not None:
         plt.savefig(filename, dpi=300)
    plt.show()





def main():
    iptacopan = Chem.MolFromSmiles("CCOC1CCN(C(C1)C2=CC=C(C=C2)C(=O)O)CC3=C(C=C(C4=C3C=CN4)C)OC")
    ms_iptacopan = Molecules([iptacopan], "IPTACOPAN")

    iptacopan_s = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/iptacopan_scaffold.sdf")
    ms_iptacopan_s = Molecules(iptacopan_s, 'Iptacopan_scaffold')

    enamine_all = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/enamine_all_filtered.sdf")
    ms_enamine_all = Molecules(enamine_all, 'Enamine')

    enamine_in = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/enamine_filtered.sdf")
    ms_enamine_in = Molecules(enamine_in, 'EnamineIn')

    CFB = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/CFB_filtered.sdf")
    ms_CFB = Molecules(CFB, 'CFB')

    InSpace = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/InSpace_filtered.sdf")
    ms_InSpace = Molecules(InSpace, 'InSpace')

    # Create a TARGET ChemSpace object
    mw = np.arange(150, 450, 1)
    logP = np.arange(0, 3, 0.1)
    nrings = np.arange(0, 4, 1)
    pfi = np.arange(0, 5, 0.1)
    hbd = np.arange(0,3,1)
    
    region = np.array(np.meshgrid(logP, mw, hbd, nrings, pfi)).T.reshape(-1, 5)


    #print("Molecules region, Matrix format:")
    ms_region = Molecules(region)
    #print(ms_region.properties)


    # CREATE CHEMSPACE OBJECT CONTAINING ALL THE SETS OF MOLECULES
    cs = ChemSpace([ms_iptacopan, ms_iptacopan_s, ms_enamine_all, ms_enamine_in, ms_CFB, ms_InSpace, ms_region])
    print(cs, end='\n\n')

    coords, model, scaler = cs.apply_model(model='PCA', nComponents=2)

    # sum all the sets of molecules in the ChemSpace object

    results = cs.to_df()
    print(results.head())


    # Compute the distance between all sets in the ChemSpace and the target region
    ipta = results[results.setName == 'IPTACOPAN']

    # Select the target region
    region = results[results.setName == ms_region.setName]
    # Select the remaining sets
    dfcopy = results[results.setName != ms_region.setName]

    d = compute_distance(ipta[['PC1', 'PC2']].to_numpy(), region[['PC1', 'PC2']].to_numpy())
    print(d)


    st = time.time()
    

    for set in dfcopy['setName'].unique():
        print(f"Current set: {set}")
        curr_group = dfcopy[dfcopy.setName == set]

        # Calculate the Euclidean distance between each molecule in the chemspace using compute_distance
        distances = compute_distance(curr_group[['PC1', 'PC2']].to_numpy(), region[['PC1', 'PC2']].to_numpy())
        results.loc[results.setName == set, 'EuclDist'] = distances


    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    # Plot the distances between the all molecules in the chemspace
    # and the "target" region of the chemspace.
    plot_distances(results, ref_region=ms_region.setName) #filename='/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/distances.png')
    



if __name__ == "__main__":
    main()
