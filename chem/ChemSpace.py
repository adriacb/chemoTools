from utils import *
from Molecules import Molecules


class ChemSpace:
    """
    Parameters
        chemspace - list of <Molecules> objects
    ===========================
    """

    def __init__(self, chemspace: list[Molecules]):

        assert type(chemspace) == list, "Chemspace must be a list of Molecules objects."
        assert len(chemspace) > 0, "Chemspace must contain at least one molecule."
        assert all(isinstance(n, Molecules) for n in chemspace), "Chemspace must be a list of Molecules objects."


        self.hash_space = dict()
        self.index_space = list()
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


        self.properties = np.vstack([x.properties for x in self.hash_space.values()])
        self.num_properties = self.properties.shape[1]
        self.num_molecules = self.properties.shape[0]
        
        # Store the PC coordinates in the Molecules objects
        self.coordinates = None
        # Store the model 
        self.model = None

        self.dataframe = None


    def _store_pc_coordinates(self):
        """
        For each set of <Molecules> in the chemspace, store the PCA coordinates in the <Molecule> object 
        by indexing the self.coordinates with the molecule index (self.index_space).
        """

        current_index = 0
        for molecules in self.hash_space.values():
            length = molecules.num_molecules
            molecules.coordinates = self.coordinates[np.r_[current_index:length+current_index, :]]
            current_index += length


    def tSNE(self, model='t-SNE', scaler=None, nComponents=2, perplexity=30,  
                early_exaggeration=12.0, learning_rate='warn', n_iter=1000, 
                n_iter_without_progress=300, min_grad_norm=1e-07, 
                metric='euclidean', metric_params=None, init='warn', 
                verbose=0, random_state=None, method='barnes_hut', angle=0.5, 
                n_jobs=None):
        """
        Parameters
        ==========
        model: str
            The model to be used for dimensionality reduction.
        scaler: sklearn.preprocessing.StandardScaler
            The scaler to be used for dimensionality reduction.
        Others:
            seehttps://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html for more information.
        
        Returns
        =======
        self.coordinates: numpy.ndarray
            The TSNE coordinates of the molecules.
        self.model: sklearn.manifold.TSNE
            The TSNE model.
        scaler: sklearn.preprocessing.StandardScaler
            The scaler used for dimensionality reduction.
        """

        self.nComponents = nComponents

        if model in ('t-SNE'):
            # Standardize data
            scaler = StandardScaler()
            descriptors_scaled = scaler.fit_transform(self.properties)

            # Apply selected model to data (PCA, TruncatedSVD, LinearDiscriminantAnalysis).
            m = TSNE(nComponents, perplexity, early_exaggeration, learning_rate, n_iter, n_iter_without_progress,
                        min_grad_norm, metric, metric_params, init, 
                        verbose, random_state, method, angle, n_jobs) # get model
            print("Applying {} model...".format(model))
            self.coordinates = m.fit_transform(descriptors_scaled)[:,:nComponents]
            print(f"Number of components: {self.coordinates.shape[1]}")

        else:
            assert type(model) in (TSNE), "Selected model must be TSNE from sklearn."
            m = model
            
            if scaler is None:
                #Standardize data
                scaler = StandardScaler()
                descriptors_scaled = scaler.fit_transform(self.properties)
            else:
                descriptors_scaled = scaler.transform(self.properties)
            # Apply selected model to data (PCA, TruncatedSVD, LinearDiscriminantAnalysis).
            self._coordinates = m.fit_transform(descriptors_scaled)[:,:nComponents]

        self.model = m

        # Store the PCA coordinates in the Molecules objects
        self._store_pc_coordinates()
        return self.coordinates, self.model, scaler

    def PCA(self, model='PCA', scaler=None, nComponents=None, copy=True, whiten=False, svd_solver='auto', 
            tol=0.0, iterated_power='auto', n_oversamples=10, power_iteration_normalizer='auto', random_state=None):
        
        """
        Parameters
        ==========
        model: str
            The model to be used for dimensionality reduction.
        scaler: sklearn.preprocessing.StandardScaler
            The scaler to be used for dimensionality reduction.
        Others:
            see https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html for more information.
        
        Returns
        =======
        self.coordinates: numpy.ndarray
            The PCA coordinates of the molecules.
        self.model: sklearn.decomposition.PCA
            The PCA model.
        scaler: sklearn.preprocessing.StandardScaler
            The scaler used for dimensionality reduction.
        """
        
        self.nComponents = nComponents
        

        if model in ('PCA'):
            # Standardize data
            scaler = StandardScaler()
            descriptors_scaled = scaler.fit_transform(self.properties)

            # Apply selected model to data (PCA).
            m = PCA(n_components=self.nComponents, copy=copy, whiten=whiten, svd_solver=svd_solver, tol=tol, 
            iterated_power=iterated_power, n_oversamples=n_oversamples, power_iteration_normalizer=power_iteration_normalizer, 
            random_state=random_state) 

            print("Applying PCA model...")
            self.coordinates = m.fit_transform(descriptors_scaled)[:,:nComponents]
            print(f"Number of components: {self.coordinates.shape[1]}")

        else:
            assert type(model) in (PCA), "Enter a pretrainded PCA model."
            m = model
            if scaler is None:
                #Standardize data
                scaler = StandardScaler()
                descriptors_scaled = scaler.fit_transform(self.properties)
            
            else:
                descriptors_scaled = scaler.transform(self.properties)

            # Apply selected model to data (PCA, TruncatedSVD, LinearDiscriminantAnalysis).
            self._coordinates = m.fit_transform(descriptors_scaled)[:,:nComponents]

        # Store the model
        self.model = m
        # Store the PCA coordinates in the Molecules objects
        self._store_pc_coordinates()

        return self.coordinates, self.model, scaler

    def __str__(self):
        
        return "<ChemSpace> object with {} sets of molecules.\nChemSpace dimensions: {}".format(len(self.hash_space.keys()), self.properties.shape)
    

    def to_df(self):
        """
        Convert the ChemSpace object to a pandas.DataFrame.
        """
        assert self.coordinates is not None, "PCA or TSNE coordinates must be computed first."

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
    

    def to_sdf(self, filename, molCol='ROMol', propsCols=None):
        """
        Convert the ChemSpace object to an sdf file.
        """
        if self.dataframe is None:
            df = self.to_df()
            propsCols = list(df.columns)
            PandasTools.WriteSDF(df , f'{filename}.sdf', molColName = molCol, properties = propsCols)
        else:
            propsCols = list(self.dataframe.columns)
            PandasTools.WriteSDF(self.dataframe , f'{filename}.sdf', molColName = molCol, properties = propsCols)
        
        print("SDF file saved to {}.sdf\n{} Molecules converted.".format(filename, self.num_molecules))
        
