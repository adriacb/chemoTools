from .utils import *

class Combine:
    """
    Combines two molecules, adding the fragment to the reference region using the indeces of the nearest atom.
    It returns the combined molecule as Chem.rdchem.Mol.
    """

    def __init__(self, ref: Chem.rdchem.Mol, frag: Chem.rdchem.Mol,
                 i_ref: int, i_frag: int = None):
        """
        Initialize the class.

        Parameters
        ----------
        ref : Chem.rdchem.Mol
            rdkit reference molecule.
        frag : Chem.rdchem.Mol
            rdkit fragment molecule.
        i_ref : int
            index of the nearest atom in the reference molecule.
        i_frag : int
            index of the nearest atom in the fragment molecule. (None if not provided)
        """
        self.ref = ref
        self.frag = frag
        self.i_ref = i_ref
        self.i_frag = None
        self.ats_ref, self.coor_ref = self._coords_atoms(ref)
        self.ats_frag, self.coor_frag = self._coords_atoms(frag)


        self.molecule = self._pipeline()
        

    def _combine(self, i: int, j: int) -> Chem.rdchem.Mol:
        """
        Combine two molecules by adding the fragment to the reference region using the indeces of the nearest atom.

        Parameters
        ----------
        ref : Chem.rdchem.Mol
            rdkit reference molecule.
        frag : Chem.rdchem.Mol
            rdkit fragment molecule.
        
        Returns
        -------
        Chem.rdchem.Mol
            rdkit molecule containing the combined molecules.
        """
        # Create a combined molecule
        combo = Chem.CombineMols(self.ref, self.frag)
        emol = Chem.EditableMol(combo)

        # num will be the index of the last atom in the reference region
        num = self.ref.GetNumAtoms()

        # Add a single bond between the two atoms using the indeces of the nearest atom
        emol.AddBond(i, num+j, order=Chem.rdchem.BondType.SINGLE)

        # Convert the combined molecule to a rdkit molecule
        mol = emol.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    

    def _coords_atoms(self, mol: Chem.rdchem.Mol) -> Tuple[np.array, np.array]:
        """
        Compute the coordinates of the atoms in a molecule.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            rdkit molecule.

        Returns
        -------
        list
            list of atom coordinates.
        np.array
            numpy array containing the coordinates of the atoms.
        """


        atoms_ = []
        coords_ = []
        conf = mol.GetConformer()
        for i, atom in enumerate(mol.GetAtoms()):
            positions = conf.GetAtomPosition(i) # Get the coordinates of the atom
            atoms_.append(atom) # append atom object
            coords_.append(np.array((positions.x, positions.y, positions.z))) # append coordinates

        coords_ = np.vstack([x for x in coords_]) # stack coordinates
        
        return atoms_, coords_

    
    def _find_nearest_coord(self) -> Tuple[int, int]:
        """
        Find the nearest coordinate in the reference region to a coordinate in the fragment region.

        Parameters
        ----------
        self.coor_frag : np.array
            numpy array containing the coordinates of the fragment region.
        self.coor_ref : np.array
            numpy array containing the coordinates of the reference region.
        
        Returns
        -------
        int
            index of the nearest coordinate in the reference region.
        np.array
            numpy array containing the coordinates of the nearest coordinate in the reference region.
        int
            index of the nearest coordinate in the fragment region.
        np.array
            numpy array containing the coordinates of the nearest coordinate in the fragment region.
        """


        distance_matrix = distance.cdist(self.coor_ref, self.coor_frag, 'euclidean')
        frag_dist = distance_matrix.min(axis=0)
        ref_dist = distance_matrix.min(axis=1)
        
        # Index of "ith" element of frag_coords having min distance
        index1 = np.where(frag_dist == np.amin(frag_dist))
        # Index of "ith" element of ref_coords having min distance
        index2 = np.where(ref_dist == np.amin(ref_dist))
        
        return index2[0][0], self.coor_ref[index2[0][0]], \
            index1[0][0], self.coor_frag[index1[0][0]]

    
    def _pipeline(self) -> Chem.rdchem.Mol:
        """
        Pipeline for combining molecules.

        Returns
        -------
        Chem.rdchem.Mol
            rdkit molecule containing the combined molecules.
        """

        if self.ats_frag is None:
            # Find the nearest coordinate in the reference region to a coordinate in the fragment region
            i_ref, coor_ref, self.i_frag, coor_frag = self._find_nearest_coord()

        try:
            # Create a combined molecule
            combo = self._combine(int(i_ref), int(self.i_frag))
        except:
            print("Error: Could not combine molecules.")
            return None
        
        return combo
    
