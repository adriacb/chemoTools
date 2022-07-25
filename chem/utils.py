import os

import pickle
import uuid
import time

from typing import Tuple, Union

import pandas as pd
import numpy as np

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, Descriptors

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE

from scipy.spatial import distance

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

def compute_distance(space_ref: np.array, space_target: np.array) -> float:
    """
    Compute min distance between all molecules in space_ref and space_comp.

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
    Plot the distances between the molecules in the chemspace
    and the "target" region of the chemspace.

    Parameters
    ----------
        df : pd.DataFrame
            pandas dataframe containing the coordinates of the molecules.
        ref_region : str
            name of the reference region.
        filename : str
            name of the file to save the plot.        
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
         plt.savefig(f"{filename}_distances.png", dpi=300)
    plt.show()

    # Histogram of distances by set name
    ax = allsets.pivot(columns='setName', values='EuclDist').plot(kind='hist', 
                              bins = 100, 
                              figsize=(12,8),
                              alpha = 0.6, grid=True)
    vline = ax.vlines(x = 2.525, ymin = 0, ymax = 80,
            colors = 'black',
            label = 'Iptacopan')
    ax.set_xlabel("Euclidean Distances")
    leg1 = ax.legend(['CFB', 'Enamine Real', 'Enamine Real - In Space', 'Iptacopan scaffold'], loc='upper right')
    leg2 = ax.legend(handles=[vline], loc='lower right')
    ax.add_artist(leg1)
    if filename is not None:
        plt.savefig(f"{filename}_hist.png", dpi=300)
    plt.show()


def plot_explained_variance(pca, dims):

    plt.bar(range(1,dims+1), pca.explained_variance_ratio_,
            alpha=0.5,
            align='center')
    plt.step(range(1,dims+1), np.cumsum(pca.explained_variance_ratio_),
            where='mid',
            color='red')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal Components')
    plt.show()


def coords_atoms(mol: Chem.rdchem.Mol) -> Tuple[list, np.array]:
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

def find_nearest_coord(frag_coords: np.array, ref_coords: np.array) -> Tuple[int, np.array, int, np.array]:
    """
    Find the nearest coordinate in the reference region to a coordinate in the fragment region.

    Parameters
    ----------
    frag_coords : np.array
        numpy array containing the coordinates of the fragment region.
    ref_coords : np.array
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


    distance_matrix = distance.cdist(ref_coords, frag_coords, 'euclidean')
    frag_dist = distance_matrix.min(axis=0)
    ref_dist = distance_matrix.min(axis=1)
    
    # Index of "ith" element of frag_coords having min distance
    index1 = np.where(frag_dist == np.amin(frag_dist))
    # Index of "ith" element of ref_coords having min distance
    index2 = np.where(ref_dist == np.amin(ref_dist))
    
    return index2[0][0], ref_coords[index2[0][0]], \
           index1[0][0], frag_coords[index1[0][0]]


def combine(ref: Chem.rdchem.Mol, frag: Chem.rdchem.Mol, i: int, j: int) -> Chem.rdchem.Mol:
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
    combo = Chem.CombineMols(ref, frag)
    emol = Chem.EditableMol(combo)

    # num will be the index of the last atom in the reference region
    num = ref.GetNumAtoms()

    # Add a single bond between the two atoms using the indeces of the nearest atom
    emol.AddBond(i,num+j, order=Chem.rdchem.BondType.SINGLE)

    # Convert the combined molecule to a rdkit molecule
    mol = emol.GetMol()
    Chem.SanitizeMol(mol)
    return mol