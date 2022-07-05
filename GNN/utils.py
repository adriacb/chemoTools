import os
import sys

import pathlib
import argparse
import pickle
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from timeit import default_timer as timer
from datetime import timedelta
import textwrap

import pandas as pd

from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw 
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import numpy as np


from torch_geometric.data.dataloader import DataLoaders
from torch_geometric.nn import MessagePassing
from torch_scatter import scatter_add
from torch_geometric.utils import add_self_loops, degree

def AutoMol(mol):
    if type(mol) == str:
        assert Chem.MolFromSmiles(mol, sanitize=False) is not None, f"Cannot convert {mol} to SMILES. Wrong format."
        mol = Chem.MolFromSmiles(mol)
    else:
        assert type(mol) == Chem.rdchem.Mol, f"Input rdkit molecule needed."
        mol = mol
    return mol

def get_atom_features(mol: Chem.rdchem.Mol):
    """
    Get atom features of a molecule.
    """
    atomic_number = list()
    num_hs = list()

    for atom in mol.GetAtoms():
        atomic_number.append(atom.GetAtomicNum())
        num_hs.append(atom.GetTotalNumHs(includeNeighbors=True))

    return torch.tensor([atomic_number, num_hs]).t()

def get_edge_index(mol: Chem.rdchem.Mol):
    """
    Get edge index of a molecule.
    """
    row, col = [], []
    
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        
    return torch.tensor([row, col], dtype=torch.long)


#function decorator to measure elapsed time
def timeit(func):
    def timed(*args, **kwargs):
        start = timer()
        print(f"Start time: {start}")
        result = func(*args, **kwargs)
        end = timer()
        print(f"{func.__name__} took {timedelta(seconds=end-start)}")
        return result
    return timed