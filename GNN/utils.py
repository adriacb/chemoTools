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


def AutoMol(mol):
    if type(mol) == str:
        assert Chem.MolFromSmiles(mol, sanitize=False) is not None, f"Cannot convert {mol} to SMILES. Wrong format."
        mol = Chem.MolFromSmiles(mol)
    else:
        assert type(mol) == Chem.rdchem.Mol, f"Input rdkit molecule needed."
        mol = mol
    return mol