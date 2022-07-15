import os
import sys
import re 

import pathlib
import argparse
import pickle
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from timeit import default_timer as timer
from datetime import timedelta
import time
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

#IPythonConsole.ipython_useSVG = True
#IPythonConsole.molSize = (400, 400)


def path_Extension(path):
    # function to return the file extension
    file_extension = pathlib.Path(path).suffix
    print("File Extension: ", file_extension)
    return file_extension

def listPath( path=os.getcwd(), extension=''):
    # list to store files
    res = []
    # Iterate directory
    for file in os.listdir(path):
        # check only text files
        if file.endswith(extension) and file.startswith('out'):
            res.append(file)
    return res

def extractMolFromSDFById(sdf: str, id: str) -> pd.DataFrame:
    # extract molecule from sdf by id
    sdf = PandasTools.LoadSDF(sdf, molColName='ROMol')
    moldf = sdf[sdf['ID'] == id]
    return moldf


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

def SDFWriter(pd_df: pd.DataFrame, out_file: str, molColName='ROMol'):
    with Chem.SDWriter(out_file) as writer:
        for index, row in pd_df.iterrows():
            mol = row[molColName]
            properties = [x for x in list(pd_df.columns) if x not in ['ROMol',molColName]]
            for p in properties:
                mol.SetProp(p, str(row[p]))
            writer.write(mol)
    writer.close()

def rDockformat(sdf: str):
    '''
    This function reads a text file and every
    line checks if a parenthesis with a number inside is present using regex
    If so, it removes the parenthesis.
    For every line, write the line to a new file.
    '''
    with open(sdf, 'r') as f:
        lines = f.readlines()
    with open('output.sdf', 'w') as f:
        for line in lines:
            if re.search(r'\(\d+\)', line):
                line = re.sub(r'\(\d+\)', '', line)
            f.write(line)
    f.close()


def SDFWriter2(pd_df: pd.DataFrame, out_file: str, molColName='ROMol'):
    path = os.path.dirname(out_file)
    i = 0
    for index, row in pd_df.iterrows():
        print(path+'s_{}.sdf'.format(i))
        with Chem.SDWriter(path+'s_{}.sdf'.format(i)) as writer:
            mol = row[molColName]
            properties = [x for x in list(pd_df.columns) if x not in ['ROMol',molColName]]
            for p in properties:
                mol.SetProp(p, str(row[p]))
            writer.write(mol)
        i += 1
        writer.close()
    
    # merge all sdf files
    os.system("cat {}s_*.sdf > {}".format(path, out_file))

# def RMSD(mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol) -> float:
#     return Chem.rdMolAlign.AlignMol(mol1, mol2, maxIters=1000)

