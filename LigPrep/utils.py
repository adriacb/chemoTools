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