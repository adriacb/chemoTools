import os
import sys

import pathlib
import argparse
import pickle

from timeit import default_timer as timer
from datetime import timedelta

import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

def GetBestDockingScores(file, inter, rest, id, out) -> pd.DataFrame:
   
    SDF = PandasTools.LoadSDF(file)
    SDF["SCORE.INTER"] = SDF["SCORE.INTER"].astype(float)
    SDF["SCORE.RESTR"] = SDF["SCORE.RESTR"].astype(float)
    SDF = SDF[(SDF["SCORE.INTER"] <= inter) & (SDF["SCORE.RESTR"]<= rest)]
    SDF = SDF.sort_values("SCORE.INTER", ascending=True).drop_duplicates(id)
    PandasTools.WriteSDF(SDF, out, molColName='ROMol', properties=list(SDF.columns))
    
    return SDF

def PlotDistributions(df):
    bins=len(df)
    sns.histplot(data=df, x="SCORE.INTER", bins=bins)
    plt.show()
    sns.histplot(data=df, x="SCORE.RESTR", bins=bins)
    plt.show()

def path_Extension(path) -> str:
    # function to return the file extension
    file_extension = pathlib.Path(path).suffix
    print("File Extension: ", file_extension)
    return file_extension

def main():

    print("rdkit version: ", rdkit.__version__, end='\n')

    parser = argparse.ArgumentParser()
    parser.add_argument("-sd", "--sdf", action="store", dest="sdf")
    parser.add_argument("-n", "--name", action="store", dest="name", type=str)
    parser.add_argument("-i", "--inter", action="store", dest="inter", type=float)
    parser.add_argument("-r", "--rest", action="store", dest="rest", type=float)
    parser.add_argument("-o", "--out", action="store", dest="out")
    args = parser.parse_args()

  
    if args.sdf:

        assert (path_Extension(args.sdf) in ('.sdf','.sd')), "Wrong input file type (.sdf, .sd)."
        start = timer()
        print(f"Start Time: {timedelta(seconds=start)}\n")

        output = args.out if args.out else 'best_docking.sdf'
        name = args.name if args.name else 'Name'
        inter = args.inter if args.inter else 5.0
        restr = args.rest if args.rest else 0.5

        bests = GetBestDockingScores(args.sdf, inter=inter, rest=restr, id=name, out=output)
        PlotDistributions(df=bests)


    else:
        sys.exit()

    end = timer()
    print(f"\nElapsed time: {timedelta(seconds=end-start)}\n")

if __name__ == "__main__":
    main()
