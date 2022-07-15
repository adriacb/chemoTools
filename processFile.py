#!usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "adria.cabello@chemotargets.com"
__date__ = "2022-07-12"

import bz2
import glob
import hashlib
import itertools
import logging
import multiprocessing
import os
import time

from argparse import ArgumentParser

import numpy as np

from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

from utils import *


class reader:

    def __init__(self, path, extension):
        self.dirname = os.path.dirname(path)
        self.extension = extension


    def __str__(self):
        return f"{self.dirname}, {self.extension}"


    def iter_infiles(self):

        # Get files
        infiles = glob.glob("{}/*{}".format(self.dirname, self.extension))
        print(infiles)
        infiles.sort()

        for infile in infiles:
            yield infile 

   
    def transform(self, chunk_size=4, cpu=14,
               outfile="out.txt"):

        processes = list()

        for f in self.iter_infiles():
            print(f)
        
            with open(f, "r") as file:
                chunk = True
                ichunk = 0

                while chunk:

                    ichunk += 1
                    #print(ichunk, f)


                    while(len(processes) >= cpu):
                        time.sleep(1)
                        processes = self.clean_processes(processes, outfile)

                    try:
                        chunk = list(itertools.islice(file, chunk_size))
                        #print(chunk)
                    except IndexError:
                        break

                    #print(chunk[0].split(b"\t")[-1])

                    if not chunk:
                        break


                    #h = hashlib.sha256(f"{outfile}{chunk_size}{threshold}".encode('utf8')).hexdigest()[:7]
                    #outfile = f"enamine.{os.path.basename(infile).split('.')[0]}.{h}.{ichunk:03d}.out"
                    if os.path.isfile(outfile):
                        continue


                    p = multiprocessing.Process(target=self.prep_lines,
                                                kwargs={
                                                    "lines": list(chunk),
                                                    "id": ichunk,
                                                    "outfile": outfile})
                    p.start()
                    processes.append(p)

        [p.join() for p in processes]
        self.clean_processes(processes, outfile)


        


    def search(self, smiles, threshold=0.5, chunk_size=1000000, cpu=14,
               outfile="out.txt"):
        """
        File based similarity search """


        self.done = 0


        fpt = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2)
        print(fpt)

        processes = []

        for infile in self.iter_infiles():
            print(infile)

            with bz2.BZ2File(infile, "r") as f:
                header = f.readline()

                chunk = True
                ichunk = 0

                while chunk:

                    ichunk += 1
                    print(ichunk, infile)


                    while(len(processes) >= cpu):
                        time.sleep(1)
                        processes = self.clean_processes(processes, outfile)

                    try:
                        chunk = list(itertools.islice(f, chunk_size))
                    except IndexError:
                        break

                    #print(chunk[0].split(b"\t")[-1])

                    if not chunk:
                        break

                    if os.path.isfile(outfile):
                        continue


                    p = multiprocessing.Process(target=self.search_lines,
                                                kwargs={
                                                    "lines": list(chunk),
                                                    "fpt": fpt,
                                                    "threshold": threshold,
                                                    "id": ichunk,
                                                    "outfile": outfile})
                    p.start()
                    processes.append(p)

        [p.join() for p in processes]
        self.clean_processes(processes, outfile)


    def clean_processes(self, processes, outfile):

        for p in processes:
            if not p.is_alive():
                processes.remove(p)
            else:
                pass
                #print(p.get())
                #with open(outfile, "aw") as f:

        return processes

    def search_lines(self, lines, fpt, threshold=0.5, id=None, outfile=None):

        found = []
        for i, line in enumerate(lines):
            print(line.split(b"\s")[0])
            fpt2 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(line.split(b"\s")[0]), 2)
            sim = DataStructs.TanimotoSimilarity(fpt,fpt2)
            if sim >= threshold:
                found.append((sim, line))

        with open(outfile, "w") as f:
            for sim, line in found:
                f.write(f"{line.decode().strip()}\t{sim:0.3f}\n")


    def prep_lines(self, lines, id=None, outfile=None):

        transformed = []
        for i, line in enumerate(lines):
            print(line.split()[0])
            fpt = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(line.split()[0]), 2)
            transformed.append((line.split()[0], fpt))
            
        with open(outfile, "w") as f:
            for smile, fps in transformed:
                f.write(f"{smile}\t{fpt}\n")

def main():

    parser = ArgumentParser()
    parser.add_argument('-i', '--input', help="input path", action="store", dest="input")
    parser.add_argument('-o', '--output', help="output file", action="store", dest="out")
    parser.add_argument('--multi', type=int, help="Multi threading", default=0)

    group = parser.add_argument_group('options')

    args = parser.parse_args()


    r = reader(args.input, '.smi')
    r.transform()


main()