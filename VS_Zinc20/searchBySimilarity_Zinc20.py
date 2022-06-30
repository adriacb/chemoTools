from timeit import default_timer as timer
from datetime import timedelta
import pandas as pd

import pickle

from rhea.db import zinc20
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw 
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
IPythonConsole.ipython_useSVG = True
IPythonConsole.molSize = (400, 400)

target_smiles = "CCOC1CCN(C(C1)C2=CC=C(C=C2)C(=O)O)CC3=C(C=C(C4=C3C=CN4)C)OC"

mol_query = Chem.MolFromSmiles(target_smiles)

fp_query = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol_query, radius=2, nBits=1024)

results = pd.DataFrame({'ids':[], 'smiles': [], 'scores': []})

threshold = 0.5


total_db = 0
total_found = 0


if __name__ == "__main__":

    start = timer()

    print(f"Start Time: {timedelta(seconds=start)}\n\nTarget: {target_smiles}\nSimilarity threshold: {threshold}\n\n")


    freq_tanim = dict()

    found = 0

    with zinc20.Zinc20() as db:

        
        db.cursor.execute("SELECT id, rdkit.mol_send(mol), mol as SMILE FROM instock.molecules")
        result = db.cursor.fetchall()

        assert len(result) > 0, "No molecules in the database"
        
        ids = ['ZINC000223246939', 'ZINC000223246892', 'ZINC001911571385', 'ZINC000223181853', 'ZINC000223181995', 'ZINC000223258195']

        # iterate over the results and calculate the tanimoto similarity
        for res in result:
            #increment the counter of the total fragments/molecules of the DB
            total_db += 1

            
            id, mol, smile = res
            mol = Chem.Mol(mol.tobytes())
            # compute fingerprint
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    mol, radius=2, nBits=1024)

            if id in ids:
                found += 1


            # calculate tanimoto similarity
            tanimoto = DataStructs.TanimotoSimilarity(fp_query, fp)
            print(f"Curr Time: {timedelta(seconds=timer())} | Curr id: {id}, tanimoto: {tanimoto}", end='\r')

            freq_tanim[round(tanimoto, 2)] = freq_tanim.get(round(tanimoto, 2), 1) + 1
            
            # if the similarity is above the threshold, add the result to the results dataframe and print it
            if tanimoto >= threshold:
                total_found += 1

                print(f"\nid: {id}, smile: {smile}, tanimoto: {tanimoto}")
                #Draw.MolToFile(mol, "./mols/%d.svg" % id)
                results = results.append({'ids': id, 'smiles': smile, 'scores': tanimoto}, ignore_index=True)
            else:
                pass
        
        # if the dataframe is not empty, save it to a csv file and print the results
        if not results.empty:
            results.to_csv("./results/similarity_results.csv", index=False)
            print(results.head(25))
        
    with open('tanim_freq.pickle', 'wb') as handle:
        pickle.dump(freq_tanim, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    end = timer()
    print(f"\nElapsed time: {timedelta(seconds=end-start)}\n")
    print(f"Total size of the DB: {total_db}")
    print(f"Number of similar fragments: {total_found}")
    print(f"{found}")
