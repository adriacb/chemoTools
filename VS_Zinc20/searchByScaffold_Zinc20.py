from timeit import default_timer as timer
from datetime import timedelta
import pandas as pd

import pickle

from rhea.db import zinc20
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw 
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
IPythonConsole.ipython_useSVG = True
IPythonConsole.molSize = (400, 400)

target_smiles = "CCOC1CCN(C(C1)C2=CC=C(C=C2)C(=O)O)CC3=C(C=C(C4=C3C=CN4)C)OC"

scaffold_query = Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(target_smiles))

target_inchiKey = Chem.MolToInchiKey(scaffold_query)
#print(scaffold_query)
#Draw.MolToFile(scaffold_query, "./test.svg")
results = pd.DataFrame({'ids':[], 'smiles': [], 'mol': [], 'inchikey': [], 'scaffold': []})


total_db = 0
total_found = 0


if __name__ == "__main__":

    start = timer()

    print(f"Start Time: {timedelta(seconds=start)}\n\nTarget InchiKey: {target_inchiKey}\n\n")
    total_failed = 0
    with zinc20.Zinc20() as db:
	#db.cursor.execute("SELECT id, rdkit.mol_send(mol), mol as SMILE FROM instock.molecules")
        db.cursor.execute("SELECT id, rdkit.mol_send(mol), mol as SMILE, inchikey, scaffold from instock.molecules")
        result = db.cursor.fetchall()
        assert len(result) > 0, "No molecules in the database"
        counter = 0
        #print(result)
        # iterate over the results and calculate the tanimoto similarity
        for res in result:
            #increment the counter of the total fragments/molecules of the DB
            total_db += 1

            id, mol, smile, inchi, scaff = res
            print(f"Curr Time: {timedelta(seconds=timer())} | Curr id: {id}, Curr InchiKey: {inchi}, curr Scaffold: {scaff}", end='\r')
            if scaffold_query == scaff or id in ('ZINC82098524', 'ZINC223246892'):
                total_found += 1
                results = results.append({'ids': id, 'smiles': smile, 'mol':mol, 'inchikey': inchi, 'scaffold': scaff}, ignore_index=True)

        # if the dataframe is not empty, save it to a csv file and print the results
        if not results.empty:
            results.to_csv("scaffoldMatch_results_v2.csv", index=False)
            print(results.head(25))
    end = timer()
    print(f"\nElapsed time: {timedelta(seconds=end-start)}\n")
    print(f"Total size of the DB: {total_db}")
    print(f"Number of similar fragments: {total_found} failed {total_failed}")
