from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

from rhea.chem.molecules import Molecules
from rhea.chem.molecule import Molecule

ref, = [m for m in Molecules(infile="/chemotargets/research/SITALA/IPTACOPAN/step1_PREPARATION_str/6rav-pdb-bundle1.JGQ.A.303.sdf")]
print(ref)


df = PandasTools.LoadSDF("/home/adria/Documents/chemoTools/LigPrep/best_docking.sdf")
print(df.info())
df["SCORE.INTER"] = df["SCORE.INTER"].astype(float)
df["rmsd"] = df.apply(lambda row: Molecule(row["ROMol"]).rmsd(ref.molecule), axis=1)
#df = df.sort_values("rmsd")
#Chem.PandasTools.WriteSDF(df, "/chemotargets/research/SITALA/IPTACOPAN/docking_rDock/Chembl/CHEMBL4594448/splited/bestRMSD_5.sdf", molColName='ROMol', properties=list(df.columns))
df = df.sort_values("SCORE.INTER", ascending=True)
print(min(df["SCORE.INTER"]))
print(df[["SCORE.INTER", "rmsd"]])


# def RMSD(mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol) -> float:
#     return Chem.rdMolAlign.AlignMol(mol1, mol2, maxIters=1000)


# df['rmsd'] = df.apply(lambda row: RMSD(row["ROMol"], ref.molecule), axis=1)
# df = df.sort_values("rmsd")
# print(df['rmsd'])

