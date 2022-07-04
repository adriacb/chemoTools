from rdkit import Chem
from rdkit.Chem import PandasTools

from rhea.chem.molecules import Molecules
from rhea.chem.molecule import Molecule

ref, = [m for m in Molecules(infile="/chemotargets/research/SITALA/IPTACOPAN/step1_PREPARATION_str/6rav-pdb-bundle1.JGQ.A.303.sdf")]
print(ref)


df = PandasTools.LoadSDF("/chemotargets/research/SITALA/IPTACOPAN/docking_rDock/Chembl/CHEMBL4594448/splited/sorted_docked_docked_s5.sdf.sd.sdf")
df["rmsd"] = df.apply(lambda row: Molecule(row["ROMol"]).rmsd(ref.molecule), axis=1)
df = df.sort_values("rmsd")
Chem.PandasTools.WriteSDF(df, "/chemotargets/research/SITALA/IPTACOPAN/docking_rDock/Chembl/CHEMBL4594448/splited/bestRMSD_5.sdf", molColName='ROMol', properties=list(df.columns))
print(df["rmsd"])