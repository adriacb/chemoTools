from rdkit import Chem
from rdkit.Chem import PandasTools

from rhea.chem.molecules import Molecules
from rhea.chem.molecule import Molecule

ref, = [m for m in Molecules(infile="/chemotargets/research/SITALA/IPTACOPAN/step1_PREPARATION_str/6rav-pdb-bundle1.JGQ.A.303.sdf")]
print(ref)


df = PandasTools.LoadSDF("/chemotargets/research/SITALA/IPTACOPAN/docking_rDock/RefLig/ref_lig.sdf")
df["rmsd"] = df.apply(lambda row: Molecule(row["ROMol"]).rmsd(ref.molecule), axis=1)
df = df.sort_values("rmsd")
df.head()