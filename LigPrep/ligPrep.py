from time import time
from utils import *

def gen3d(mol: Chem.rdchem.Mol, n_conf: int) -> pd.DataFrame:
    '''
    Generate 3D conformers for a molecule.
    -----------------------------------------------------------------
    Parameters:
    mol: rdkit.Chem.rdchem.Mol
        molecule to generate 3D conformers for
    n_conf: int
        number of 3D conformers to generate
    -----------------------------------------------------------------
    Returns:
    df: pd.DataFrame
        dataframe with 3D conformers


        Based on:
        https://iwatobipen.wordpress.com/2021/01/31/generate-conformers-script-with-rdkit-rdkit-chemoinformatics/ 
    '''
    m = Chem.AddHs(mol, addCoords=True)
    refmol = Chem.AddHs(Chem.Mol(m))   
    df = pd.DataFrame({'ROMol': [], 'Energy': [], 'CID': []})
    
    p = Chem.rdDistGeom.ETKDGv2()
    p.pruneRmsThresh = 0.1
    cids = Chem.rdDistGeom.EmbedMultipleConfs(m, n_conf, p)
    
    mp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0, mmffVariant='MMFF94s')

    Chem.rdMolAlign.AlignMolConformers(m)
    
    
    for cid in cids:

        ff = AllChem.MMFFGetMoleculeForceField(m, mp, confId=cid)
        energy = ff.CalcEnergy()
        mol = Chem.Mol(refmol, True)
        mol.AddConformer(m.GetConformer(cid))
        df = df.append({'ROMol': mol, 'Energy': energy, 'CID': cid}, ignore_index=True)

    return df


def iterSMILES(smiles: str, n_conf=1) -> pd.DataFrame:
    '''
    It will generate 3D conformers for each molecule in the sdf file.
    -----------------------------------------------------------------
    Parameters:
    sdf: str
        path to the smi file
    n_conf: int
        number of 3D conformers to generate
    -----------------------------------------------------------------
    Returns:
    df: pd.DataFrame
        dataframe with "n_conf" x 3D conformers for each molecule
    '''
    
    df = pd.DataFrame({'ROMol': [], 'Smile': [], 'Name': []})
    
    with open(smiles) as f:
        for line in f:
            line = line.split()
            
            if len(line) == 1:
                assert Chem.MolFromSmiles(line[0], sanitize=False) is not None, f"Cannot convert {line[0]} to SMILES."
                mol = Chem.MolFromSmiles(line[0])
                df = df.append({'ROMol': mol, 'Smile': line, 'Name': ''}, ignore_index=True)
            if len(line) == 2:

                assert Chem.MolFromSmiles(line[0], sanitize=False) is not None, f"Cannot convert {line[0]} to SMILES."
                mol = Chem.MolFromSmiles(line[0])
                df = df.append({'ROMol': mol, 'Smile': line[0], 'Name': line[1]}, ignore_index=True)

    # empty list to store all de dataframes
    dfs = list()

    # iterate over all molecules in sdf
    for index, row in df.iterrows():
        curr_mol_df = pd.DataFrame([row.to_dict()])
        currdf = gen3d( row['ROMol'], n_conf)

        # add the previous properties to the current dataframe
        for index, row2 in currdf.iterrows():
            to_add = pd.concat([pd.DataFrame([row2.to_dict()]), curr_mol_df.loc[:, curr_mol_df.columns != 'ROMol']], axis=1)
            dfs.append(to_add)

    # concat all the molecules in one dataframe
    new_dfs = pd.concat(dfs)
    return new_dfs

def iterSDF(sdf: str, n_conf=1) -> pd.DataFrame:
    '''
    It will generate 3D conformers for each molecule in the sdf file.
    -----------------------------------------------------------------
    Parameters:
    sdf: str
        path to the sdf file
    n_conf: int
        number of 3D conformers to generate
    -----------------------------------------------------------------
    Returns:
    df: pd.DataFrame
        dataframe with "n_conf" x 3D conformers for each molecule
    '''
    
    #from sdf to pd.DataFrame
    sdf = PandasTools.LoadSDF(sdf, molColName='ROMol')
    
    # empty list to store all de dataframes
    dfs = list()

    # iterate over all molecules in sdf
    for index, row in sdf.iterrows():
        curr_mol_df = pd.DataFrame([row.to_dict()])
        currdf = gen3d( row['ROMol'], n_conf)

        # add the previous properties to the current dataframe
        for index, row2 in currdf.iterrows():
            to_add = pd.concat([pd.DataFrame([row2.to_dict()]), curr_mol_df.loc[:, curr_mol_df.columns != 'ROMol']], axis=1)
            dfs.append(to_add)

    # concat all the molecules in one dataframe
    new_dfs = pd.concat(dfs)
    return new_dfs
        




@timeit
def main():


    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''\
            Generate 3D conformers for a molecule
            =====================================
                Usage: $ python3 Gen3dConformers.py (-s <smi> or -sd <sdf>) -o <output> -n <number of conformers>
        '''),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-s", "--smi", action="store", dest="smi")  
    parser.add_argument("-sd", "--sdf", action="store", dest="sdf")
    parser.add_argument("-n", "--nconf", action="store", dest="nconf", default=1)
    parser.add_argument("-o", "--out", action="store", dest="out")
    args = parser.parse_args()


    if not args.smi and not args.sdf:
        parser.print_help()
        sys.exit(1)
    
    else:

        assert ((args.sdf and not args.smi) or (args.smi and not args.sdf)), "Please, enter just ONE input file."

        if args.sdf:
            file = args.sdf
            assert (path_Extension(file) in ('.sdf','.sd')), "Wrong input file type (.sdf, .sd)."
            assert (os.path.exists(file)), "Input file does not exist."

            n_conf = int(args.nconf)
            df = iterSDF(file, n_conf)

            if not args.out:
                out = os.path.basename(args.sdf).split('.')[0]+'_conf.sdf'
            else:
                out = args.out


        elif args.smi:
            file = args.smi
            assert (path_Extension(file) in ('.smi','.smile', '.smiles')), "Wrong input file type (.smi, .smile, .smiles)."
            assert (os.path.exists(file)), "Input file does not exist."

            n_conf = int(args.nconf)
            df = iterSMILES(file, n_conf)

            if not args.out:
                out = os.path.basename(file).split('.')[0]+'_conf.sdf'
            else:
                out = args.out
        
        Chem.PandasTools.WriteSDF(df, out, molColName='ROMol', properties=list(df.columns))


if __name__ == '__main__':
    main()

    


