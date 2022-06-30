from time import time
from utils import *

def gen3d(mol: rdkit.Chem.rdchem.Mol, n_conf: int) -> pd.DataFrame:
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
    '''
    m = Chem.AddHs(mol)   
    df = pd.DataFrame({'ROMol': []})
    
    for i in range(n_conf):
        AllChem.EmbedMolecule(m)
        #AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        block = Chem.MolToMolBlock(m)
        df = df.append({'ROMol': Chem.MolFromMolBlock(block)}, ignore_index=True)
    return df


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

        for i in tqdm(range(n_conf), total=n_conf):
            name = f"conf{i}"
            mol = row['ROMol']
            currdf = gen3d(mol, n_conf)
            curr_mol_df['ROMol'] = currdf['ROMol']
            curr_mol_df['conformer_id'] = name
            dfs.append(curr_mol_df)
    
    new_dfs = pd.concat(dfs)
    return new_dfs
        


def iterSMILES(smi: str, n_conf=1) -> pd.DataFrame:
    '''
    It will generate 3D conformers for each molecule in the smi file.
    -----------------------------------------------------------------
    Parameters:
    smi: str
        path to the smi file
    n_conf: int
        number of 3D conformers to generate
    -----------------------------------------------------------------
    Returns:
    df: pd.DataFrame
        dataframe with "n_conf" x 3D conformers for each molecule
    '''

    df = pd.DataFrame({'SMILES': [], 'ROMol': [], 'conforer_id': []})

    #iterate over smiles
    for smile in smi:
        #for each smile, generate n_conf conformers
        for i in tqdm(range(n_conf), total=n_conf):
            name = f"conf{i}"
            mol = Chem.MolFromSmiles(smile)
            currdf = gen3d(mol, n_conf)
            df = df.append({'ROMol': currdf['ROMol'], 'conformer_id': name}, ignore_index=True)

        



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
            Chem.PandasTools.WriteSDF(df, out, molColName='ROMol', properties=list(df.columns))

        elif args.smi:
            file = args.smi
            assert (path_Extension(file) in ('.smi','.smile', '.smiles')), "Wrong input file type (.smi, .smile, .smiles)."
            assert (os.path.exists(file)), "Input file does not exist."

            n_conf = int(args.nconf)
            df = iterSDF(file, n_conf)

            if not args.out:
                out = os.path.basename(file).split('.')[0]
            else:
                out = args.out
            Chem.PandasTools.WriteSDF(df, out, molColName='ROMol', properties=list(df.columns))


if __name__ == '__main__':
    main()

    


