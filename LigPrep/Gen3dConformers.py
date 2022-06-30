from utils import *

def gen3d(mol, n_conf):
    m = Chem.AddHs(mol)   
    df = pd.DataFrame({'3DBlock': [], 'ROMol': []})
    
    for i in range(n_conf):
        AllChem.EmbedMolecule(m)
        block = Chem.MolToMolBlock(m)
        df = df.append({'ROMol': Chem.MolFromMolBlock(block)}, ignore_index=True)
    return df


def iterSDF(sdf: str, n_conf=1) -> pd.DataFrame:

    #assert "Check extension file"
    
    #from sdf to pd.DataFrame
    sdf = PandasTools.LoadSDF(sdf, molColName='ROMol')
    
    # empty list to store all de dataframes
    dfs = list()


    for index, row in tqdm(sdf.iterrows(), total=sdf.shape[0]):
        curr_mol_df = pd.DataFrame([row.to_dict()])

        for i in range(n_conf):
            name = f"conf{i}"
            mol = row['ROMol']
            currdf = gen3d(mol, n_conf)
            curr_mol_df['ROMol'] = currdf['ROMol']
            curr_mol_df['conforer_id'] = name
            dfs.append(curr_mol_df)
    
    new_dfs = pd.concat(dfs)
    return new_dfs
        

# def main():
#     m = Chem.MolFromSmiles('CCO[C@H]1CCN([C@@H](C1)C2=CC=C(C=C2)C(=O)O)CC3=C(C=C(C4=C3C=CN4)C)OC')    
#     df = gen3d(m, n_conf=10)
#     Chem.PandasTools.WriteSDF(df, 'test_conf.sdf', molColName='ROMol', properties=list(df.columns))

def main():
    file = '/chemotargets/research/SITALA/IPTACOPAN/docking_rDock/Chembl/best_unique.sdf'
    dfs = iterSDF(file, n_conf=10)
    Chem.PandasTools.WriteSDF(dfs, 'test_conf.sdf', molColName='ROMol', properties=list(dfs.columns))


if __name__ == '__main__':
    main()
    


