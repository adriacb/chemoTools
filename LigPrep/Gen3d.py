
def Gen3dFromSDF(file, output='prepared_ligands.sdf'):
    
    SDF = PandasTools.LoadSDF(file)
    SDF['ROMol'] = SDF['ROMol'].apply(lambda x: Chem.AddHs(x) )
    SDF['ROMol'].apply(lambda x: AllChem.EmbedMolecule(x, randomSeed=0xf00d) )
    SDF['ROMol'].apply(lambda x: AllChem.MMFFOptimizeMolecule(x, maxIters=200) )
    #SDF['ROMol'].apply(lambda x: Chem.rdmolops.SanitizeMol(x) )
    SDF['3DMol'] = SDF['ROMol'].apply(lambda x: Chem.MolFromMolBlock((Chem.MolToMolBlock(x)) ) )
    for i in range( len(SDF)):
        for col, item in SDF.iloc[i].items():
            if col not in ('ROMol', '3DMol'):
                SDF['3DMol'].iloc[i].SetProp(col, item)

    SDF['ROMol'] = SDF['3DMol']
    SDF.drop(columns=['3DMol'], axis=1, inplace=True)
    SDF.drop(columns=['ID'], axis=1, inplace=True)

    Chem.PandasTools.WriteSDF(SDF, output, molColName='ROMol', properties=list(SDF.columns))
    print(f"Output file: {output}")
    #os.system("awk '/[A-Z][\-]/{next} 1' test.sdf > rm_dashes.sdf")



def Gen3dFromSMILES(file):
    pass


