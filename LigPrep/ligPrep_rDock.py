
import sys
  
# setting path
sys.path.append('../scripts')

from utils import *

def prepareLigands(path, ext, output):
    # generate 3D coordinates and save into mol2 (required for prepare_ligand)
    #os.system(f"obabel -i {ext} {path} -o mol2 -O out.mol2 -m -h --gen3d")
    os.system(f"obabel -i {ext} {path} -o sdf -O out.sdf -h --gen3d")
    Gen3dFromSDF("out.sdf", output=output)
    #Gen3dFromSDF(path)
    os.system("rm out.sdf")
    
    # os.system(f"obabel -i sdf test.sdf -o pdb -O out.pdb -m -h --property")
    # mol2_gen3d = listPath(extension='.pdb')
    # # for each .mol2 ligand file prepare the ligand with autodock vina (MGLTools)
    # for file in mol2_gen3d:
    #     os.system(f"./prepare_ligand -l {file}")
    
    # pdbqt_gen3d = listPath(extension='.pdbqt')

    # # the output of prepare_ligand is a list of 
    # for f in pdbqt_gen3d:
    #     os.system(f"obabel -ipdbqt {f} -osdf -O prep_{f}.sdf --property")
    
    # os.system("cat prep_*.sdf > prepared_ligands.sdf")
    # os.system("rm *.pdb *.pdbqt prep_out* test.sdf")
    


def main():



    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--smi", action="store", dest="smi")  
    parser.add_argument("-sd", "--sdf", action="store", dest="sdf")
    parser.add_argument("-o", "--out", action="store", dest="out")
    args = parser.parse_args()

    if args.sdf and args.smi:
        print("Please, enter just ONE input file.")
    
    if args.sdf:

        assert (path_Extension(args.sdf) in ('.sdf','.sd')), "Wrong input file type (.sdf, .sd)."
        start = timer()
        print(f"Start Time: {timedelta(seconds=start)}\n")
        prepareLigands(args.sdf, ext='sdf', output=args.out)

    if args.smi:

        assert (path_Extension(args.smi) in ('.smi','.smiles', '.smile')), "Wrong input file type (.smi, .smiles, .smile)."
        start = timer()
        print(f"Start Time: {timedelta(seconds=start)}\n")
        prepareLigands(args.smi, ext='smi')

    else:
        sys.exit()

    end = timer()
    print(f"\nElapsed time: {timedelta(seconds=end-start)}\n")

if __name__ == "__main__":
    main()