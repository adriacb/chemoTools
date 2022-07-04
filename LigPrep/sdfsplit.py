from time import time
from utils import *

def sdfsplit(sdf: pd.DataFrame, n_splits: int) -> list:
    '''
    Split sdf file into n_splits files.
    -----------------------------------------------------------------
    Parameters:
    sdf: pd.DataFrame
        dataframe with molecules in sdf format
    n_splits: int
        number of splits to create
    -----------------------------------------------------------------
    Returns:
    dfs: list
        list of dataframes with molecules in sdf format
    '''
    sdf = PandasTools.LoadSDF(sdf, molColName='ROMol')


    assert n_splits > 0, "n_splits must be greater than 0"
    assert n_splits <= len(sdf), "n_splits must be less than or equal to the number of molecules in the sdf file"
    assert len(sdf) % n_splits == 0, "n_splits must divide the number of molecules in the sdf file"

    dfs = list()

    # split dataframe into n_splits dataframes
    for i in range(n_splits):
        dfs.append(sdf.iloc[i::n_splits])
 
    return dfs




@timeit
def main():


    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''\
            Split SDF file into n_splits files.
            =====================================
                Usage: $ python3 sdfsplit.py --sdf <sdf> -o <output> -n <number of splits>
        '''),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-i", "--sdf", action="store", dest="sdf")
    parser.add_argument("-n", "--nsplit", action="store", dest="nsplit", default=1)
    parser.add_argument("-o", "--out", action="store", dest="out")
    args = parser.parse_args()


    if not args.sdf:
        parser.print_help()
        sys.exit(1)
    
    else:

        assert (args.sdf), "Please, enter an .sdf input file."

        if args.sdf:
            file = args.sdf
            assert (path_Extension(file) in ('.sdf','.sd')), "Wrong input file type (.sdf, .sd)."
            assert (os.path.exists(file)), "Input file does not exist."

            nsplit = int(args.nsplit)
            df = sdfsplit(file, nsplit)

            for i in range(len(df)):

                if not args.out:
                    out = os.path.basename(args.sdf).split('.')[0]+'_{}_split.sdf'.format(i+1)
                else:
                    out = args.out+'_{}.sdf'.format(i+1)
                    
                Chem.PandasTools.WriteSDF(df[i], out, molColName='ROMol', properties=list(df[i].columns))


if __name__ == '__main__':
    main()

    


