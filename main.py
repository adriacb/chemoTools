from chem.ChemSpace import ChemSpace
from chem.Molecules import Molecules
from chem.utils import *


def main():
    st = time.time()
    iptacopan = Chem.MolFromSmiles("CCOC1CCN(C(C1)C2=CC=C(C=C2)C(=O)O)CC3=C(C=C(C4=C3C=CN4)C)OC")
    ms_iptacopan = Molecules([iptacopan], "IPTACOPAN")

    iptacopan_s = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/iptacopan_scaffold.sdf")
    ms_iptacopan_s = Molecules(iptacopan_s, 'Iptacopan_scaffold')

    enamine_all = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/enamine_all_filtered.sdf")
    ms_enamine_all = Molecules(enamine_all, 'Enamine')

    enamine_in = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/enamine_filtered.sdf")
    ms_enamine_in = Molecules(enamine_in, 'EnamineIn')

    CFB = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/CFB_filtered.sdf")
    ms_CFB = Molecules(CFB, 'CFB')

    InSpace = Chem.SDMolSupplier("/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/InSpace_filtered.sdf")
    ms_InSpace = Molecules(InSpace, 'InSpace')
    

    # Create a TARGET ChemSpace object
    mw = np.arange(150, 450, 1)
    logP = np.arange(0, 3, 0.1)
    nrings = np.arange(0, 4, 1)
    pfi = np.arange(0, 5, 0.1)
    hbd = np.arange(0,3,1)
    
    region = np.array(np.meshgrid(logP, mw, hbd, nrings, pfi)).T.reshape(-1, 5)


    ms_region = Molecules(region)
    #print(ms_region.properties)

    #ms_region
    # CREATE CHEMSPACE OBJECT CONTAINING ALL SETS OF MOLECULES
    cs = ChemSpace([ms_iptacopan, ms_iptacopan_s, ms_enamine_all, ms_enamine_in, ms_CFB, ms_InSpace])
    print(cs, end='\n\n')
    #print(cs.properties)


    coords, model, scaler = cs.PCA(nComponents=2)

    results = cs.to_df()
    print(results.head())

    scaled_region = scaler.transform(ms_region.properties)
    region_coords = model.transform(scaled_region)

    
    # Select the target region
    # region = results[results.setName == ms_region.setName]
    # Select the remaining sets
    dfcopy = results[results.setName != ms_region.setName]



    # IPTACOPAN distance
    ipta = results[results.setName == 'IPTACOPAN']
    d = compute_distance(ipta[['PC1', 'PC2']].to_numpy(), region_coords)
    print(d)
    

    for set in dfcopy['setName'].unique():
        print(f"Current set: {set}")
        curr_group = dfcopy[dfcopy.setName == set]

        # Calculate the Euclidean distance between each molecule in the chemspace using compute_distance
        distances = compute_distance(curr_group[['PC1', 'PC2']].to_numpy(), region_coords)
        results.loc[results.setName == set, 'EuclDist'] = distances


    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    # Plot the distances between the all molecules in the chemspace
    # and the "target" region of the chemspace.
    #plot_distances(results, ref_region=ms_region.setName) #filename='/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/distances.png')
    
   
    fig, ax = plt.subplots(1,1, figsize=(10, 6))
    plt.grid()

    ax1 = ax.scatter(x=region_coords[:,0], y=region_coords[:,1], marker='o', s=125, color="whitesmoke", alpha=1)
    ax2 = ax.scatter(data=results, x='PC1', y='PC2', marker='+', c=results['EuclDist'], cmap='plasma')
    ax3 = plt.scatter(data=results[results.setName == "IPTACOPAN"],
                    x='PC1',
                    y='PC2',
                    marker='X',
                    s=223)

    cbar = fig.colorbar(ax2, ax=ax)
    cbar.set_label('Euclidean Distance')

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    
    # LEGEND
    patch1 = mpatches.Patch(color='whitesmoke', label='Target Chemical Space')

    ipt = mlines.Line2D([], [], color='indigo', marker='X', linestyle='None',
                            markersize=10, label='Iptacopan', markeredgewidth=4)

    plt.legend(handles=[patch1, ipt], frameon=False, loc='upper left')
    plt.tight_layout()

    plt.show()

    # Histogram of distances by set name
    r = results.pivot(columns='setName', values='EuclDist')
    
    for s in results['setName'].unique():
        ax = r[s].plot(kind='hist', bins = 100, 
                             figsize=(12,8),
                             alpha = 0.6, grid=True)

        ax.set_xlabel("Euclidean Distances")
        #leg1 = ax.legend(['CFB', 'Enamine Real', 'Enamine Real - In Space', 'Iptacopan scaffold'], loc='upper right')
        #leg2 = ax.legend(handles=[vline], loc='lower right')
        #ax.add_artist(leg1)
        plt.show()


if __name__ == "__main__":
    main()