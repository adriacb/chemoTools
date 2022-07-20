from chem.ChemSpace import ChemSpace
from chem.Molecules import Molecules
from chem.utils import *
import operator
# import Kmeans, DBScan
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN

def main2():
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
    
   
    # fig, ax = plt.subplots(1,1, figsize=(10, 6))
    # plt.grid()

    # ax1 = ax.scatter(x=region_coords[:,0], y=region_coords[:,1], marker='o', s=125, color="whitesmoke", alpha=1)
    # ax2 = ax.scatter(data=results, x='PC1', y='PC2', marker='+', c=results['EuclDist'], cmap='plasma')
    # ax3 = plt.scatter(data=results[results.setName == "IPTACOPAN"],
    #                 x='PC1',
    #                 y='PC2',
    #                 marker='X',
    #                 s=223)

    # cbar = fig.colorbar(ax2, ax=ax)
    # cbar.set_label('Euclidean Distance')

    # ax.set_xlabel("PC1")
    # ax.set_ylabel("PC2")
    
    # # LEGEND
    # patch1 = mpatches.Patch(color='whitesmoke', label='Target Chemical Space')

    # ipt = mlines.Line2D([], [], color='indigo', marker='X', linestyle='None',
    #                         markersize=10, label='Iptacopan', markeredgewidth=4)

    # plt.legend(handles=[patch1, ipt], frameon=False, loc='upper left')
    # plt.tight_layout()

    # plt.show()

    # # Histogram of distances by set name
    # r = results.pivot(columns='setName', values='EuclDist')
    
    # for s in results['setName'].unique():
    #     ax = r[s].plot(kind='hist', bins = 100, 
    #                          figsize=(12,8),
    #                          alpha = 0.6, grid=True)

    #     ax.set_xlabel("Euclidean Distances")
    #     #leg1 = ax.legend(['CFB', 'Enamine Real', 'Enamine Real - In Space', 'Iptacopan scaffold'], loc='upper right')
    #     #leg2 = ax.legend(handles=[vline], loc='lower right')
    #     #ax.add_artist(leg1)
    #     plt.show()


    # apply K-means to the results
    kmeans = KMeans(n_clusters=4, random_state=0).fit(results[['PC1', 'PC2']].to_numpy())
    results['cluster'] = kmeans.labels_
    print(results.head())

    # apply DBScan to the results
    dbscan = DBSCAN(eps=0.5, min_samples=10).fit(results[['PC1', 'PC2']].to_numpy())
    results['cluster2'] = dbscan.labels_

    fig, ax = plt.subplots(1,1, figsize=(10, 6))
    plt.grid()

    ax1 = ax.scatter(x=region_coords[:,0], y=region_coords[:,1], marker='o', s=125, color="whitesmoke", alpha=1)
    ax2 = ax.scatter(data=results, x='PC1', y='PC2', marker='+', c=results['cluster'], cmap='plasma')
    ax3 = plt.scatter(data=results[results.setName == "IPTACOPAN"],
                    x='PC1',
                    y='PC2',
                    marker='X',
                    s=223)
    plt.show()

    fig, ax = plt.subplots(1,1, figsize=(10, 6))
    plt.grid()

    ax1 = ax.scatter(x=region_coords[:,0], y=region_coords[:,1], marker='o', s=125, color="whitesmoke", alpha=1)
    ax2 = ax.scatter(data=results, x='PC1', y='PC2', marker='+', c=results['cluster2'], cmap='plasma')
    ax3 = plt.scatter(data=results[results.setName == "IPTACOPAN"],
                    x='PC1',
                    y='PC2',
                    marker='X',
                    s=223)
    plt.show()    





# Check Properties
def checkprops(row):    
    
    # Add properties if not present
    if "MolWt" not in row.index:
        row["MolWt"] = Chem.Descriptors.MolWt(row["ROMol"])
    
    if "LogP" not in row.index:
        row["LogP"] = Descriptors.MolLogP(row["ROMol"])
        
    if "NumHDonors" not in row.index:
        row["NumHDonors"] = Chem.Lipinski.NumHDonors(row["ROMol"])
        
    if "NumAromaticRings" not in row.index:
        row["NumAromaticRings"] = Chem.Lipinski.NumAromaticRings(row["ROMol"])
    
    if "PFI" not in row.index:
        row["PFI"] = row["LogP"] + row["NumAromaticRings"]
        
    
    props = {
        "MolWt": ("<", 450),
        "LogP":  ("<", 3.),
        "NumHDonors": ("<=", 2),
        "NumAromaticRings": ("<=", 3),
        "PFI": ("<=", 5.),
    }
    
    operator_name =  {
        "<": operator.lt,
        "<=": operator.le
    }
    
    failed = []
    for prop in props:
        if not operator_name[props[prop][0]](row[prop], props[prop][1]):
            failed.append(prop)
        
    row.loc["FailedProperties"] = failed
    return row
 





def main():
    cfb = PandasTools.LoadSDF("/chemotargets/research/SITALA/IPTACOPAN/03.Known/clarity.sdf")
    cfb['LogP'] = cfb.apply(lambda row: Descriptors.MolLogP(row["ROMol"]), axis=1)
    cfb["MolWt"] = cfb.apply(lambda row: Chem.Descriptors.MolWt(row["ROMol"]), axis=1)
    cfb['NumHDonors'] = cfb.apply(lambda row: Chem.Lipinski.NumHDonors(row["ROMol"]), axis=1)
    cfb['NumAromaticRings'] = cfb.apply(lambda row: Chem.Lipinski.NumAromaticRings(row["ROMol"]), axis=1)
    cfb['PFI'] = cfb['LogP'] + cfb['NumAromaticRings']
    cfb = cfb.apply(checkprops, axis=1)
    cfb["iptacopan"] = cfb["inchikey"] == "RENRQMCACQEWFC-UGKGYDQZSA-N"
    cfb[cfb["iptacopan"]]
    cfb["peptide"] = cfb["peptide"] == "True"
    cfb["iptacopan_scaffold"] = cfb["iptacopan_scaffold"] == "True"
    descriptors = cfb[['MolWt', 'LogP', 'NumHDonors', 'NumAromaticRings', 'PFI']].values

    # Scale values
    scaler = StandardScaler()
    descriptors_scaled = scaler.fit_transform(descriptors)

    # PCA
    pca = PCA()
    descriptors_pca = pca.fit_transform(descriptors_scaled)

    print(sum(pca.explained_variance_ratio_[:2]))

    mw = np.arange(150, 450, 1)
    logP = np.arange(0, 3, 0.1)
    nrings = np.arange(0, 4, 1)
    pfi = np.arange(0, 5, 0.1)
    hbd = np.arange(0,3,1)
    
    region = np.array(np.meshgrid(mw, logP, hbd, nrings, pfi)).T.reshape(-1, 5)
    ms_region = Molecules(region)
    region = pca.transform(scaler.transform(ms_region.properties))
    region.shape

    outliers2, = np.where((descriptors_pca[:,1] > 4) == True)
    outliers1, = np.where((descriptors_pca[:,0] > 6) == True)
    outliers = np.unique(np.concatenate([outliers1, outliers2]))
    cfb["outlier"] = False
    cfb.loc[outliers, "outlier"] = True
    cfb[cfb["outlier"]]
    cfb["activity_ct_numeric"] = pd.to_numeric(cfb["activity_ct_numeric"])
    df = cfb.loc[ (~cfb["outlier"]) & (cfb["activity_ct_numeric"] > 5)]
    print(df.shape)
    descriptors = df[['MolWt', 'LogP', 'NumHDonors', 'NumAromaticRings', 'PFI']].values

    # Scale values
    scaler = StandardScaler()
    descriptors_scaled = scaler.fit_transform(descriptors)

    # PCA
    pca = PCA()
    descriptors_pca = pca.fit_transform(descriptors_scaled)

    print(sum(pca.explained_variance_ratio_[:2])) 


    #In Space
    ms_inspace = Molecules(descriptors[df["FailedProperties"].apply(len) == 0], "InSpace")
    
    # Iptacopan
    ms_iptacopan = Molecules(descriptors[df["iptacopan"]], "Iptacopan")

    # iptacopan_scaffold
    ms_iptacopan_scaffold = Molecules(descriptors[df["iptacopan_scaffold"]], "Iptacopan_Scaffold")

    # CFB
    ms_cfb = Molecules(descriptors, "CFB")


    cs = ChemSpace([ms_inspace, ms_iptacopan, ms_iptacopan_scaffold, ms_cfb])
    print(cs, end='\n\n')

    coords, model, scaler = cs.PCA(nComponents=2)
    results = cs.to_df2()
    print(results.head())
    print(results.shape)

    # scaled_region = scaler.transform(ms_region.properties)
    # region_coords = model.transform(scaled_region)

    # dfcopy = results[results.setName != ms_region.setName]

    # # IPTACOPAN distance
    # ipta = results[results.setName == 'Iptacopan']
    # d = compute_distance(ipta[['PC1', 'PC2']].to_numpy(), region_coords)
    # print(d)
    

    # for set in dfcopy['setName'].unique():
    #     print(f"Current set: {set}")
    #     curr_group = dfcopy[dfcopy.setName == set]

    #     # Calculate the Euclidean distance between each molecule in the chemspace using compute_distance
    #     distances = compute_distance(curr_group[['PC1', 'PC2']].to_numpy(), region_coords)
    #     results.loc[results.setName == set, 'EuclDist'] = distances

    # # Plot the distances between the all molecules in the chemspace
    # # and the "target" region of the chemspace.
    # #plot_distances(results, ref_region=ms_region.setName) #filename='/chemotargets/research/SITALA/IPTACOPAN/PCA_Analysis/Distances/distances.png')
    # fig, ax = plt.subplots(1,1, figsize=(10, 6))
    # plt.grid()

    # ax1 = ax.scatter(x=region_coords[:,0], y=region_coords[:,1], marker='o', s=125, color="whitesmoke", alpha=1)
    # ax2 = ax.scatter(data=results, x='PC1', y='PC2', marker='+', c=results['EuclDist'], cmap='plasma')
    # ax3 = plt.scatter(data=results[results.setName == "Iptacopan"],
    #                 x='PC1',
    #                 y='PC2',
    #                 marker='X',
    #                 s=223)

    # cbar = fig.colorbar(ax2, ax=ax)
    # cbar.set_label('Euclidean Distance')

    # ax.set_xlabel("PC1")
    # ax.set_ylabel("PC2")
    
    # # LEGEND
    # patch1 = mpatches.Patch(color='whitesmoke', label='Target Chemical Space')

    # ipt = mlines.Line2D([], [], color='indigo', marker='X', linestyle='None',
    #                         markersize=10, label='Iptacopan', markeredgewidth=4)

    # plt.legend(handles=[patch1, ipt], frameon=False, loc='upper left')
    # plt.tight_layout()

    # plt.show()
    # print(results.pivot(columns='setName', values='EuclDist'))
    # # Histogram of distances by set name
    # ax = results.pivot(columns='setName', values='EuclDist').plot(kind='hist', 
    #                           bins = 100, 
    #                           figsize=(12,8),
    #                           alpha = 0.6, grid=True)
    
    # vline = ax.vlines(x = d[0], ymin = 0, ymax = 80,
    #         colors = 'black',
    #         label = 'Iptacopan')
    # ax.set_xlabel("Euclidean Distances")
    # leg1 = ax.legend(['CFB', 'In Space', 'Iptacopan', 'Iptacopan scaffold'], loc='upper right')
    # leg2 = ax.legend(handles=[vline], loc='lower right')
    # ax.add_artist(leg1)

    # plt.show()


if __name__ == "__main__":
    main()