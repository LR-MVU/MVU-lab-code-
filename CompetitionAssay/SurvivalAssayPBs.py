
"""

"""


import os # allows to parse images from folder
from skimage.io import imread #to read files
from skimage.filters import threshold_otsu #to compute threshold values / autodetect nuc
from skimage.measure import label, regionprops #find and label individual nuc
from skimage.morphology import remove_small_objects, remove_small_holes, closing, disk #needed to avoid fuzzy outlines
from sklearn.mixture import GaussianMixture
from scipy import ndimage
import matplotlib.pyplot as plt #to show images 
import numpy as np #for array operations
import pandas as pd #better manipulation in dataframe
import matplotlib.pyplot as plt
import pandas as pd

# Root directory that contains one subfolder per condition ("DOX- cntl", "DOX- 2h", etc.)
#path_root = input("Enter path to the main experiment folder: ").strip()
path_root = '/Users/george/Desktop/Academic/Lab/DATAanalysis/CompetitionAssay_Shruti'
#check for conditions 
Conditions = [f for f in os.listdir(path_root) if os.path.isdir(os.path.join(path_root, f))]
print("\n1: condition folders found:")
for c in Conditions:
    print("___", c)

# Suffixes for each channel
suffix_dapi  = '_395S'   # DAPI channel
suffix_green = '_470S'   # Green channel
suffix_red   = '_740S'   # Red channel (note that this will be called ch640 in the code)

# Minimum number of nearest-neighbor P-bodies a cell needs to be considered "PB positive"
pb_thresh = 2



# Start analysis loop for each condition
###################################################################################################
summary_data = []  # will store count summary per condition

for C in Conditions:
    print(f"\n2: Processing condition: {C}")
    os.makedirs(f"{C}/Results", exist_ok=True)

    path_input = os.path.join(path_root, C)

    # find first matching file (just use first field of view for manual steps)
    all_files = os.listdir(path_input)
    dapi_files = [f for f in all_files if suffix_dapi in f and f.endswith('.tif')]
    green_files = [f for f in all_files if suffix_green in f and f.endswith('.tif')]
    red_files = [f for f in all_files if suffix_red in f and f.endswith('.tif')]
    
    if not dapi_files or not green_files or not red_files:
        print("Missing channel, skipping condition.")
        continue
    
    dapi_files.sort()            # sort to match file order
    green_files.sort()
    red_files.sort()

    condition_all_cells = []    # process each image set in condition
    
    # Build arrays for later
    N_PBpositive = np.zeros(len(dapi_files))
    N_PBnegative = np.zeros(len(dapi_files))
    filename = []
    
    # Loop through all images in condition folder
    for i in range(len(dapi_files)):
        dapi_path  = os.path.join(path_input, dapi_files[i])
        ch470_path = os.path.join(path_input, green_files[i])
        ch640_path = os.path.join(path_input, red_files[i])

        print(f"3:Analyzing {os.path.basename(dapi_path)}")

        # read images and max project
        dapi  = imread(dapi_path).max(axis=0)
        ch470 = imread(ch470_path).max(axis=0)
        ch640 = imread(ch640_path).max(axis=0)

        # ACTUAL PROCESSING STARTS HERE - SEGMENT NUCLEI
        ###################################################################################################
        #after loading files and max proj, auto segment nuclei
        th = threshold_otsu(dapi)           # compute thresh
        nuc_mask = dapi > th                # mask for nuc that can be sorted 
        # prevent tiny blobs and nuc at edges (super important!)
        nuc_mask = closing(nuc_mask, disk(1))
        nuc_mask = remove_small_holes(nuc_mask, area_threshold=200)
        nuc_mask = remove_small_objects(nuc_mask, min_size=500)
        labeled_nuc = label(nuc_mask)       #label each mask with an integer 

        # plot initial nuc outlines with IDs
        plt.figure(figsize=(7, 7))          #set size coorindates so integer labels work
        plt.imshow(dapi, cmap='gray')
        plt.contour(nuc_mask, colors='r')   #draw nuc outlines
        for r in regionprops(labeled_nuc):  #attach nuc ID at each nuc 
            y, x = r.centroid
            plt.text(x, y, str(r.label), color ='yellow', fontsize=10, ha='center')
        plt.title(f"Nuclei detected in {C} ({os.path.basename(dapi_path)})")
        plt.show()
    
        # MANUAL REMOVE BAD NUCLEI
        ###################################################################################################
        #prompts to list any nuc to discard (there is no loop, you need to select properly the first time)
        #new cleaner nuc collection
        confirmed = False #for comfirming correct nuc loop
        clean_mask = np.copy(nuc_mask) #copy of nuc mask but excludes any that were listed
        while not confirmed: #loop so that you can select additional nuc, but MUST carry whole list of nuc in new line 
            #remove_input = input("4:List nuclei IDs to remove (comma-separated), you can add more after! : ").strip()
            #remove_input = '93, 95, 102, 98, 148, 66'
            remove_input = []
            removed_labels = [int(x) for x in remove_input.split(',')] if remove_input else []
            clean_mask = np.copy(nuc_mask)
            for lbl in removed_labels:
                clean_mask[labeled_nuc == lbl] = False #sets excluded nuc back to background

            # preview of kept nuclei- plot with only selected nuc, can comment out but its fast 
            plt.figure(figsize=(7, 7))
            bright_470 = 7.0   # green (470)
            bright_640 = 7.0  # red (640) - brighter
            bright_dapi = 5.0  # blue (DAPI)
            


        #    happy = input("5:Done removing nuclei? (y/n):")
            
        #    if happy == 'y':
        #        confirmed = True
        #        removed_nuclei = removed_labels
        #    elif happy == 'n':
        #        continue
        #    else:
        #        print("Please type 'y' or 'n'. Repeating selection loop.")
            confirmed = True      
        
        clean_labeled = label(clean_mask)   #uses dapi channel to check eacgh cleaned nuc across other channels
        green_thr = float(input("6:Enter green (470) intensity threshold (~3500): "))
        #green_thr = 3500 # GW: for debugging so it doesn't stop every time (maybe this doesn't matter since detection depends on PBs anyways)
        green_nuc = [r.label for r in regionprops(clean_labeled) if np.mean(ch470[clean_labeled == r.label]) > green_thr]
        green_nuc = np.array(green_nuc)
        
        
        
        # DETECT P-BODIES 
        ###################################################################################################
        # Flatten image to (n_pixels, 1)
        pixels = ch640.reshape(-1, 1)

        gmm = GaussianMixture(n_components=3, random_state=0)
        gmm.fit(pixels)

        # Predict class for each pixel
        pb_classes = gmm.predict(pixels)

        # Reshape back to image
        
        pb_mask = pb_classes.reshape(ch640.shape) == 2
        
        #gmm = GaussianMixture(n_components=3).fit(ch640)
        #pb_mask = gmm.predict(ch640).reshape(ch640.shape) == 2
        pb_mask_label = label(pb_mask)
        filtered_pb_mask = np.zeros_like(pb_mask_label)
        max_area = 30
        props = regionprops(pb_mask_label)
        keep_labels = [r.label for r in props if r.area <= max_area]

        # Mask out large components
        pb_mask_label = label(pb_mask)
        filtered_pb_mask = np.where(np.isin(pb_mask_label, keep_labels), 1, 0)
        filtered_pb_mask = label(filtered_pb_mask)
                
        nuc_binary = clean_labeled > 0
        dist, indices = ndimage.distance_transform_edt(~nuc_binary, return_indices=True)
        
        
        pbs = regionprops(pb_mask_label)

        # Plot detected PBs
        plt.figure(figsize=(7, 7))
        plt.imshow(ch640, cmap='gray')
        plt.contour(filtered_pb_mask > 0, colors='cyan', linewidths=0.6)
        plt.title("P-body detection (740 channel)")
        plt.axis('off')
        plt.show()        


        pb_to_nuc = []
        # Get PB centroids as arrays (instead of using a long for loop)
        pb_centroids = np.array([r.centroid for r in regionprops(filtered_pb_mask)])
        pb_centroids = np.round(pb_centroids).astype(int)
        
        # ensure within image bounds 
        pb_centroids[:, 0] = np.clip(pb_centroids[:, 0], 0, clean_labeled.shape[0] - 1)
        pb_centroids[:, 1] = np.clip(pb_centroids[:, 1], 0, clean_labeled.shape[1] - 1)

        # remove 'orphan' PBs that are too far away from any nuclei         
        pb_dists = dist[pb_centroids[:, 0], pb_centroids[:, 1]]
        #orphanVals = [x for x in pb_dists if x <= 70]
        orphanIdx = [i for i in range(len(pb_dists)) if pb_dists[i]>=70]
        pb_centroids_refn = np.delete(pb_centroids, orphanIdx, axis=0)
        
        # Vectorized lookup of nearest nucleus
        nearest_rows = indices[0, pb_centroids_refn[:, 0], pb_centroids_refn[:, 1]]
        nearest_cols = indices[1, pb_centroids_refn[:, 0], pb_centroids_refn[:, 1]]
        
        pb_to_nuc = clean_labeled[nearest_rows, nearest_cols]

        # number of PBs for each cell 
        pb_to_nuc = np.array(pb_to_nuc)
        # Ignore background (label 0), if present
        pb_to_nuc = pb_to_nuc[pb_to_nuc > 0]
        pb_counts = np.bincount(pb_to_nuc)   
        
        # Build per-cell table
        cell_ids = np.arange(len(pb_counts))
        cell_df = pd.DataFrame({
            'cell_id': cell_ids,
            'pb_count': pb_counts
        })

        cell_df = cell_df[cell_df.cell_id > 0]
        cell_df['PB_positive'] = cell_df.pb_count >= pb_thresh


        ## Plot results 
        pb_counts = np.bincount(pb_to_nuc[pb_to_nuc > 0])
        nuc_labels = np.unique(clean_labeled)
        nuc_labels = nuc_labels[nuc_labels > 0]

        pb_pos_nuclei = nuc_labels[  # cells w/ threshold nearest-neighbor PB spots 
            pb_counts[nuc_labels] > pb_thresh 
            ]
        
        pb_positive_nuclei = np.union1d(pb_pos_nuclei, green_nuc)

        pb_negative_nuclei = nuc_labels[  # cells with confidently no P bodies
            pb_counts[nuc_labels] == 0
            ]
        
        pb_uncertain_nuclei = nuc_labels[  # cells with uncertainty (# of P bodies >0 but <threshold)
            (pb_counts[nuc_labels] <= pb_thresh) &
            (pb_counts[nuc_labels] != 0)
            ]

        pb_pos_mask = np.isin(clean_labeled, pb_positive_nuclei)
        pb_neg_mask = np.isin(clean_labeled, pb_negative_nuclei)
        pb_unc_mask = np.isin(clean_labeled, pb_uncertain_nuclei)


        plt.figure(figsize=(8, 8))



        # Background channels
        nChannels = 3
        rgb = np.zeros((*dapi.shape, nChannels), dtype=float) # Base RGB image

        #(ch640-np.mean(ch640))
        rgb[..., 0] = np.clip(ch640 / ch640.max(), 0, 1)*10 # PBs (might be actually 740)
        rgb[..., 1] = np.clip(ch470 / ch470.max(), 0, 1)*1    # GFP (green)
        rgb[..., 2] = np.clip(dapi  / dapi.max(), 0, 1)     # DAPI (blue)


        plt.imshow(rgb)

        # PB-positive nuclei (yellow outlines)
        plt.contour(pb_pos_mask, colors='yellow', linewidths=0.8, linestyles='--', label='PB-positive nuclei')

        # PB-negative nuclei (blue outlines)
        plt.contour(pb_neg_mask, colors='cyan', linewidths=0.8, linestyles='--', label='PB-negative nuclei')

        # uncertain nuclei (purple outlines)
        plt.contour(pb_unc_mask, colors='xkcd:eggshell', linewidths=1.2, label='uncertain nuclei (counted as PB-negative)')


        # P-body outlines
        plt.contour(filtered_pb_mask > 0, colors='red', linewidths=0.7, label='P-bodies')

        
        plt.axis('off')

        # Legend (manual)
        plt.plot([], [], 'yellow', linestyle='--', label='PB-positive nuclei')
        plt.plot([], [], 'cyan', linestyle='--', label='PB-negative nuclei')
        plt.plot([], [], 'xkcd:eggshell', label='uncertain nuclei (counted as PB-negative)')
        plt.plot([], [], 'red', label='P-bodies')
        plt.legend(loc="lower center", bbox_to_anchor=(0.5, 1.05), ncol=1)

       
        
  
        # SAVE RESULTS (final image and cell count table)
        ########################################################
        filename.append(os.path.splitext(os.path.basename(dapi_path))[0])

        plt.savefig(f"{C}/Results/COUNTEDIMAGE_{filename[i]}.png", dpi=300, bbox_inches="tight")
        
        plt.show()
        
        # Get final cell counts: 
        N_PBpositive[i] = len(pb_positive_nuclei)
        N_PBnegative[i] = len(pb_negative_nuclei) + len(pb_uncertain_nuclei)
        
    save_df = pd.DataFrame({
        'Name': filename,
        'N PBpositive': N_PBpositive,
        'N PBnegative': N_PBnegative
        })
    
    save_df.to_csv(f"{C}/Results/{C}_cellCountsTable.csv", index=False)

        