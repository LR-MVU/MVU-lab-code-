/*
  Image J macro to detect all spots above a threshold value and
  perform 3D spot segmentation using Gaussian fit local thresholding
  for IF P-Body z-stack images. 
 To Use:
 - Ensure Conds matches subdirectories to loop through for each condition
 - Set thresholds values determined from 3D_spot_seg procedure.
 - Run program
 - You will be prompted to select the working directory
 
 Output:
 - z-stack of segmentation mask as .tiff file
 - max projecton of segmentation mask as .tiff file
 */

wds = getDir("Choose Directory");
// Enter list of folder names to find image data
Conds = newArray("MEFs WT_100ng_15H_Dox+/", "MEFs WT_100ng_24H_Dox-/", "NBDY_100ng_24H_Dox-/")
//Ctl/", "HS/", "HS_4H_R/", "HS_6H_R/", "HS_8H_R/", "HS_10H_R/", "HS_12H_R/"); 
	//"HS_10H_R/", "HS_12H_R/");
// Enter minimum intensity threshold for each folder
thresholds = newArray(8131, 9158, 9060);

//Loop through each folder
for (j=0; j<Conds.length; j++){
	wd = wds + Conds[j];
	od = wd + "/3dSpotSeg/"; //output directory
	File.makeDirectory(od);

	Files = getFileList(wd);
	
	thresh = thresholds[j];
	print(Conds[j]);
	print(thresh);
	
	// filter for file names containing identifier
	Files640 = newArray();
	for(i=0; i<Files.length; i++){
		if (matches(Files[i], ".*555S.*")){
			Files640=Array.concat(Files640, Files[i]);
		}
	}
	
	// loop through all files
	for (i=0; i<Files640.length; i++){
		open(wd + Files640[i]);
		img = getTitle();
		selectImage(img);
		spotimg = substring(Files640[i],0,lengthOf(Files640[i])-4);
		// run 3D Spot Segmentation macro
		run("3D Spot Segmentation", "seeds_threshold=&thresh local_background=0 local_diff=0 radius_0=1.50 radius_1=3 radius_2=6 weigth=0.50 radius_max=6.00 sd_value=2.00 local_threshold=[Gaussian fit] seg_spot=Block watershed volume_min=3 volume_max=10000 seeds=Automatic spots=img radius_for_seeds=8 output=[Label Image]");
		
		//Save outputs
		selectImage("Index");
		saveAs("Tiff", od + Files640[i]);
		run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", od + "max_" + Files640[i]);
		close("*");
	}
}

