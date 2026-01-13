/*
 Image J macro to detect all spots above a threshold value for
 all images meeting a specifc naming criteria within the selected
 directory. 
 To Use:
 - Ensure Conds matches subdirectories to loop through for each condition
 - Set thresholds to a reasonably low value (If too low program will take 
 	quite long to perform)
 - Run program
 - You will be prompted to select the working directory
 
 Output:
 - .csv file containing spot position and intensity
*/

wds = getDir("Choose Directory");
// Enter list of folder names to find image data
Conds = newArray("CT/", "HS_24H_R/"');
//, "HS+ActD_3HR/","HS+ActD_6HR/", "HS+ActD_9HR/", "HS+ActD_12HR/", "HS+ActD_HS/" 
//"ActD+HS_3HR/", "ActD+HS_6HR/", "ActD+HS_9HR/", "ActD+HS_12HR/", "ActD+HS_HS/", 
//("Ctl/", "HS/", "HS_2H_R/", "HS_3H_R/", "HS_4H_R/", "HS_6H_R/", "HS_8H_R/", "HS_10H_R/", "HS_12H_R/");
// Enter minimum intensity threshold for each folder

//Loop through each folder
for (j=0; j<Conds.length; j++){
	wd = wds + Conds[j];
	od = wd + "/SpotDetect/"; //output directory name
	File.makeDirectory(od)
	
	Files = getFileList(wd);
	print(Conds[j]);
	
	// filter for file names with containing identifier
	Files740 = newArray();
	for(i=0; i<Files.length; i++){
		if (matches(Files[i], ".*740S.*")){ 
			Files740=Array.concat(Files740, Files[i]);
		}
	}
	
	// Loop through i files. (5 should be sufficient)
	for (i=0; i<Files740.length; i++){
		if (matches(Files740[i], ".*740S.*")){
			open(wd + Files740[i]);
			img = getTitle();
			selectImage(img);
			//Run peak finder plugin
			run("3D Maxima Finder", "minimmum=&thresh radiusxy=8 radiusz=2 noise=&thresh");
			//Save output
			saveAs("Results", od + substring(Files740[i],0,lengthOf(Files740[i])-4) + "PB_spots.csv");
			close("Results");
			close("*");
		}
	}
}

