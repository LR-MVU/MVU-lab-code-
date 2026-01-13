%%%%%%%%%%%%%%%%%%
% Code that takes in tiff stacks in nucleus, mRNA, and IF channels along
% with cell masks from FISH-quant then crops the images to the bounding box
% of the cell mask. 
% Change the paths to where the images are located.
% Keep the cell masks in a folder called 'SkeletonAndPrint' within the
% image directory
% The output folder is call bigFISH_In and will be made automatically.
% Input:
%   Nucleus channel tiff stacks, mRNA channel tiff stacks, PBody channel
%   tiff stacks, cell masks for each FOV.
% Output:
%   Cropped Nucleus, mRNA, and PBody tiff stacks and cell mask for each outlined cell
%
% Author: Lydia Hodgins
%%%%%%%%%%%%%%%%%%

clear; clc;

%%%%% Define Paths %%%%%%%
%Change to your working directory
main_dir = 'Z:\Lilli\Ex 4_sept 22 NB rec smFISH (rep1)';
img_dir = 'Z:\Lilli\Ex 4_sept 22 NB rec smFISH (rep1)';
mask_dir = strcat(main_dir, filesep, 'SkeletonAndPrint');
out_dir = strcat(main_dir, filesep, 'bigFISH_In');
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make output directory if it does not yet exist
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%Change the string to match the correct channels
nuc_imgs = dir(fullfile(img_dir,'*395S.tif'));
mRNA1_imgs = dir(fullfile(img_dir,'*555S.tif'));
mRNA2_imgs = dir(fullfile(img_dir,'*640S.tif'));
pbody_imgs = dir(fullfile(img_dir,'*740S.tif'));

for i = 1:length(pbody_imgs)
   nuc = tiffreadVolume(fullfile(img_dir, nuc_imgs(i).name)); 
   mRNA1 = tiffreadVolume(fullfile(img_dir, mRNA1_imgs(i).name));
   mRNA2 = tiffreadVolume(fullfile(img_dir, mRNA2_imgs(i).name));
   pbody = tiffreadVolume(fullfile(img_dir, pbody_imgs(i).name));
   
   expName = extractBefore(pbody_imgs(i).name, "740S");
   disp(expName)
   cell_masks = dir(fullfile(mask_dir, expName + "*_print*.gif"));
   nuc_masks = dir(fullfile(mask_dir, expName + "*_nuc*.gif"));
   
   
   for j=1:length(cell_masks)
      mask = imread(fullfile(mask_dir, cell_masks(j).name));
      nucmask = imread(fullfile(mask_dir, nuc_masks(j).name));
      stats = regionprops(mask, 'BoundingBox');
      x1 = round(stats.BoundingBox(1));
      y1 = round(stats.BoundingBox(2));
      x2 = x1 + stats.BoundingBox(3) -1;
      y2 = y1 + stats.BoundingBox(4) -1;
      nuc_crop = nuc(y1:y2, x1:x2, :);
      mRNA1_crop = mRNA1(y1:y2, x1:x2, :);
      mRNA2_crop = mRNA2(y1:y2, x1:x2, :);
      pbody_crop = pbody(y1:y2, x1:x2, :);
      mask_crop = mask(y1:y2, x1:x2);
      nucmask_crop = nucmask(y1:y2, x1:x2);

      nuc_name = fullfile(out_dir, extractBefore(nuc_imgs(i).name, ".tif") + "_00" + int2str(j) + ".tif");
      mRNA1_name = fullfile(out_dir, extractBefore(mRNA1_imgs(i).name, ".tif") + "_00" + int2str(j) + ".tif");
      mRNA2_name = fullfile(out_dir, extractBefore(mRNA2_imgs(i).name, ".tif") + "_00" + int2str(j) + ".tif");
      pbody_name = fullfile(out_dir, extractBefore(pbody_imgs(i).name, ".tif") + "_00" + int2str(j) + ".tif");
      mask_name = fullfile(out_dir, expName + "_print_00" + int2str(j) + ".tif");
      nucmask_name = fullfile(out_dir, expName + "_nuc_00" + int2str(j) + ".tif");
        
      imwrite(nuc_crop(:,:,1), nuc_name);
      imwrite(mRNA1_crop(:,:,1), mRNA1_name);
      imwrite(mRNA2_crop(:,:,1), mRNA2_name);
      imwrite(pbody_crop(:,:,1), pbody_name);
      imwrite(mask_crop, mask_name);
      imwrite(nucmask_crop, nucmask_name)

      for x = 2 : 41
        imwrite(nuc_crop(:,:,x),nuc_name, 'WriteMode', "append");
        imwrite(mRNA1_crop(:,:,x),mRNA1_name, 'WriteMode', "append");
        imwrite(mRNA2_crop(:,:,x),mRNA2_name, 'WriteMode', "append");
        imwrite(pbody_crop(:,:,x),pbody_name, 'WriteMode', "append");
      end
      
   end
end