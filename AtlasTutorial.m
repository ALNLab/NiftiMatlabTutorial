%% Atlas Analysis
% 1.Overview and requirements:
% This tutorial teaches you how to manipulate nifti files in matlab. It's
% the first tutorial in a series of 3. What we learn here will be applied
% to create a script that automatically tags all clusters in a nifti image
% by cross-referencing each clusters' overlap with ROIs from an atlas.
%
% Here we will learn to put a 0.5mm MNI space atlas prepared in one fMRI 
% processing software (freesurfer) into 2mm MNI space associated with a 
% different fMRI processing software (FSL). The atlas has already been
% converted from TAL to MNI space for us, but this is trivial to do on
% your own. You can try this using MNI2TAL conversions found in niftiManip 
% and niftiUtilities repositiories. 
%
% The primary issue we will have to work around is the interpolation 
% associated with downsampling images. In the process of solving this
% problem we will learn to: extract voxels from nifti files and compute
% their associated coordinates, perform registration with FLIRT, extract
% specific areas of the brain, and read/write NIFTI files. 
% 
% To follow this tutorial you will need to get:
% i) FSL (at least get flirt; doesn't need to be
% You will also need to get/have the following scripts in your matlab path:
% i) Glasser's (2016) atlas (you can find this on the leadDB website)
% ii) SPM (or just the scripts for loading in NIFTI files; this is
% for working with the atlas)
% ivii) NFITTools by Jimmy Shen (you can find this on file exchange for
% matlab, this is for working with FSL's nifti files; also it's better than
% spm for reading stuff into matlab so you should get it anyway! We just
% need SPM for this specific tutorial)
% 
% 2.Introductory registration step -- registering 0.5mm MNI freesurfer 
% template to 2mm MNI FSL template:
% Conveniently, Glasser's (2016) functional surface space atlas has already 
% been converted to 0.5mm MNI space but we want to get it into 2mm MNI space.
% Further, the voxel offset/positioning of the 0.5mm MNI space atlas file
% doesn't line up with the 0.5 mm MNI standard space brain. To get this 0.5mm
% atlas file into the proper 2mm MNI template space, we will use FLIRT to register 
% the two images together. To use FLIRT we will first need to remove the skull in the 
% 0.5mm MNI freesurfer template file by masking voxels outside the brain with fslmaths. 
% Luckily the atlas comes with a mask of the brain so this will be easy. 
% 
% To do these two steps go to terminal and run the following (you need to 
% specify your own unique file paths here so take note):
% fslmaths mni_icbm152_t1_tal_nlin_asym_09a.nii -mas mni_icbm152_t1_tal_nlin_asym_09a_mask.nii mni_icbm152_t1_tal_nlin_asym_09a_brain.nii.gz
% flirt -in mni_icbm152_t1_tal_nlin_asym_09a_brain.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -omat mni_icbm152_t1_TO_MNI_2mm.mat
%
% 3. Why we can't just apply the transformation matrix to the atlas
% You may be thinking that since we've registered the 0.5mm MNI file 
% to a 2mm template, now we can just apply the matrix to the atlas file 
% (HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii) and be done with this tutorial. 
% Unfortunately it's not that simple. Since we are downsampling, voxel 
% intensities will be interpolated. This means our nice whole integer 
% indices for each ROI that we see in the atlas file now will be impossible
% to make out from each other! Go ahead and try to do this with applyxfm, 
% I'll wait. 
%
% To do this run the following in terminal:
% flirt -in HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -applyxfm -init mni_icbm152_t1_TO_MNI_2mm.mat -o mni_icbm152_t1_tal_nlin_asym_09a_brain_2mm_MNI.nii.gz
% 
% Now use fslview to visualize the map you've just created. You should be 
% seeing a bunch of decimals that are very different for almost every single 
% voxel. 
%
% See, not so simple! You can't tell which ROI is which. 
% Rounding wouldn't help because the ROIs bleed into each 
% other once we downsample. To get around this we can go
% back to the original atlas file, and first split up each ROI into a
% unique file. We can then convert each of those ROIs and treat any
% non-zero voxel in the resulting image as a voxel corresponding to that
% ROI. Now this approach is a little bit messy because, again, multiple ROIs
% in a higher res image can load onto the same lower res voxels. If an
% overlap between ROIs makes you nervous, you can pick some threshold for the 
% interpolated intensities that makes for less overall overlap between ROIs. 
% You'll have to play around with this threshold to find one you feel 
% comfortable with. I recommend using a percentile; it is common to threshold
% out the bottom 20% when intensities reflect overlap across subjects in a voxel
% (just for reference). 
%
% 4. Now lets get into the matlab tutorial! This step is all about seperating
% our one atlas file into many files. Each file will correspond to a different
% ROI in the atlas.
%
% First, make sure all of the necessary scripts are in your matlab path. For
% example, if I can't call SPM, I would run the following command in matlab to
% add all subfolders in a folder to my path: 
addpath(genpath('/Users/ateghipc/Projects/scripts/spm12'))

% Now lets also setup file paths to files we want to manipulate
baseDir = '/Users/ateghipc/Projects/atlases/Glasser2016/'; % this is the directory of our project, which will contain all of the atlas files we're interested in
atlasFile = [baseDir '/HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii']; % within the project directory, here is the path to our atlas file which we'll need to extract ROIs from
keyFile = [baseDir '/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt']; % within the project directory, here is the path to the key associated with the atlas file, which tells us what the ROI names are. This will make it easy to name ROI files.
writeDir = [baseDir '/IndividualROIs']; % We should also pick a directory to write our ROI files to 
mkdir(writeDir) % and we should make sure this directory exists!

% Okay now lets start by reading in the key we defined above since
% we'll need it to name our ROIs.
fileID = fopen(keyFile,'r'); 
formatSpec = '%d %s'; % the key has a number refering to the intensity value corresponding to the ROI (I will refer to these values as ROI indices from now on), followed by a string referring to the ROI name
key = textscan(fileID,formatSpec);

keyNum = key{1}; % these are the numbers  we just read in(roi indices)
keyName = key{2}; % these are the strings we just read in(roi names)

% Next, we need to actually load in our atlas file
atlasNii = spm_vol(atlasFile); % this loads in the header of the file
atlas = spm_read_vols(atlasNii); % this loads in the image (i.e. intensities) from the file
atlas = round(atlas); % the intensities could be floats. If we don't round, we won't be able to find exact ROI intensities

% Now that we have our atlas image loaded in, lets find the voxels for ea. ROI and write them out as seperate files!
for roiNum = 1:max(max(max(max(atlas)))) % loop over every roi index; there are 3 dimensions so we need to get the max in each dimension to get the overall max roi index
    outMat = zeros(size(atlas)); % create an empty image based on the size of the atlas image
    idx = find(atlas == roiNum); % find the atlas voxels that correspond to the roi we are working with on this loop
    outMat(idx) = 1; % make those voxels equal to 1 in the empty image we just made
    roiName = keyName{roiNum}; % get the ROI name associated with this roi index
    atlasNii.fname = [writeDir '/' roiName '.nii']; % give the new image file the ROI name and specify directory to save it in
    spm_write_vol(atlasNii,outMat); % save the new image
end

% 5. Converting ROI files into 2mm MNI space
% Great! Now we need to apply the transformation we got earlier to each
% of the files we've generated. You should do this in FSL by running the following
% commands. But first, here are some notes on what this code is assuming. If you 
% didn't alter the filenames in step 2, you should only need to adjust your baseDir
% to match the one you used above.
%
% Assumptions:
% i) you put your transformation matrix into the baseDir directory
% as defined above
% ii) you named the transformation file
% mni_icbm152_t1_TO_MNI_2mm.mat 
% iii) you want the transformed files to be
% appended with the string, '_MNI_2mm'.

% baseDir=/Users/ateghipc/Projects/atlases/Glasser2016/
% cd $baseDir
% for files in *; do flirt -in ${files} -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -applyxfm -init $baseDir/mni_icbm152_t1_TO_MNI_2mm.mat -o ${files%????}_MNI_2mm.nii; done
% for files in *MNI_2mm.nii.gz; do gunzip ${files}; done 

% 6. Binarize ROIs and split left from right hemispheres
% Now that we have our ROIs in 2mm space, we need to binarize them and 
% split them into left and right hemisphere images (the ROIs in the atlas
% are bilateral by default). To do this we can take everything to the left 
% and right of coordinate 44 in the x plane (i.e. saggital slice). Coordinate
% 44 bisects the brain perfectly into two hemisphere --check it out for 
% yourself by looking at the 2MM mni brain in FSL. If you are working in some 
% other space you will need to adjust this identity. Briefly, and as a 
% warning, these coordinate are in system space, they are not true mni coordinates. 
% True MNI coordinates have negative and positive values. We won't worry about that
% for now. I'll wait for you to take a look at that coordinate in fslview.
%
% See? Bisects the hemispheres perfectly. You may notice that I use an x-coordinate
% of 45 instead of 44 below. This is because coordinates in FSL start with 0, but
% in matlab they can only start with 1. Voxels are not coordinates! (see below)
%
% Okay, so we've been working with voxels but we need to get the coordinates
% associated with those voxels. Turns out these coordinates are just the indices
% of the 3 dimensions our coordinate is in. So, to get our coordinates we can 
% just convert the linear indices that are our images into subscripts.

% First, lets get all of the files we've just generated
roiFiles = dir([writeDir '/*MNI_2mm.nii']); 

% Now lets load up the 2mm brain we've registered all of our files to. We
% need this to figure out which voxels are in the left hemsiphere and which
% are in the right. 
template = load_untouch_nii('/Users/ateghipc/Desktop/mPMC/Spt/ROI/ALE/MNI152_T1_2mm_brain.nii');

% From this template lets convert all linear indices to subscripts. On
% second thought, lets only treat voxels that are not zeros so that we are
% only converting voxels in the brain and not empty space.
[coord(:,1),coord(:,2),coord(:,3)] = ind2sub(size(tempMNIFile.img),find(tempMNIFile.img ~= 0)); 
lhVox = find(coord(:,1) >= 45); % from these nonzero voxel coordinates we can get everything in the left hemisphere by looking at the center of the brain on saggital slices (i.e. x-plane)
rhVox = find(coord(:,1) <= 45); % and right hemisphere  

% Okay so we have our left and right hemsiphere voxel coordinates, but since
% our images are in linear indices by default, lets get these coordinates
% back into linear indices. Why? Well if we don't, we will have to convert all 
% of our ROI images to subscripts in order to figure out which of our ROI voxels
% are in the left hemsiphere or the right hemisphere. That sounds terrible.
lhIdx = sub2ind(size(tempMNIFile.img),vox(lhVox,1),vox(lhVox,2),vox(lhVox,3)); 
rhIdx = sub2ind(size(tempMNIFile.img),vox(rhVox,1),vox(rhVox,2),vox(rhVox,3)); 

% And you know, since we already have this MNI brain here, lets use that to
% write out our new files. To use it in this way, we need to change some infof
% in the header of the file since these images won't be anatomical.
template.hdr.dime.datatype = 16;
template.hdr.bitpix = 32;
template.hdr.dime.cal_max = 1;
template.hdr.dime.cal_min = 0;

% Finally lets load in our ROI images and split the hemispheres!
for i = 1:length(roiFiles) % for each roi file
    roiNifti = load_untouch_nii([writeDir '/' roiFiles(i).name]); % load up the file
    roiIdx = find(roiNifti.img ~= 0); % find the indices that correspond to the ROI
    lhRoiIdx = lhIdx(ismember(lhIdx,roiIdx,'rows')); % find out which of those indices are in the list of left hemsiphere indices in the brain
    rhRoiIdx = rhIdx(ismember(rhIdx,roiIdx,'rows')); % find out which of those indices are in the list of right hemsiphere indices in the brain
    
    outMat = zeros(size(roiNifti.img)); % now lets make an empty brain image in the same space
    outMat(lhRoiIdx) = 1; % and lets put an intensity of 1 into each voxel in that brain
    template.img = outMat; % write our new image into the template file 
    save_untouch_nii(template,[writeDir '/' roiFiles(i).name(1:end-4) '_BIN_LH.nii']) % save the template file
    
    outMat = zeros(size(roiNifti.img)); % repeat for the right hemisphere
    outMat(rhRoiIdx) = 1;
    template.img = outMat;
    save_untouch_nii(template,[writeDir '/' roiFiles(i).name(1:end-4) '_BIN_RH.nii'])
end

% You did it! Now go learn to get clusters from your nifti maps using the native matlab image processing toolbox function 
% bwconncomp. You have all the tools you need. Once you've done this you can combine what you've learned to automatically tag
% all clusters in your nifti maps.
