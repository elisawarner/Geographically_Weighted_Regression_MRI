% Shariq Mohammed <shariqm@bu.edu>
% Last edited: Feb 7, 2022

clear all;
close all;

% set working directory
directory = 'C:/Users/shariqm/Documents/Fall 2019/GWR-LGG/LGG-Spyros/WS_MC'
cd(directory);
% set directory to save final images for largest tumor connected component (CC)
target_dir = 'C:\Users\shariqm\Documents\Fall 2019\GWR-LGG\LargestConnected\MC\';

% List the files in the directory
wd_cont = dir();
% Remove the first two and last files
ind = 3:1:size(wd_cont,1);
modality = {'flair','t1','t1Gd','t2','tumor_mask'};

% 12:10'TCGA-CS-6188'; 71:69 'TCGA-EZ-7265A' %issue with whitestripe
ind([10 69]) = []; 

for i = ind
    wd_temp = dir(wd_cont(i).name);
    
    % Set file path for the mask and all four modalities
    for j=1:5
        assignin('base',[modality{j},'_file'],...
            strcat(wd_cont(i).name,'/',[modality{j},'.nii.gz']));
        assignin('base', ['tmr_',modality{j}], niftiread(eval([modality{j},'_file'])));
    end
    
    % number of axial slices
    nslice = size(tmr_tumor_mask,3);
    % initialize a vector to store largest tumor CC size in each axial slice
    largestblobsize = [];
    for l = 1:nslice
        % identify tumor region in slice l
        BW = tmr_tumor_mask(:,:,l)~=0;
        % find the largest tumor CC in slice l
        largestBW = bwareafilt(BW, 1);
        % compute largest tumor CC size in slice l
        largestblobsize(l) = sum(sum(largestBW));
    end
    
    % find the slice number with the largest tumor CC size
    l_final = find(largestblobsize == max(largestblobsize));
    % more than one slice with the same tumor CC size then choose the first one
    if size(l_final,2)>1
        l_final = l_final(1);
    end
    % identify the non tumor region in the slice with largest tumor CC
    BW = tmr_tumor_mask(:,:,l_final)~=0;
    % find the pixels of largest tumor CC in the slice with largest tumor CC
    largestBW = bwareafilt(BW, 1);
    
    % identify the slice across all sequences and mask with the largest tumor CC
    t1 = tmr_t1(:,:,l_final);
    t1Gd = tmr_t1Gd(:,:,l_final);
    t2 = tmr_t2(:,:,l_final);
    flair = tmr_flair(:,:,l_final);
    tumor_mask = tmr_tumor_mask(:,:,l_final);
    
    % assign all pixels not belonging to the largest tumor CC as NaN
    t1(~largestBW) = NaN;
    t1Gd(~largestBW) = NaN;
    t2(~largestBW) = NaN;
    flair(~largestBW) = NaN;
    tumor_mask(~largestBW) = NaN;
    
    % save final images
    save([target_dir wd_cont(i).name '.mat'],...
        't1', 't1Gd', 't2', 'flair', 'tumor_mask');
end