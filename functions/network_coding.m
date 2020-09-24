function out_obj = network_coding(fmri_obj)
% 
%
% :: Usage
%       out_obj = network_coding(fmri_obj)
%
% :: Input
%       - fmri_obj: fmri_data object you want to compare
%
% :: Output
%       - out_obj.get_wh_images(index; see below)
%          
%               ----------------------------------------------------------
%               1. 'Visual network'          2. 'Somato sensory ',
%               3. 'doral Attention',        4. 'Ventral Attention',
%               5. 'Limbic',                 6. 'Fronto Parietal',
%               7. 'Default mode',           8. 'Thalamus',
%               9. 'Hippocampus/Amygdala',  10. 'Brainstem'
%               ----------------------------------------------------------
%
% :: Examples 
%
% mask = which('gray_matter_mask.nii');
% dat = fmri_data(nii_fullpath, mask);
% out_obj = network_coding(fmri_obj);
% othviews(out_obj.get_wh_image(1)); % visual (see index)
%
% See also, radialplot_network_overlap
%   
% Suhwan Gim (suhwan.gim.psych@gmail.com)
% 2020. 05. 12

%% LOAD data
% load(which('buckner_networks_atlas_object.mat'),'atlas_obj');
img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');
img_names = {'visual','somato','dorsal attention','ventral attention', 'Limbic', ...
    'Fronto Parietal','Default','Thalamus','Hippocampus', 'Amygdala','Brainstem'};
%atlas_obj = load_atlas('yeo17networks');
% this image is in CanlabCore/canlab_canonical_brains/Combined_multiatlas_ROI_masks/rBucknerlab_7clusters_SPMAnat_Other_combined.img

% check imgs
if isempty(img)
    warning('No file: rBucknerlab_7clusters_SPMAnat_Other_combined.img');
    error('Please check whether CanlabCore is in your path!');
end
%mask = fmri_data(img, which('brainmask_canlab.nii'));
mask = fmri_data(img, which('gray_matter_mask.nii')); 
%mask = fmri_dat(img, fmri_obj.mask.dat_descrip); % using mask imgs of fmri_obj
%%
dat = [mask.dat==1 | mask.dat==8 | mask.dat==15 ...
    mask.dat==2 | mask.dat==9 | mask.dat==16 ...
    mask.dat==3 | mask.dat==10 | mask.dat==17 ...
    mask.dat==4 | mask.dat==11 | mask.dat==18 ...
    mask.dat==5 | mask.dat==12 | mask.dat==19 ...
    mask.dat==6 | mask.dat==13 | mask.dat==20 ...
    mask.dat==7 | mask.dat==14 | mask.dat==21 ...
    mask.dat>=22 & mask.dat<=35 ...
    mask.dat>=36 & mask.dat<=47 ... % including 48 (hypothalamus)
    mask.dat==49];
%% 
for index = 1:10 
    mask.dat = double( dat(:,index));
    conj_dat(:,index) = double(mask.dat >0 & fmri_obj.dat ~=0);
end
%orthviews_multiple_objs({mask fmri_obj});

%% output
res_mask = fmri_data(which('gray_matter_mask.nii'), which('gray_matter_mask.nii'));
out_obj = res_mask;
out_obj.dat = conj_dat;

end






