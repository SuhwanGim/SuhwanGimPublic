function [out_obj,dat_obj] = network_coding(fmri_obj_path,mask)
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

fmri_obj = fmri_data(fmri_obj_path, mask);
burkner7 = fmri_data(img, mask);


isdiff = compare_space(fmri_obj, burkner7);

if isdiff == 1 || isdiff == 2 % diff space, not just diff voxels
    % == 3 is ok, diff non-empty voxels
    
    % Both work, but resample_space does not require going back to original
    % images on disk.
    %mask = resample_to_image_space(mask, dat);
    burkner7 = resample_space(burkner7, fmri_obj);
    
    % tor added may 1 - removed voxels was not legal otherwise
    %mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    % resample_space is not *always* returning legal sizes for removed
    % vox? maybe this was updated to be legal
    
    if length(burkner7.removed_voxels) == burkner7.volInfo.nvox
        disp('Warning: resample_space returned illegal length for removed voxels. Fixing...');
        burkner7.removed_voxels = burkner7.removed_voxels(burkner7.volInfo.wh_inmask);
    end
    
end

%mask = fmri_data(img, which('brainmask_canlab.nii'));
%mask = fmri_data(img, which('gray_matter_mask.nii')); 
%mask = fmri_data(img, fmri_obj.mask); % using mask imgs of fmri_obj
%%
dat = [burkner7.dat==1 | burkner7.dat==8 | burkner7.dat==15 ...
    burkner7.dat==2 | burkner7.dat==9 | burkner7.dat==16 ...
    burkner7.dat==3 | burkner7.dat==10 | burkner7.dat==17 ...
    burkner7.dat==4 | burkner7.dat==11 | burkner7.dat==18 ...
    burkner7.dat==5 | burkner7.dat==12 | burkner7.dat==19 ...
    burkner7.dat==6 | burkner7.dat==13 | burkner7.dat==20 ...
    burkner7.dat==7 | burkner7.dat==14 | burkner7.dat==21 ...
    burkner7.dat>=22 & burkner7.dat<=35 ...
    burkner7.dat>=36 & burkner7.dat<=47 ... % including 48 (hypothalamus)
    burkner7.dat==49];
%% 
conj_dat = [];
for index = 1:10 
    burkner7.dat = double( dat(:,index));
    conj_dat(:,index) = double(burkner7.dat >0 & fmri_obj.dat ~=0);
end
%orthviews_multiple_objs({burkner7 fmri_obj});

%% output
dat_obj = fmri_obj; %input
out_obj = dat_obj;
out_obj.dat = conj_dat; %results

end






