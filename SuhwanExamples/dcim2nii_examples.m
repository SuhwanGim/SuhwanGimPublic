% help dicm2nii

basedir = '/Users/suhwan/Documents/COIL_COMBINED_TEST/COCOAN_SUHWAN_GIM_20210304_155925_168000';

dcmSource = filenames(fullfile(basedir,'*'));
niiFolder = fullfile(basedir,'HEAD_CNIR_NII');
% for i = 1:length(dcmSource)
%     dicm2nii(dcmSource{i}, niiFolder, 0)
% end
%%
dicm2nii(dcmSource{16}, niiFolder, 0)

%% SIMPLE PREPROC
func_bold_files = filenames(fullfile(niiFolder,'*.nii'));

func_bold_files{1} % CMRR _ MB4
func_bold_files{3} % CMRR _ MB4

i = 3;
[~, ~, ~, ~, outputname] = fmri_mask_thresh_canlab(char(func_bold_files{i}),...
    fullfile(niiFolder, 'implicit_mask.nii'));
implicit_mask_file = outputname;
mask = fmri_data(implicit_mask_file,implicit_mask_file);
dat = fmri_data(func_bold_files{i}, implicit_mask_file);
orthviews(mean(dat),'overlay',implicit_mask_file);
dat = preprocess(dat, 'outliers', 'plot');  % Spike detect and globals by slice
subplot(5, 1, 5);
dat = preprocess(dat, 'outliers_rmssd', 'plot');  % RMSSD Spike detect


%% Comparison
semic_mask = '/Volumes/sein/data_SEMIC/imaging/preprocessed/sub-semic003/implicit_mask.nii';
mdat2 = fmri_data(semic_mask,semic_mask );
