
%% SIMPLE PREPROCESS FOR MUTLI-ECHO IMAGES 
clear; close;
%% SET PATH
username = char(java.lang.System.getProperty('user.name'));
project_name = 'sync_SEP';
basedir = ['/Users/' username sprintf('/Dropbox/Projects/%s',project_name)];
addpath(genpath(basedir)); % add path
%%
datdir = fullfile(basedir, 'data','coil_test','20210821_multi_echo_RAs','GEONWOO_ME_TEST','COCOAN_MULTI-ECHO_20210821_134827_358000');
datdir = fullfile(basedir, 'data','coil_test','20210821_multi_echo_RAs','GEONWOO2_ME2','COCOAN_MULTI-ECHO_20210821_151325_111000');
datdir = fullfile(basedir, 'data','coil_test','20210821_multi_echo_RAs','MYUNGEUN2_ME','COCOAN_MULTI-ECHO_20210821_154512_629000');

%datdir = fullfile(basedir, 'data','coil_test','20210803_multi_echo','COCOAN_SUHWAN_GIM_20210803_111531_408000');
%datdir = fullfile(basedir, 'data','coil_test','20210803_multi_echo_JJ','COCOAN_SUHWAN_GIM_20210803_093830_891000');
niiFolder = fullfile(basedir,'data','coil_test','20210821_multi_echo_RAs','GEONWOO_ME_TEST','preprocessed');
niiFolder = fullfile(basedir,'data','coil_test','20210821_multi_echo_RAs','GEONWOO2_ME2','preprocessed');
niiFolder = fullfile(basedir,'data','coil_test','20210821_multi_echo_RAs','MYUNGEUN2_ME','preprocessed');
%niiFolder = fullfile(basedir,'data','coil_test','20210803_multi_echo','preprocessed');
%niiFolder = fullfile(basedir,'data','coil_test','20210803_multi_echo_JJ','preprocessed');

%mkdir(niiFolder)
%niiFolder = fullfile(basedir,'data','coil_test','20210803_multi_echo','results_nii');

dcmSource = filenames(fullfile(datdir,'*'));
%% DICOM 2 NIIFT 
% dcm2niix /Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/20210803_multi_echo/COCOAN_SUHWAN_GIM_20210803_111531_408000/1_2_3ME_HALF_CMRR_2_7ISO_PA_64CH+FLEX_MB4_IA2_0
for i = 1:length(dcmSource)
    %addpath    
    [~,b] = fileparts(dcmSource{i});
    outputdir = fullfile(niiFolder,b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % no compress images     
    %1. should be saved using diary
    %2. recommed to save .mat files using json outputs
end
%% SIMPLE PREPROC
%
mask = []; dat = []; 
spike_covariates =[]; 
func_bold_files = []; 
%func_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','*CMRR*mb4_4D.nii'));
func_bold_files = filenames(fullfile(niiFolder,'*MB*00*','2*.nii'));
func_bold_files(contains(func_bold_files,'SBREF')) = []; 
for i = 1:length(func_bold_files)
    clf;
    [a,b]=fileparts(func_bold_files{i});    
    % implicit_mask
    [~, ~, ~, ~, outputname] = fmri_mask_thresh_canlab(char(func_bold_files{i}),...
        fullfile(a, sprintf('implicit_mask_%02d.nii',i)));
    implicit_mask_file = outputname;
    mask{i}= fmri_data(implicit_mask_file,implicit_mask_file);
    dat{i} = fmri_data(func_bold_files{i}, implicit_mask_file);
    dat{i}.images_per_session = size(dat{i}.dat,2);
    % spike id
    dat{i} = preprocess(dat{i}, 'outliers', 'plot');  % Spike detect and globals by slice
    subplot(5, 1, 5);
    dat{i} = preprocess(dat{i}, 'outliers_rmssd', 'plot');  % RMSSD Spike detect
    sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
    set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
    drawnow;
    
    qcspikefilename = fullfile(a, ['qc_spike_plot_' b '.png']); % Scott added some lines to actually save the spike images
    saveas(gcf,qcspikefilename);
    spike_covariates{i} = dat{i}.covariates;
end
%% Slice time correction 
%tr = 1;                 % will be determined 
%mbf = 4;

% func_json_files = filenames(fullfile(datdir, '*MB*00*','*.json'));
% func_json_files (contains(func_json_files ,'SBREF')) = []; 
for i = 1:length(func_bold_files)    
    slice_timing_job = []; 
    json_read = []; json_file = []; 
    % Read Json file 
    [a,b] = fileparts(func_bold_files{i});        
    json_file = fullfile(a,[b '.json']);
    fid = fopen(json_file);
    raw = fread(fid, inf); 
    str = char(raw');
    fclose(fid); 
    json_read = jsondecode(str);
    
    % set parameter 
    slice_time = json_read.SliceTiming;
    tr = json_read.RepetitionTime;
    mbf = json_read.MultibandAccelerationFactor;
    %% DATA
    slice_timing_job{1}.spm.temporal.st.scans{1} = spm_select('expand', func_bold_files(i)); % individual 4d images in cell str
    %% 1. nslices
    Vfirst_vol = spm_vol([func_bold_files{i} ',1']);
    num_slices = Vfirst_vol(1).dim(3);
    slice_timing_job{1}.spm.temporal.st.nslices = num_slices; % number of slices
    %% 2. tr
    slice_timing_job{1}.spm.temporal.st.tr = tr;   
    %% 3. ta: acquisition time
    slice_timing_job{1}.spm.temporal.st.ta = tr - tr * mbf / num_slices; % if not multi-band, mbf = 1;
    %% 4. so: Slice order
    
    slice_timing_job{1}.spm.temporal.st.so = slice_time;
%     if ~exist('custom_slice_timing', 'var')
%         slice_timing_job{1}.spm.temporal.st.refslice = find(slice_time==0, 1, 'first');
%     else
%         slice_timing_job{1}.spm.temporal.st.refslice = find(slice_time==min(slice_time), 1, 'first');
%     end
    if min(slice_time) >= 0 && max(slice_time) <= tr % Time-based
        slice_timing_job{1}.spm.temporal.st.refslice = min(slice_time);
    end
    slice_timing_job{1}.spm.temporal.st.prefix = 'a';
    
    %% Saving slice time correction job    
    %josb_slice_timing_job = slice_timing_job{1};    
    %% RUN    
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run', slice_timing_job);
    %spm_jobman('interactive', slice_timing_job);
    
end
%% Check the affine matrix (?)
temp_t = [];
for i =1:4
    %temp_t{i} = spm_vol([func_bold_files{i} ',2']); % raw images
    temp_t{i} = spm_vol([afunc_bold_files{i} ',2']); % after slice-timing correction 
    %temp_t{i} = spm_vol([rafunc_bold_files{i} ',2']); % after realinment
end
%% Motion correction 
% The consensus is to do 1) estimate realigment parameter using
% before-slice timing correction images and 2) alignment using after-slice
% timing correction images
%
% By J.J
%
adat = []; 
afunc_bold_files = filenames(fullfile(niiFolder,'*MB*00*','a*.nii'));
afunc_bold_files(contains(afunc_bold_files,'SBREF')) = []; 
temp_files = filenames(fullfile(niiFolder,'*MB*00*','*.nii'));
temp_files(~contains(temp_files,'SBREF')) = []; 
sbref_bold_files = temp_files; temp_files = []; 
for i = 1:length(afunc_bold_files)
    def = spm_get_defaults('realign');
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = def.estimate;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions = def.write;
    
    % change a couple things
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % do not register to mean (twice as long)
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0; % do not mask (will set data to zero at edges!)
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0]; % do not output mean image
    
    data = []; data_all = []; 
    data = afunc_bold_files(i);
    data_all = [sbref_bold_files(i); data];
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = data_all;
    
    %% RUN
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run', {matlabbatch});
    %adat{i} = fmri_data(afunc_bold_files{i},mask{i}.fullpath);
end
%% ========================================================================
%                        Do run TEDANA for ME-EPI images (only one image)
% =========================================================================
combined_bold_files = [];
if sum(contains(afunc_bold_files, '3ME')) > 1
    temp_3ME=afunc_bold_files(contains(afunc_bold_files, '3ME'));
    for i =1:length(temp_3ME)
        [a,b] = fileparts(temp_3ME{i});
        combined_bold_files{i} = fullfile(a,'t2smap_outputs','desc-optcom_bold.nii'); %from T2S_workflow in tedana    
        if ~isfile(combined_bold_files{i}) 
            gunzip(fullfile(a,'t2smap_outputs','desc-optcom_bold.nii.gz')); %from T2S_workflow in tedana    
        end
    end
end
%     ==========================    =================================================
%     Filename                      Content
%     ==========================    =================================================
%     T2starmap.nii.gz              Limited estimated T2* 3D map or 4D timeseries.
%                                   Will be a 3D map if ``fitmode`` is 'all' and a
%                                   4D timeseries if it is 'ts'.
%     S0map.nii.gz                  Limited S0 3D map or 4D timeseries.
%     desc-full_T2starmap.nii.gz    Full T2* map/timeseries. The difference between
%                                   the limited and full maps is that, for voxels
%                                   affected by dropout where only one echo contains
%                                   good data, the full map uses the single echo's
%                                   value while the limited map has a NaN.
%     desc-full_S0map.nii.gz        Full S0 map/timeseries.
%     desc-optcom_bold.nii.gz       Optimally combined timeseries.
%    ==========================    =================================================
%% AFTER TEDANA
% Distortion correction 
%% add fsl path 
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

distort_pa_dat = fullfile(niiFolder,'DISTORTION_CORR_64CH_PA_0002','20210821154513_distortion_corr_64ch_pa_2.nii');
distort_ap_dat = fullfile(niiFolder,'DISTORTION_CORR_64CH_PA_POLARITY_INVERT_TO_AP_0003','20210821154513_distortion_corr_64ch_pa_polarity_invert_to_ap_3.nii');
%dicomheader_files = filenames(fullfile(niiFolder,'*CMRR_2*MB4_00*','dcmHeaders.mat'));
rafunc_bold_files = []; 
rafunc_bold_files = filenames(fullfile(niiFolder,'*MB*00*','ra*.nii'));
rafunc_bold_files(contains(rafunc_bold_files,'SBREF')) = []; 

rafunc_bold_files = unique(combined_bold_files)'; %^ {rafunc_bold_files{1}; combined_bold_files};

epi_enc_dir = 'pa';
distortion_correction_out = fullfile(niiFolder, 'Distortion_combied','dc_combined.nii');
%system(['fslmerge -t ', distortion_correction_out, ' ', distort_ap_dat, ' ', distort_pa_dat]);
system(['fslmerge -t ', distortion_correction_out, ' ', distort_pa_dat, ' ', distort_ap_dat]);
%%
% distort_ap_info = nifti(distort_ap_dat);
% distort_pa_info = nifti(distort_pa_dat);
distort_info = nifti(distort_ap_dat);
distort_num = distort_info.dat.dim(4);

rdtime = [0.0397 0.0198]; % manually exported from Json files .TotalReaduoutTime; Three echo files' are same 
rdtime = 0.0198;
for i = 1:length(rafunc_bold_files)     
    readout_time = rdtime(1);
    [a,~ ] = fileparts(rafunc_bold_files{i});
    dc_param = fullfile(a, ['dc_param_', epi_enc_dir, '.txt']);    
    fileID = fopen(dc_param, 'w');
    
    %distort_param_dat = [repmat([0 -1 0 readout_time], distort_num, 1); repmat([0 1 0 readout_time], distort_num, 1)];
    distort_param_dat = [repmat([0 1 0 readout_time], distort_num, 1); repmat([0 -1 0 readout_time], distort_num, 1)];
    fprintf(fileID, repmat([repmat('%.4f\t', 1, size(distort_param_dat, 2)), '\n'], 1, size(distort_param_dat, 1)), distort_param_dat');
    fclose(fileID);

    % Running topup
    disp('Running topup....');
    topup_out = fullfile(a, 'topup_out');
    topup_fieldout = fullfile(a, 'topup_fieldout');
    topup_unwarped = fullfile(a, 'topup_unwarped');
    topup_config = '/usr/local/fsl/src/topup/flirtsch/b02b0.cnf';
    system(['topup --imain=', distortion_correction_out, ' --datain=', dc_param, ' --config=', topup_config, ' --out=', topup_out, ...
        ' --fout=', topup_fieldout, ' --iout=', topup_unwarped]);
    
    
    system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' distortion_correction_out]);
    system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' topup_unwarped '.nii.gz']);
    
    
%     topup_unwarped_png{1} = fullfile(a, 'topup_unwarped_dir-ap_epi.png');
%     topup_unwarped_png{2} = fullfile(a, 'topup_unwarped_dir-pa_epi.png');
    topup_unwarped_png{1} = fullfile(a, 'topup_unwarped_dir-pa_epi.png');
        topup_unwarped_png{2} = fullfile(a, 'topup_unwarped_dir-ap_epi.png');
    
    for top_i = 1:numel(topup_unwarped_png)
        topup_before_list = cellstr(strcat(distortion_correction_out, ',', num2str([2*top_i-1;2*top_i])));
        topup_after_list = cellstr(strcat([topup_unwarped '.nii'], ',', num2str([2*top_i-1;2*top_i])));
        canlab_preproc_show_montage([topup_before_list; topup_after_list], topup_unwarped_png{top_i});
        drawnow;
    end
    close all;
    
    %% Applying topup on BOLD files        
    input_dat = rafunc_bold_files{i};
    [a, b] = fileparts(input_dat);
    dcr_func_bold_files{i,1} = fullfile(a, ['dc' b '.nii']);
    system(['applytopup --imain=', input_dat, ' --inindex=1 --topup=', topup_out, ' --datain=', dc_param, ...
        ' --method=jac --interp=spline --out=', dcr_func_bold_files{i}]);
    
    % removing spline interpolation neg values by absolute
    system(['fslmaths ', dcr_func_bold_files{i}, ' -abs ', dcr_func_bold_files{i}]);
    
    % unzip
    system(['gzip -d -f ' dcr_func_bold_files{i} '.gz']);
    
    
%     %% Applying topup on SBREF files
%     input_dat = sbref_bold_files{i};
%     [a, b] = fileparts(input_dat);
%     dc_func_sbref_files{i,1} = fullfile(a, ['dc_' b '.nii']);
%     system(['applytopup --imain=', input_dat, ' --inindex=1 --topup=', topup_out, ' --datain=', dc_param, ...
%         ' --method=jac --interp=spline --out=', dc_func_sbref_files{i}]);
%     
%     % removing spline interpolation neg values by absolute
%     system(['fslmaths ', dc_func_sbref_files{i}, ' -abs ', dc_func_sbref_files{i}]);
%     
%     % unzip
%     system(['gzip -d -f ' dc_func_sbref_files{i} '.gz']);
%     % system(['gzip -d ' PREPROC.dc_func_sbref_files{i} '.gz']);
    
    %% save mean image across all runs
    dat = fmri_data(char(dcr_func_bold_files{i}), implicit_mask_file);
    mdat = mean(dat);
    [a, b] = fileparts(dcr_func_bold_files{i});
    mdat.fullpath = fullfile(a, ['mean_' b '.nii']);
    try
        write(mdat);
    catch
        write(mdat, 'overwrite');
    end
    %% save mean_r_func_bold_png
    
%     mean_dcr_func_bold_png = fullfile(a, 'mean_dcr_func_bold.png'); % Scott added some lines to actually save the spike images
%     canlab_preproc_show_montage( mdat.fullpath , mean_dcr_func_bold_png);
%     drawnow;
    
    dat = [];
end
%%
%copyfile(SOURCE,DESTINATION,MODE)
dcr_opt_func_bold_files = []; 
for i =1:length(dcr_func_bold_files)
    [a,~] = fileparts(dcr_func_bold_files{i});
    d = fileparts(a);
    dcr_opt_func_bold_files{i} = fullfile(d, sprintf('dcr_OPT_RUN_%02d_combined.nii',i));
    copyfile(dcr_func_bold_files{i},fullfile(d, sprintf('dcr_OPT_RUN_%02d_combined.nii',i)));
    
end
%% Coregi
use_dc = true;
%use_sbref = true;
dcrafunc_bold_files = []; 
%dcrafunc_sbref_files = filenames(fullfile(niiFolder,'*CMRR_2*','dc_*CMRR*mb4_SBRef_4D.nii'));
t1_files = fullfile(niiFolder, 'T1_MPRAGE_SAG_0_7ISO_0008','20210821134827_T1_mprage_sag_0.7iso_8.nii');
t1_files = '/Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/20210821_multi_echo_RAs/GEONWOO_ME_TEST/preprocessed/T1_MPRAGE_SAG_0_7ISO_0008/20210821134827_T1_mprage_sag_0.7iso_8.nii';
t1_files  = '/Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/20210726_ME/results_nii/T1_MPRAGE_SAG_0_7ISO_0012/T1_mprage_sag_0_7iso.nii';
%t1_files = '/Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/JIHOON2_JIHOON2/results_nii/T1_MPRAGE_SAG_0_7ISO_0008/T1_mprage_sag_0_7iso.nii'; % 64-ch full coil
for i = 1:length(rafunc_bold_files)
%     [a,b,c] = fileparts(rafunc_bold_files{i});
%     dcrafunc_bold_files{i} = fullfile(a,['dc' b c]);
    dcrafunc_bold_files{i} = dcr_opt_func_bold_files{i};
    matlabbatch = []; 
    def = spm_get_defaults('coreg');
    %matlabbatch{1}.spm.spatial.coreg.estimate.ref = dcrafunc_sbref_files(i);
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[dcrafunc_bold_files{1} ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {t1_files};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions = def.estimate;
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run', {matlabbatch});
    
end
%spm_check_registration(t1_files, dcrafunc_bold_files{i});
%%
temp_t = fmri_data(dcrafunc_bold_files{1},dcrafunc_bold_files{1});

%% Normalization
matlabbatch = [];

load(which('segment_job.mat'));

%for i = 1:length(dcrafunc_bold_files)
for j = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(j).tpm{1} = [which('TPM.nii') ',' num2str(j)];
end
matlabbatch{1}.spm.spatial.preproc.channel.vols{1} = t1_files;

[b,c] = fileparts(t1_files);
deformation_nii = fullfile(b, ['y_' c '.nii']);
matlabbatch{2}.spm.spatial.normalise.write.subj.def = {deformation_nii};
% use dc
matlabbatch{2}.spm.spatial.normalise.write.subj.resample = dcrafunc_bold_files';

matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78  -112   -70
    78    76    85];
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run', {matlabbatch});

%
% warping anatomical image
clear matlabbatch;

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deformation_nii};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {t1_files};

matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78  -112   -70
    78    76    85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run', {matlabbatch});

%%
%wdcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','wdcra*CMRR*mb4_4D.nii'));
wdcrfunc_bold_files  = []; 
for i = 1:length(dcrafunc_bold_files)
    [a,b,c] = fileparts(dcrafunc_bold_files{i});
    wdcrfunc_bold_files{i} =fullfile(a,[sprintf('w%s',b) c]);
end
spm_check_registration(which('gray_matter_mask.nii'), [wdcrfunc_bold_files{1} ',1']);

%% Smoothing
fwhm = 5; % default fwhm
%wdcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','wdcra*CMRR*mb4_4D.nii'));
matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
matlabbatch{1}.spm.spatial.smooth.dtype = 0; % data type; 0 = same as before
matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit mask; 0 = no
matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(fwhm, 1, 3); % override whatever the defaults were with this
matlabbatch{1}.spm.spatial.smooth.data = wdcrfunc_bold_files';
run_num = [];

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
%%
spm_check_registration(which('gray_matter_mask.nii'), [swdcrfunc_bold_files{1} ',1']);
%%
swdcrafunc_bold_files = [];
for i = 1:length(dcrafunc_bold_files)
    [a,b,c] = fileparts(wdcrfunc_bold_files{i});
    swdcrfunc_bold_files{i} = fullfile(a,[sprintf('s%s',b) c]);
end

%swdcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','swdcra*CMRR*mb4_4D.nii'));
temp_dat = fmri_data(swdcrfunc_bold_files{1},which('gray_matter_mask.nii'));
temp_dat2 = fmri_data(swdcrfunc_bold_files{2},which('gray_matter_mask.nii'));

%%
plot(temp_dat2)