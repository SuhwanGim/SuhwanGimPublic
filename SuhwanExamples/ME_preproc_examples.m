%% ADDpath
addpath(genpath('/Users/suhwan/Dropbox/github/CanlabTools'))
addpath(genpath('/Users/suhwan/Dropbox/github/CocoanTools/cocoanCORE'));
addpath(genpath('/Users/suhwan/Dropbox/github/CocoanTools/humanfmri_preproc_bids'));
addpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12'))
addpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12_updates_r7487'))
addpath(genpath('/Users/suhwan/Dropbox/Projects/SEMIC2/scripts/mediation_dream_toolbox'))
rmpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12/external/fieldtrip/compat'));
rmpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12_updates_r7487/external/fieldtrip/compat'))
rmpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12_updates_r7487'))
%% SIMPLE PREPROCESS FOR MUTLI-ECHO IMAGES 
clear; close;
%% SET PATH
username = char(java.lang.System.getProperty('user.name'));
project_name = 'sync_SEP';
project_name = 'EXCOOL';
%basedir = ['/Users/' username sprintf('/Dropbox/Projects/%s',project_name)];
basedir = '/Volumes/homeo/EXCOOL';
%addpath(genpath(basedir)); % add path
rawdir = fullfile(basedir,'raw');
%% subject directory 
fundir = [];
anadir = [];
fmapdir = []; 
for sub_i = 1:7
    fundir{sub_i} = filenames(fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'dicom','func','func*'));
    anadir{sub_i} = filenames(fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'dicom','anat','T1*'));
    fmapdir{sub_i} = filenames(fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'dicom','fmap','*ISO*'));
end
%% DICOM 2 NIIFT for anat
sub_i = 1;
for i = 1:length(anadir{sub_i})
    %addpath    
    [~,b] = fileparts(anadir{sub_i}{i});
    outputdir = fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'anat',b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,anadir{sub_i}{i})); % on sein
    %system(sprintf('/Users/admin/dcm2niix/build/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % on homeo    
end
%% DICOM 2 NIIFT for fmap
sub_i = 1;
for i = 1:length(fmapdir{sub_i})
    %addpath    
    [~,b] = fileparts(fmapdir{sub_i}{i});
    outputdir = fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'fmap',b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,fmapdir{sub_i}{i})); % on sein
    %system(sprintf('/Users/admin/dcm2niix/build/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % on homeo    
end
%% DICOM 2 NIIFT for func
% dcm2niix /Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/20210803_multi_echo/COCOAN_SUHWAN_GIM_20210803_111531_408000/1_2_3ME_HALF_CMRR_2_7ISO_PA_64CH+FLEX_MB4_IA2_0
sub_i = 1;
for i = 1:length(fundir{sub_i})
    %addpath    
    [~,b] = fileparts(fundir{sub_i}{i});
    outputdir = fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'func',b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,fundir{sub_i}{i})); % on sein
    %system(sprintf('/Users/admin/dcm2niix/build/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % on homeo    
end
%% MAKE DIRECTORY
sub_i =1;
prodir = fullfile(basedir,'preprocessed',sprintf('sub-EXCOOL%03d',sub_i)); 
preproc_func_dir = fullfile(prodir, 'func');
if ~exist(preproc_func_dir, 'dir'), mkdir(preproc_func_dir); end
preproc_mean_func_dir = fullfile(prodir, 'mean_func');
if ~exist(preproc_mean_func_dir, 'dir'), mkdir(preproc_mean_func_dir); end
preproc_anat_dir = fullfile(prodir, 'anat');
if ~exist(preproc_anat_dir, 'dir'), mkdir(preproc_anat_dir); end
preproc_fmap_dir = fullfile(prodir, 'fmap');
if ~exist(preproc_fmap_dir, 'dir'), mkdir(preproc_fmap_dir); end
qcdir = fullfile(prodir, 'qc_images');
if ~exist(qcdir, 'dir'), mkdir(qcdir); end
%% COPY FILES FROM RAW TO PREPROCESSED
sub_i = 1;
subjanatdir = []; 
for i = 1:length(anadir{sub_i})
    [~,b] = fileparts(anadir{sub_i}{i});    
    copyfile([fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'anat',b) '/*'], fullfile(preproc_anat_dir,b));
    subjanatdir{sub_i}{i} = fullfile(preproc_anat_dir,b);
end
subjfmapdir=[]; 
for i = 1:length(fmapdir{sub_i})
    [~,b] = fileparts(fmapdir{sub_i}{i});    
    copyfile([fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'fmap',b) '/*'], fullfile(preproc_fmap_dir,b));
    subjfmapdir{sub_i}{i} = fullfile(preproc_fmap_dir,b);
end
subjfmapdir{sub_i} = subjfmapdir{sub_i}';

subjdir = [];
for i = 1:length(fundir{sub_i})
    [~,b] = fileparts(fundir{sub_i}{i});    
    copyfile([fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'func',b) '/*'], fullfile(preproc_func_dir,b));
    subjdir{sub_i}{i} = fullfile(preproc_func_dir,b);
end
subjdir{sub_i} = subjdir{sub_i}';
%% SIMPLE PREPROCESSTING
mask = []; dat = [];
spike_covariates =[];
sub_i = 1;
subjdat=subjdir{sub_i};
subjdat(contains(subjdat,'sbref')) = [];

for i = 1:length(subjdat) % each task
    
    temp_list = filenames(fullfile(subjdat{i},'func_*.nii'));
    mask = []; 
    dat = []; 
    for e_i = 1:length(temp_list) % each echo
        clf;
        [a,b]=fileparts(temp_list{e_i});
        % implicit_mask
        [~, ~, ~, ~, outputname] = fmri_mask_thresh_canlab(char(temp_list{e_i}),...
            fullfile(a, sprintf('implicit_mask_e%02d.nii',e_i)));
        implicit_mask_file = outputname;
        mask{e_i} = fmri_data(implicit_mask_file,implicit_mask_file);
        dat{e_i} = fmri_data(temp_list{e_i}, implicit_mask_file);
        dat{e_i}.images_per_session = size(dat{e_i}.dat,2);
        % spike id
        diary(fullfile(a, ['qc_diary_' b '.txt']));
        dat{e_i} = preprocess(dat{e_i}, 'outliers', 'plot');  % Spike detect and globals by slice
        subplot(5, 1, 5);
        dat{e_i} = preprocess(dat{e_i}, 'outliers_rmssd', 'plot');  % RMSSD Spike detect
        diary off;
        
        sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
        set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
        drawnow;
        
        qcspikefilename = fullfile(a, ['qc_spike_plot_' b '.png']); % Scott added some lines to actually save the spike images
        saveas(gcf,qcspikefilename);
        spike_covariates{e_i} = dat{e_i}.covariates;        
    end
    save(fullfile(a,'spike_covariates_all.mat','spike_covariates'));
    
end

%% Slice time correction 
sub_i = 1;
subjdat=subjdir{sub_i};
subjdat(contains(subjdat,'sbref')) = [];

for i = 1:length(subjdat) % each task
    
    temp_list = filenames(fullfile(subjdat{i},'func_*.nii'));
    a = [];
    b = []; 
    for e_i = 1:length(temp_list) % each echo
        
        slice_timing_job = [];
        json_read = []; json_file = [];
        % Read Json file
        [a,b] = fileparts(temp_list{e_i});
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
        slice_timing_job{1}.spm.temporal.st.scans{1} = spm_select('expand', temp_list(e_i)); % individual 4d images in cell str
        %% 1. nslices
        Vfirst_vol = spm_vol([temp_list{e_i} ',1']);
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
end
%% Check the affine matrix (?)
temp_t = [];
for i =1:4
    %temp_t{i} = spm_vol([func_bold_files{i} ',2']); % raw images
    %temp_t{i} = spm_vol([afunc_bold_files{i} ',2']); % after slice-timing correction 
    %temp_t{i} = spm_vol([rafunc_bold_files{i} ',2']); % after realinment
end
%% REALIGNMENT: Motion correction 
% The consensus is to do 1) estimate realigment parameter using
% before-slice timing correction images and 2) alignment using after-slice
% timing correction images
%
% By J.J
%
sub_i = 1;
subjdat=subjdir{sub_i};
subjdat(contains(subjdat,'sbref')) = [];


for i = 1:length(subjdat) % each task
    afunc_bold_files = filenames(fullfile(subjdat{i},'afunc_*e*.nii'));
    a = [];
    b = [];
    for e_i = 1:length(afunc_bold_files) % each echo
        [a,b]=fileparts(afunc_bold_files{e_i});
        def = spm_get_defaults('realign');
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = def.estimate;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions = def.write;
        
        % change a couple things
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % do not register to mean (twice as long)
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0; % do not mask (will set data to zero at edges!)
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0]; % do not output mean image
        
        data = []; 
        data_all = afunc_bold_files(e_i);        
        matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = data_all;
        
        %% RUN
        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_jobman('run', {matlabbatch});
        
        
        %% Save realignment parameter
        
        [d, f] = fileparts(data_all{1});
        
        tempcpfile = fullfile(d, ['rp_' f '.txt']);
        temp_mvmt = textread(tempcpfile);
        
        %temp_mvmt(1:(images_per_session+kk-1),:) = [];
        
        % save plot
        create_figure('mvmt', 2, 1)
        subplot(2,1,1);
        plot(temp_mvmt(:,1:3));
        legend('x', 'y', 'z');
        
        subplot(2,1,2);
        plot(temp_mvmt(:,4:6));
        legend('pitch', 'roll', 'yaw');
        
        sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
        set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
        drawnow;
        
        %[~,a] = fileparts(afunc_bold_files{e_i});
        [a,b]=fileparts(afunc_bold_files{e_i});
        mvmt_qcfile = fullfile(a, ['qc_mvmt_' b '.png']); % Scott added some lines to actually save the spike images
        saveas(gcf,mvmt_qcfile);
        close all;
    end
end
%% ========================================================================
%
%                 RUN TEDANA for ME-EPI images
%            1) inputs: realigned three echos images 
%            2) outputs: several images, but combined one 
% =========================================================================
%
%
%     see: https://tedana.readthedocs.io/en/stable/outputs.html
%     ===============================   =================================================
%     Filenames                         Explanation
%     ===============================   =================================================
%     desc-optcom_bold.nii.gz           Optimally combined timeseries. 
%     desc-optcomDenoised_bold.nii.gz   (Recommended) Denoised & optimally combined time series. 
%     ===============================   =================================================

% ======================================================================
% you can access your denoised results 
% -> https://rica-fmri.netlify.app
% ======================================================================
% ABOUT the order of preprocessing
% -> https://tedana.readthedocs.io/en/stable/multi-echo.html#processing-multi-echo-fmri
%%
sub_i = 1;
subjdat=subjdir{sub_i};
subjdat(contains(subjdat,'sbref')) = [];
a = [];
b = [];
e_i = 2;
rafunc_bold_files = [];    % Single-echo image
combined_bold_files1 = []; % Optimally combined
combined_bold_files2 = []; % ME-ICA denoised
for i = 1:length(subjdat)    
    
    rafunc_bold_files{i} = filenames(fullfile(subjdat{i},'rafunc_*e2.nii'));    
    [a,b]=fileparts(rafunc_bold_files{i});
                
    combined_bold_files1{i} = fullfile(a,'tedana_output','desc-optcom_bold.nii');
    if ~isfile(combined_bold_files1{i})
        gunzip(fullfile(a,'tedana_output','desc-optcom_bold.nii.gz')); %from T2S_workflow in tedana
    end
    
    combined_bold_files2{i} = fullfile(a,'tedana_output','desc-optcomDenoised_bold.nii');
    if ~isfile(combined_bold_files2{i})
        gunzip(fullfile(a,'tedana_output','desc-optcomDenoised_bold.nii.gz')); %from T2S_workflow in tedana
    end
end

%% Three conditions
% 1) second echo images
% 2) optimally combined
% 3) optimally combeind + ICA (denoised)
%% AFTER TEDANA
% Distortion correction 
%% add fsl path 
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
sub_i = 1;

distort_pa_dat = filenames(fullfile(subjfmapdir{sub_i}{1},'*.nii'));
distort_ap_dat = filenames(fullfile(subjfmapdir{sub_i}{2},'*.nii'));

epi_enc_dir = 'pa';
distort_info = nifti(distort_ap_dat);
distort_num = distort_info.dat.dim(4);


[a,~] = fileparts(subjfmapdir{1}{1});
output_dc_combined = fullfile(a,'Distortion_combied');
if ~isfolder(output_dc_combined)
    mkdir(output_dc_combined)
end
distortion_correction_out = fullfile(output_dc_combined,'dc_combined.nii');
system(['fslmerge -t ', distortion_correction_out, ' ', distort_pa_dat{1}, ' ', distort_ap_dat{1}]);
%%
distort_json{1} = filenames(fullfile(subjfmapdir{sub_i}{1},'*.json'));
distort_json{2} = filenames(fullfile(subjfmapdir{sub_i}{2},'*.json'));
% Manually exported from Json files .TotalReaduoutTime; Three echo files' are same

% Read Json file
rdtime = NaN(1,2) ;
for i =1:2
    fid = fopen(distort_json{i}{1});
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    json_read = jsondecode(str);
    rdtime(1,i) = json_read.TotalReadoutTime;
end




readout_time = rdtime(1);
[a,~ ] = fileparts(rafunc_bold_files{1});
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
combined_bold_files{1} = rafunc_bold_files{1}{1};
combined_bold_files{2} = combined_bold_files1{1};
combined_bold_files{3} = combined_bold_files2{1};
dcr_func_bold_files = [];
for i = 1:length(combined_bold_files)
    
    input_dat = combined_bold_files{i};
    [a, b] = fileparts(input_dat);
    dcr_func_bold_files{i,1} = fullfile(a, ['dc' b '.nii']);
    system(['applytopup --imain=', input_dat, ' --inindex=1 --topup=', topup_out, ' --datain=', dc_param, ...
        ' --method=jac --interp=spline --out=', dcr_func_bold_files{i}]);
    
    % removing spline interpolation neg values by absolute
    system(['fslmaths ', dcr_func_bold_files{i}, ' -abs ', dcr_func_bold_files{i}]);
    
    % unzip
    system(['gzip -d -f ' dcr_func_bold_files{i} '.gz']);
    
    
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
    
    mean_dcr_func_bold_png = fullfile(a, 'mean_dcr_func_bold.png'); % Scott added some lines to actually save the spike images
    canlab_preproc_show_montage( mdat.fullpath , mean_dcr_func_bold_png);
    drawnow;
    
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
%% Coregister
use_dc = true;
dcrafunc_bold_files = []; 
sub_i = 1;
t1_files = filenames(fullfile(subjanatdir{sub_i}{1},'*.nii'),'char');

for i = 1:length(dcr_func_bold_files)
%     [a,b,c] = fileparts(rafunc_bold_files{i});
%     dcrafunc_bold_files{i} = fullfile(a,['dc' b c]);
    dcrafunc_bold_files{i} = dcr_func_bold_files{i};
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
%spm_check_registration(t1_files, [dcrafunc_bold_files{i} ', 1']);

%% Normalization
matlabbatch = [];

load(which('segment_job.mat'));

%for j = 1:length(dcrafunc_bold_files)
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
wdcrfunc_bold_files  = []; 
for i = 1:length(dcrafunc_bold_files)
    [a,b,c] = fileparts(dcrafunc_bold_files{i});
    wdcrfunc_bold_files{i} =fullfile(a,[sprintf('w%s',b) c]);
end
spm_check_registration(which('gray_matter_mask.nii'), [wdcrfunc_bold_files{3} ',1']);

%% Smoothing
fwhm = 5; % default fwhm
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
swdcrafunc_bold_files = [];
for i = 1:length(dcrafunc_bold_files)
    [a,b,c] = fileparts(wdcrfunc_bold_files{i});
    swdcrfunc_bold_files{i} = fullfile(a,[sprintf('s%s',b) c]);
end
%%
spm_check_registration(which('gray_matter_mask.nii'), [swdcrfunc_bold_files{3} ',1']);
%%
temp_dat = fmri_data(swdcrfunc_bold_files{1},which('gray_matter_mask.nii'));
temp_dat2 = fmri_data(swdcrfunc_bold_files{2},which('gray_matter_mask.nii'));
temp_dat3 = fmri_data(swdcrfunc_bold_files{3},which('gray_matter_mask.nii'));


%%
plot(temp_dat2)
