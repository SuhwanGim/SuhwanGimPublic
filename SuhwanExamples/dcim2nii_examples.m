clear; close;
%% SET PATH
username = char(java.lang.System.getProperty('user.name'));
project_name = 'sync_SEP';
basedir = ['/Users/' username sprintf('/Dropbox/Projects/%s',project_name)];
addpath(genpath(basedir)); % add path
% JJ
% datdir = fullfile(basedir, 'data','coil_test','JJ_JJ_face_test_with_distortion_correction','COCOAN_SUHWAN_GIM_20210708_195321_871000');
% niiFolder = fullfile(basedir,'data','coil_test','JJ_JJ_face_test_with_distortion_correction','results_nii');
% Jihoon
datdir = fullfile(basedir, 'data','coil_test','JIHOON_JIHOON','COCOAN_SUHWAN_GIM_20210726_130239_993000');
niiFolder = fullfile(basedir,'data','coil_test','JIHOON_JIHOON','results_nii');

% ME (July 26)
datdir = fullfile(basedir, 'data','coil_test','20210726_ME','COCOAN_SUHWAN_GIM_20210726_130239_993000');
niiFolder = fullfile(basedir,'data','coil_test','20210726_ME','results_nii');

dcmSource = filenames(fullfile(datdir,'*'));
%% DICOM 2 NIIFT 
for i = 1:length(dcmSource)
    if ~contains(dcmSource{i},'_unused')
        
        if contains(dcmSource{i},'T1_')
            [~,b] = fileparts(dcmSource{i});
            outputfold = fullfile(niiFolder,b);
            dicm2nii(dcmSource{i}, outputfold , 4);
        elseif contains(dcmSource{i},'AAHEAD')
            % do nothing
        elseif contains(dcmSource{i}, 'LOCALIZER')
        else
            [~,b] = fileparts(dcmSource{i});
            outputfold = fullfile(niiFolder,b);
            dicm2nii(dcmSource{i}, outputfold , 4);
            out = load(fullfile(outputfold, 'dcmHeaders.mat'));
            f = fields(out.h);
            %
            
            %cd(outputfold);
            
            nifti_3d = filenames(fullfile(outputfold,[f{1} '*.nii']));
            if contains(dcmSource{i},'SBREF_') | contains(dcmSource{i},'DISTORTION')
                disdaq = 0;
            else
                disdaq = 18;
            end
            
            disp('Converting 3d images to 4d images...');
            output_4d_fnames = fullfile(outputfold, sprintf('%s_4D.nii',f{1}));
            spm_file_merge(nifti_3d((disdaq+1):end), output_4d_fnames);
            
            %system(['cd ' outputfold '; rm ' f{1} '*nii']);
            %system(['cd ' outputfold '; rm ' nifti_3d]);
            delete(nifti_3d{:})
            
        end
        %dicm2nii(dcmSource{i}, niiFolder, 4);
        
    end
end
%% SIMPLE PREPROC
%
mask = []; dat = []; 
spike_covariates =[]; 
func_bold_files = []; 
func_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','*CMRR*mb4_4D.nii'));
for i = 1:length(func_bold_files)
    clf;
    [a,b]=fileparts(func_bold_files{i});    
    % implicit_mask
    [~, ~, ~, ~, outputname] = fmri_mask_thresh_canlab(char(func_bold_files{i}),...
        fullfile(a, 'implicit_mask.nii'));
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
tr = 1; 
mbf = 4;

for i = 1:length(func_bold_files)
    slice_timing_job = []; 
    [a,b] = fileparts(func_bold_files{i});
    dicomheader = load(fullfile(a,'dcmHeaders.mat'));
    f = fields(dicomheader.h);
    eval(['slice_time = dicomheader.h.' f{1} '.MosaicRefAcqTimes;']);    
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
    if ~exist('custom_slice_timing', 'var')
        slice_timing_job{1}.spm.temporal.st.refslice = find(slice_time==0, 1, 'first');
    else
        slice_timing_job{1}.spm.temporal.st.refslice = find(slice_time==min(slice_time), 1, 'first');
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
%% Motion correction 
adat = []; 
afunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','a*CMRR*mb4_4D.nii'));
sbref_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','*CMRR*mb4_SBRef_4D.nii'));
for i = 1:length(afunc_bold_files)
    def = spm_get_defaults('realign');
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = def.estimate;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions = def.write;
    
    % change a couple things
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % do not register to mean (twice as long)
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0; % do not mask (will set data to zero at edges!)
    
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
%% Distortion correction 
%% add fsl path 
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

distort_pa_dat = fullfile(niiFolder,'DISTORTION_CORR_64CH_PA_0002','distortion_corr_64ch_pa_4D.nii');
distort_ap_dat = fullfile(niiFolder,'DISTORTION_CORR_64CH_PA_POLARITY_INVERT_TO_AP_0003','distortion_corr_64ch_pa_polarity_invert_to_ap_4D.nii');
dicomheader_files = filenames(fullfile(niiFolder,'*CMRR_2*MB4_00*','dcmHeaders.mat'));
rafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','ra*CMRR*mb4_4D.nii'));
epi_enc_dir = 'pa';
distortion_correction_out = fullfile(niiFolder, 'Distortion_combied','dc_combined.nii');
%system(['fslmerge -t ', distortion_correction_out, ' ', distort_ap_dat, ' ', distort_pa_dat]);
system(['fslmerge -t ', distortion_correction_out, ' ', distort_pa_dat, ' ', distort_ap_dat]);


%%
for i = 1:length(rafunc_bold_files) 
    % calculate and write the distortion correction parameter

    dicomheader = load(dicomheader_files{i});
    %readout_time = dicomheader.h.CMRR_2_7iso_ap_64ch_Flex_mb4.ReadoutSeconds;
    eval(sprintf('readout_time = dicomheader.h.%s.ReadoutSeconds;',char(fieldnames(dicomheader.h))));
    distort_info = nifti(distort_ap_dat);
    distort_num = distort_info.dat.dim(4);

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
    
    
    %% Applying topup on SBREF files
    input_dat = sbref_bold_files{i};
    [a, b] = fileparts(input_dat);
    dc_func_sbref_files{i,1} = fullfile(a, ['dc_' b '.nii']);
    system(['applytopup --imain=', input_dat, ' --inindex=1 --topup=', topup_out, ' --datain=', dc_param, ...
        ' --method=jac --interp=spline --out=', dc_func_sbref_files{i}]);
    
    % removing spline interpolation neg values by absolute
    system(['fslmaths ', dc_func_sbref_files{i}, ' -abs ', dc_func_sbref_files{i}]);
    
    % unzip
    system(['gzip -d -f ' dc_func_sbref_files{i} '.gz']);
    % system(['gzip -d ' PREPROC.dc_func_sbref_files{i} '.gz']);
    
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
    
    
end

%% Coregi
use_dc = true;
use_sbref = true;
dcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','dcra*CMRR*mb4_4D.nii'));
dcrafunc_sbref_files = filenames(fullfile(niiFolder,'*CMRR_2*','dc_*CMRR*mb4_SBRef_4D.nii'));
t1_files = fullfile(niiFolder, 'T1_MPRAGE_SAG_0_7ISO_0012','T1_mprage_sag_0_7iso.nii');
%t1_files = '/Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/JIHOON2_JIHOON2/results_nii/T1_MPRAGE_SAG_0_7ISO_0008/T1_mprage_sag_0_7iso.nii'; % 64-ch full coil
for i = 1:length(dcrafunc_bold_files)
    matlabbatch = []; 
    def = spm_get_defaults('coreg');
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = dcrafunc_sbref_files(i);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {t1_files};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions = def.estimate;
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run', {matlabbatch});
    
end
spm_check_registration(t1_files, dcrafunc_sbref_files{i});
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
matlabbatch{2}.spm.spatial.normalise.write.subj.resample = dcrafunc_bold_files;

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
wdcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','wdcra*CMRR*mb4_4D.nii'));
spm_check_registration(which('gray_matter_mask.nii'), [wdcrafunc_bold_files{1} ',1']);

%% Smoothing
fwhm = 5; % default fwhm
wdcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','wdcra*CMRR*mb4_4D.nii'));
matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
matlabbatch{1}.spm.spatial.smooth.dtype = 0; % data type; 0 = same as before
matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit mask; 0 = no
matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(fwhm, 1, 3); % override whatever the defaults were with this
matlabbatch{1}.spm.spatial.smooth.data = wdcrafunc_bold_files;
run_num = [];

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
%%
swdcrafunc_bold_files = filenames(fullfile(niiFolder,'*CMRR_2*','swdcra*CMRR*mb4_4D.nii'));
temp_dat = fmri_data(swdcrafunc_bold_files{1},which('gray_matter_mask.nii'));

%%
