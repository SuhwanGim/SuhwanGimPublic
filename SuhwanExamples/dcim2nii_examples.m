clear; close;
%% SET PATH
username = char(java.lang.System.getProperty('user.name'));
project_name = 'sync_SEP';
basedir = ['/Users/' username sprintf('/Dropbox/Projects/%s',project_name)];
addpath(genpath(basedir)); % add path
datdir = fullfile(basedir, 'data','coil_test','KSH_COIL_TEST_20210513','COCOAN_SUHWAN_GIM_20210513_113857_829000');
dcmSource = filenames(fullfile(datdir,'*'));
niiFolder = fullfile(basedir,'data','coil_test','KSH_COIL_TEST_20210513','results_nii');
%% DICOM 2 NIIFT 
for i = 1:length(dcmSource)
    if ~contains(dcmSource{i},'_unused')
        
        if contains(dcmSource{i},'T1_')
            [~,b] = fileparts(dcmSource{i});
            outputfold = fullfile(niiFolder,b);
            dicm2nii(dcmSource{i}, outputfold , 4);
        else
            [~,b] = fileparts(dcmSource{i});
            outputfold = fullfile(niiFolder,b);
            dicm2nii(dcmSource{i}, outputfold , 4);
            out = load(fullfile(outputfold, 'dcmHeaders.mat'));
            f = fields(out.h);
            %
            
            %cd(outputfold);
            
            nifti_3d = filenames(fullfile(outputfold,[f{1} '*.nii']));
            if contains(dcmSource{i},'SBREF_')
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
func_bold_files = filenames(fullfile(niiFolder,'CMRR_2*','CMRR*mb4_4D.nii'));
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
%%
adat = []; 
afunc_bold_files = filenames(fullfile(niiFolder,'CMRR_2*','aCMRR*mb4_4D.nii'));
for i = 1:length(afunc_bold_files)
    adat{i} = fmri_data(afunc_bold_files{i},mask{i}.fullpath);
end
plot(adat{1});