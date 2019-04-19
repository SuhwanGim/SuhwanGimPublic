function resample_FSL(in_file,out_file,varargin)
%

vox_size = 4;
do_verbose = true;
do_ref = false;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'vox_size'}
                vox_size = varargin{i+1};                            
            case {'no_verbose'}
                do_verbose = false;
            case {'ref'}
                temp_ref = varargin{i+1};
                do_ref = true;
        end
    end
end
%% ====================================================================== %
%                      Set environment for FSL                            %
% ======================================================================= %
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

setenv('FSLDIR', '/usr/local/fsl');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath; 
%% check input
% variable type
if ~iscell(in_file) && ~iscell(out_file)
    in_file{1} = in_file;
    out_file{1} = out_file;
elseif (~iscell(in_file) && iscell(out_file)) || (iscell(in_file) && ~iscell(out_file))
    error('in_file and out_files are not same, check the variables types')
else
    in_file = in_file;
    out_file = out_file;
end
% comparing input and output name
id_idx1 = []; id_idx2 = [];
for file_i = 1:numesl(in_file)
    for file_ii = 1:numesl(out_file)
        if strcmp(in_file{file_i}, out_file{file_ii})
            id_idx1 = [id_idx1 file_i];
            id_idx2 = [id_idx2 file_i];
        end   
    end
    
end
% if 
if sum(id_idx1) > 0 || sum(id_idx2) >0
    disp('See list');
    disp('=============================================================== ')
    disp(['index of in_file: ' id_idx1]);
    disp(['index of out_file: ' id_idx2]);
    disp('=============================================================== ')
    error('There is identical file name');
end

%% resample using FSL
for file_i = 1:numel(in_file)    
    % SETUP: imgs and reference 
    temp_imgs = in_file{file_i};       
    if do_ref 
        temp_ref = inpu_ref;
    else %if not, using own image as reference/
        temp_ref = in_file{file_i};
    end
    out_name = out_file{file_i};
    %% VERBOSE
    if do_verbose
        tic;
        disp('---------------------------------------------------------------');
        disp(['Original file name: ' temp_imgs ]);
        disp(['Processed file name: ' out_file ]);
    end
    %% RUN 
    %       FSL: flirt 
    system(['export FSLOUTPUTTYPE=NIFTI;'  ...
        ...
        'flirt' ...
        ' -in ' temp_imgs ...
        ' -ref ' temp_ref ...
        ' -applyisoxfm ' num2str(vox_size) ' '  ...
        ' -interp trilinear'  ...
        ' -noresampblur'  ...
        ' -out ' out_name]);
    
    if do_verbose
        toc;
    end
end
disp('_____________________________________________________        END   ');
end
