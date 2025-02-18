%
% https://canlab.github.io/_pages/first_level_design_matrix_exploration/first_level_design_matrix_exploration.html
%% Building HRF Build model: Onsets (with durations) convolve with HRF (need SPM on path)
basedir = '/Users/WIWIFH/Dropbox/Projects/EXCOOL/data';
%% Load RUN DATA
sub_i = 2; % subject number
run_i = 1; % run number 
loadnames = fullfile(basedir,[sprintf('EXC_task_Sub-EX%03d_run%02d',sub_i,run_i) '.csv']);
T = readtable(loadnames);
%% 1. SIMPLE MODEL (ALL TRIALS ARE IN ONE EVENT REGRESSOR)
TR = 1;             % repeation time in pulse sequence 
numberOfTR = 350;   % it could be difference accross participant 
len = TR * numberOfTR;

% set
events = []; 
events{1} = [T.onsets_cue T.durations_cue];         % First event regressor 
events{2} = [T.onsets_ratings T.durations_ratings]; % Second event regressor

X = onsets2fmridesign(events, 1, len, spm_hrf(1));  % Making HRF convolved event regressors
%% Plot it
create_figure('X');
h = plot_matrix_cols(zscore(X(:, 1:2)), 'horiz');   % Ignore last column 

%vertical%
%
% Customize
set(h, 'LineWidth', 3);
drawnow, snapnow
%% Another plot
create_figure('X 2nd plot');
h = plot_matrix_cols(zscore(X(:, 1:2)), 'horiz');
% Customize
set(h(1), 'LineWidth', 3, 'Color', [0    0.4470    0.7410]);
set(h(2), 'LineWidth', 3, 'Color', [0.8500    0.3250    0.0980]);
hh = plot_vertical_line(1 - .25); set(hh, 'LineStyle', '--');
hh = plot_vertical_line(2 - .25); set(hh, 'LineStyle', '--');
axis tight
axis off
drawnow, snapnow
%% 2. Single-trial cue model
% This model treats each cue event as a distinct event regressor.
%   - 18 event regressors for each cue event
%   - 1 event regressor for all rating events

% 1) Load data
sub_i = 2; % subject number
run_i = 1; % run number 
loadnames = fullfile(basedir,[sprintf('EXC_task_Sub-EX%03d_run%02d',sub_i,run_i) '.csv']);
T = readtable(loadnames);
%%
% 2) Set parameters
TR = 1;
numberOfTR = 350;
len = TR * numberOfTR;
events = []; 
for i = 1:length(T.onsets_cue)
    events{i} = [T.onsets_cue(i) T.durations_cue(i)];    
end
events{end+1} = [T.onsets_ratings T.durations_ratings]; % ratings

% 3) convolve HRF based on information of events 
X = onsets2fmridesign(events, 1, len, spm_hrf(1));
%%
create_figure('X 2nd plot');
%h = plot_matrix_cols(zscore(X(:, 1:2)), 'horiz');
h = plot_matrix_cols(zscore(X(:, 1:size(X,2)-1)), 'horiz');

%% 2. Single-trial rating model
% This model treats each rating event as a distinct event regressor.
%   - 1 event regressors for all cue events
%   - 18 event regressor for each rating event

% 1) Load data
sub_i = 2; % subject number
run_i = 1; % run number 
loadnames = fullfile(basedir,[sprintf('EXC_task_Sub-EX%03d_run%02d',sub_i,run_i) '.csv']);
T = readtable(loadnames);
%% 
% 2) Set parameters
TR = 1;
numberOfTR = 350;
len = TR * numberOfTR;
events = []; 
events{i} = [T.onsets_cue T.durations_cue];    
for i = 2:length(T.onsets_ratings)
    events{i} = [T.onsets_ratings(i) T.durations_ratings(i)]; % ratings
end
% 3) convolve HRF based on information of events 
X = onsets2fmridesign(events, 1, len, spm_hrf(1));
%% 
create_figure('X 2nd plot');
%h = plot_matrix_cols(zscore(X(:, 1:2)), 'horiz');
h = plot_matrix_cols(zscore(X(:, 1:size(X,2)-1)), 'horiz');