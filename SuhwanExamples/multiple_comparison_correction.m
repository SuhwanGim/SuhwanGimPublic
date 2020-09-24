%% BASIC SETTINGS

rng('default');
p_matrix = rand(10000,1)./100;


q= 0.05;
%% External Toolbox
% 1. Multiple Testing Toolbox
%    REF1: https://www.mathworks.com/matlabcentral/fileexchange/70604-multiple-testing-toolbox
%    REF2: http://www.gib.tel.uva.es/members.php?lang=en_EN#martinezcagigal_v)
%[c_pvalues, c_alpha, h] = fwer_bonf(p_matrix, 0.05);
%[c_pvalues, c_alpha, h] = fdr_BY(p_matrix,q,'unknown',true);
[qvalues, fdr, h] = fdr_storey(p_matrix,q,true);
%% In-house toolbox
% 1. CanlabCore
%    REF1: github.com/canlab/canlabcore
corrected_p_thresh = FDR(pr_array(:),q);
if isempty(corrected_p_thresh), corrected_p_thresh = -Inf;end
fprintf('Corrected threshold p-value: %03d \n\n',corrected_p_thresh);

% 2. CocoanCore
%    REF1: github.com/cocoanlab/cocoanCore




p_thresh=getFDR(pr_array,q);
fprintf('toTl: %03d \n\n',length(find(pr_array < p_thresh)));

%% Matlab implemented Function
mafdr