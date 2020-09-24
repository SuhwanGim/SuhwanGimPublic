function [r_idx, r_out]=get_region_idx(r, k)
%% 
% :: Usage
%       - [r_idx, r_out] = get_region_idx(r, k);
% :: Input 
%       - r: region_object
%       - k: number of voxels
% :: Output
%       - r_idx: region index of under k 
%       - r_out: eliminated region object using r_idx
%
% Written by Suhwan Gim (2019.09.23);

%%
r_del_idx = [];
for i=1:length(r)
    if (r(i).numVox <= k)
        r_del_idx = [r_del_idx i];
    end
end
r(r_del_idx) = [];

r_idx = r_del_idx; 
r_out = r;
%orthviews(rr);
return