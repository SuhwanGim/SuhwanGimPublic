function cluster_surf_suhwan_four(r,depth, c)
rr = r;
surface_name_L = 'fsavg_left';
surface_name_R = 'fsavg_right';
surface_style ='veryinflated';
% surface_name_L = which('surf_workbench_inflated_32k_Left.mat');
% surface_name_R = which('surf_workbench_inflated_32k_Right.mat');
% surface_style = 'inflated'
%% left
%depth = 4;
axes_positions = {[0.02 0.5 .46 .5], [0.52 0.5 .46 .5], [0.02 0.1 .46 .5], [0.52 0.1 .46 .5]};
axes('Position', axes_positions{1});
%cluster_surf_jj(region(schaefer_mask,'unique_mask_values'), 2, which('surf_workbench_inflated_32k_Left.mat'),cols);
cluster_surf(rr, depth,surface_name_L,'colormaps', c, [], 'heatmap');
out.h = get(gca, 'children');
set(out.h(2), 'BackFaceLighting', 'lit')
camlight(-90,-20);
axis vis3d;
view(-90, 0);

% Right
axes('Position', axes_positions{2});
%cluster_surf_jj(region(schaefer_mask,'unique_mask_values'), 2, which('surf_workbench_inflated_32k_Right.mat'),cols);
cluster_surf(rr, depth, surface_name_R,'colormaps', c, [], 'heatmap');

if strcmp(surface_style,'veryinflated')
    %surface_light(gca);
    out.h = get(gca, 'children');
    set(out.h(2), 'BackFaceLighting', 'lit')
    camlight(-90,-20);
    axis vis3d;
else
    camlight(-90,-20); axis vis3d;
end
view(90, 0);

% left medial
axes('Position', axes_positions{3});
%cluster_surf_jj(region(schaefer_mask,'unique_mask_values'), 2, surface_name_L,cols);
cluster_surf(rr, depth, surface_name_L,'colormaps', c, [], 'heatmap');
camlight(-90,-20); axis vis3d;
view(90, 0);

% Right medial
axes('Position', axes_positions{4});
%cluster_surf_jj(region(schaefer_mask,'unique_mask_values'), 2, which('surf_workbench_inflated_32k_Right.mat'),cols);
cluster_surf(rr, depth, surface_name_R,'colormaps', c, [], 'heatmap');
%out.h = get(gca, 'children');
%set(out.h(2), 'BackFaceLighting', 'lit')
%camlight(-90,-20); axis vis3d;
if strcmp(surface_style,'veryinflated')
    camlight(-90,-20); axis vis3d;
else
    out.h = get(gca, 'children');
    set(out.h(2), 'BackFaceLighting', 'lit')
    camlight(-90,-20);
    axis vis3d;
end
view(-90, 0);
end
%%
% % Right medial
% axes('Position', axes_positions{4});
% %cluster_surf_jj(region(schaefer_mask,'unique_mask_values'), 2, which('surf_workbench_inflated_32k_Right.mat'),cols);
% cluster_surf(rr, depth, surface_name_R,'colormaps', c, [], 'heatmap');
% %out.h = get(gca, 'children');
% %set(out.h(2), 'BackFaceLighting', 'lit')
% % camlight(-90,-20); axis vis3d;
% %view(-90, 0);