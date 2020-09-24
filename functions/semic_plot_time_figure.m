
function semic_plot_time_figure(brain_idx, rating_idx)
%%
% :: semic_plot_time_figure(brain_idx, rating_idx)
% 
%    brain_idx : 1~45
%    rating_idx: 1~32
%

%%
ii=2; % middle color
ribbon_color2 = [213,62,79
    244,109,67
    253,174,97]./255;
lnwidth = 6;
create_figure('time_plot');
%set(gcf, 'position', [1        1155         259         190]);
set(gcf,'position',[-1301         639         162         121]); % for examples
axis([-0.3, 1.3, -0.3, 1])
set(gca, 'YDir', 'reverse');
axis off;

river = riverplot_draw_ribbon([0 (min(brain_idx))/45], [0.8 min(rating_idx/45)+0.14], [0 (max(brain_idx))/45], [0.8 max(rating_idx/45)+0.14], 'steepness', .05, 'color', ribbon_color2(4-ii,:));
river.patchh.EdgeAlpha = 1;
river.patchh.EdgeColor = ribbon_color2(4-ii,:);%ribbon_color;
river.patchh.FaceAlpha= 1;
river.patchh.LineWidth = 0.4;
river.line1.h.LineStyle = 'none';
river.line2.h.LineStyle = 'none';


% DRAW BAR
line([0;0], [1/45;14/45], 'linewidth', lnwidth, 'color', [.7 .7 .7]);
line([0;0], [14/45;31/45], 'linewidth', lnwidth, 'color', [.5 .5 .5]);
line([0;0], [31/45;45/45], 'linewidth', lnwidth, 'color', [.3 .3 .3]);


line([0.8;0.8], [1/45;4/45]+0.14, 'linewidth', lnwidth, 'color', 'k');
line([0.8;0.8], [4/45;14/45]+0.14, 'linewidth', lnwidth, 'color', [78 182 230]/255);
line([0.8;0.8], [14/45;23/45]+0.14, 'linewidth', lnwidth, 'color', [117 107 177]/255);
line([0.8;0.8], [23/45;32/45]+0.14, 'linewidth', lnwidth, 'color', [232 126 114]/255);


end