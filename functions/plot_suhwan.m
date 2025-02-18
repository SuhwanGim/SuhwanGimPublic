function h = plot_suhwan(x,y,varargin)
%%
% create plot with white backgroud and bigger size
%
% :: Usage
%       plot_suhwan(x,y,'color',[0,0,1],'linwidth',2);
%
% :: Optional input
%       1. color
%       2. linewidth
%
% :: Output (it works from R2014b)
%       1. plot: plot axis handle       
% 
% :: Examples
%       h = plot_suhwan(x,y,'color',[0,0,1],'linwidth',2);
%       plot_suhwan(x,y);
%
% Written by Suhwan Gim (19. 03. 22) 
%%
% Parse varargin
color = [0.3333, 0.6588, 1.0000];
linew = 2;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'color'}
                color = varargin{i+1};
            case {'linewidth'}
                linew = varargin{i+1};                
        end
    end
end
%%
facecolor = [0.2 0.2 0.35;
    0.5 0.5 1;
    1 0.5 0.5;
    1 0.7333 1];
figure;
set(gcf, 'Position', [-1919        -279        1920         984]);
set(gcf, 'color', 'w');

h = plot(x,y,'-','color',color','linewidth',linew);
xlabel('Time point', 'fontsize', 22);
ylabel('Effect magnitude', 'fontsize', 22);

set(gca, 'FontSize',22);
set(gca, 'linewidth', 1.5);%,'xlim',[0 N+0.5],'XTick',[1,N],'XTickLabel', [0,14.5],'ylim',[-0.15 1.5]);% ,'XTickLabelRotation',360-45);
end
