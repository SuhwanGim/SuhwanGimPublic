function [mkr,h] = sepplot2(x, y, prop, varargin)
% Draw points and a shorter line plots between points. This function almost
% similar with sepplot in CanlabCore. 
%
% :Usage:
% ::
%
%    [mkr,h] = sepplot(x, y, prop, varargin)
%
% :Inputs:
%
%   **x, y:**
%        The function plots vector Y against vector X
%
%   **prop:**
%        The proportion of the lines: prop can be between 0 and 1
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'linecolor':**
%        followed by color (e.g., 'color', [.5 .5 .5]) (default = black)
%
%   **'linewidth':**
%        followed by a number for linewidth (e.g., 'linewidth', 2) (default = .5)
%
%   **'linestyle':**
%        linestyle, e.g., followed by '-', '--', ':' (default = '-')
%
% :Output:
%
%   **mkr:**
%       graphic handles for markers
%  
%   **h:**
%        graphic handles for lines
%
% :Examples: you can see the output in 
% http://wagerlab.colorado.edu/wiki/doku.php/help/core/figure_gallery
% ::
%
%    x = 1:5; % x values
%    y = [32 40 55 84 130]; % mean
%    e = [6 6 6 6 6]; % standard error of the mean
%
%    create_figure(y_axis);
%    set(gcf, 'Position', [1   512   268   194]);
%    col = [0.3333    0.6588    1.0000];
%    markercol = col-.2;
%
%    h = errorbar(x, y, e, 'o', 'color', 'k', 'linewidth', 1.5, 'markersize', 7, 'markerfacecolor', col);
%    hold on;
%    sepplot(x, y, .75, 'color', col, 'linewidth', 2);
%    errorbar_width(h, x, [0 0]);
%
%    set(gca, 'xlim', [.5 5.5], 'linewidth', 1.5);
%
%    try
%        pagesetup(gcf);
%        saveas(gcf, 'example.pdf');
%    catch
%        pagesetup(gcf);
%        saveas(gcf, 'example.pdf');
%    end
%
% ..
%    Copyright (C) 2014  Wani Woo
% ..
dcol = [1 1 1]; % Marker color default (white)
lcol = [0 0 0]+0.1; % Line color default (black)
linew = .5; % linewidth default
lines = '-'; % linestyle default
edgeCol = 'none';
mSize = 36;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'markercolor'
                dcol = varargin{i+1};
            case 'linecolor'
                lcol = varargin{i+1};
            case 'linewidth'
                linew = varargin{i+1};
            case 'linestyle'
                lines = varargin{i+1};
            case 'edgecolor'
                edgeCol = varargin{i+1};
            case 'markersize'
                mSize = varargin{i+1};
            otherwise
                error('Unknown inputs');
        end
    end
end
hold on;
%mkr=plot(x,y,'Marker','o','MarkerFaceColor',dcol,'MarkerEdgeColor', [edgeCol 0.8],'LineStyle','none');%,'MarkerSize',mSize);
mkr = scatter(x,y,mSize,dcol ,'o','filled','MarkerEdgeColor',[edgeCol 0.8]);
for i = 1:(numel(x)-1)
    xstep = (x(i+1)-x(i)).*((1-prop)/2);
    newx = [x(i)+xstep x(i+1)-xstep];
    
    ystep = (y(i+1)-y(i)).*((1-prop)/2);
    newy = [y(i)+ystep y(i+1)-ystep];
    
    h(i) = plot(newx, newy, 'color', lcol, 'linewidth', linew, 'linestyle', lines);
end

end
