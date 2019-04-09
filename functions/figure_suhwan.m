function figure_Handle = figure_suhwan(varargin)
%%
% create figure with white backgroud and bigger size.
%
% :: Usage
%   figure_Handle = figure_suhwan(n);
%
% :: Optional input
%       1. n: figure number
% 
% :: Output (it works from R2014b)
%       1. figure_Handle: figure handle        
% 
% :: Examples
%       figure_suhwan(2);
%       figure_suhwan;
%
% Written by Suhwan Gim (19. 03. 22) 
%%
fn = 1; %figure number
sz = [-1919        -279        1920         984];
%title_str=[];

for i = 1:length(varargin)
%     if isnumeric(varargin{i})
%         fn = varargin{i};
%     end
    
    if ischar(varargin{i})
        switch lower(varargin{i})
            % functional commands
            case {'number'}
                fn = varargin{i+1};
            %    title_str = varargin{i+1};
            case {'size','position'}
                sz = varargin{i+1};
        end
    end
    

end
%%
figure_Handle = figure(fn);
set(gcf, 'Position', sz);
set(gcf, 'color', 'w');

end
