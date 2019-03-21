function figure_suhwan(varargin)

fn = 1;
for i = 1:length(varargin)   
    if isnumeric(varargin{i})
        fn = varargin{i};
    else
        warning('Invalid input for figure index (not numeric value)');
        fn = 1;
    end
end

figure(fn)
set(gcf, 'Position', [-1919        -279        1920         984]);
set(gcf, 'color', 'w');

end
