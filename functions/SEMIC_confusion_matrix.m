function [lss, percetage_matrix] = SEMIC_confusion_matrix(Y, Yfit, varargin)
%: Calculate frequency and draw confustion matrix.
%
% :: Usage: [lss, percetage_matrix] = SEMIC_confusion_matrix(Y, Yfit, varargin)
%
% :: Input 
%       1) Y and Yfit is from resulsts of fmri_dat.predict (especially,
%       algorithm 'CV_SVM'). 
% 
% :: Output
%       1)               lss : frequency matrix
%       2) percentage_matrix : percentage matrix (i.e., lss./number of data)
%
% :: Option
%       1) 'plot': draw confision matrix using imagesc (colormap: 'hot')
% 
% 
%
% Original scripts form confusion_matrix.m (in SpiderToolbox))
%
% Written by Suhwan Gim (28, Jan, 2019)
%% Set parameters
x_org = Y;
y_org = Yfit;
do_plot = false;
%% varargin
for i = 1:length(varargin)
    if ischar(varargin{i})
        % functional commands
        switch varargin{i}            
            case {'plot'}
                do_plot = true;                
        end
    end
end

%x_org = pred_model{1}.stats.Y;
%y_org = pred_model{1}.stats.yfit;
%%
nClasses = size( y_org, 2); %
% =============================== %
%            multi-class .        %
% =============================== %
if nClasses > 1
    x = ( x_org+1)./2;       % convert to class-numbers
    y = ( y_org+1)./2;
    % sanity-check
    if ( sum( sum( x, 2) ~= ones( size( x, 1), 1)) ~= 0 | ...
            sum( sum( y, 2) ~= ones( size( y, 1), 1)) ~= 0)
        error( '?? multi-class labels corrupted!')
    end
    labels = [ 1:nClasses]';
    x = x * labels;       % convert to class-numbers
    y = y * labels;
    
    for ii = 1:nClasses
        for jj = 1:nClasses
            lss( ii, jj) = sum( x==jj & y==ii);
        end
    end
    % =============================== %
    %            Binary-class.        %
    % =============================== %
else
    lss(1,1)= sum(x_org==1  & y_org==1);
    lss(2,1)= sum(x_org==1  & y_org==-1);
    lss(1,2)= sum(x_org==-1  & y_org==1);
    lss(2,2)= sum(x_org==-1  & y_org==-1);
end
percetage_matrix=lss./histcounts(x);
%%
if do_plot
    imagesc(percetage_matrix); colorbar; colormap('hot'); title('x = Y, y = yfit');
end
%dat=data([get_name(dat) ' -> confusion_matrix [tp, fp; fn, tn]' ],[],lss);