function  p = perm_test(x,y,num_rounds)
%
% References: permutation_test in mlxtend

rng('shuffle');
observeddifference = nanmean(x(:)) - nanmean(y(:));
xy = [x(:) y(:)];
randomdifferences = zeros(1,num_rounds);
for i=1:num_rounds
    all_perm=randperm(length(xy));
    
    randomSample1 = xy(all_perm(1:length(x)));
    randomSample2 = xy(all_perm(length(x)+1:length(all_perm)));
    
    % saving differences between the two samples
    randomdifferences(i) = nanmean(randomSample1) - nanmean(randomSample2);
end

p = (length(find(abs(randomdifferences) > abs(observeddifference)))+1) / (num_rounds+1);
end

