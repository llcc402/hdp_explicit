% Input:
%     prob       a row vector. The probs of each bin.
%     n          a scalar. The number of observations.
% Output:
%     x          a row vector of length n.
function x = discreternd(prob, n)
if sum(prob) ~= 1
    prob = prob / sum(prob);
end

r = rand(1, n);
prob = [0, cumsum(prob)];
[~, ~, x] = histcounts(r, prob);
end
