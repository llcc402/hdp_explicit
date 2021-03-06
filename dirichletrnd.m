% Input:
%     alpha      a positive matrix.
% Output:
%     x          a matrix of the same order as alpha. Each row is a
%                dirichlet distributed vector whose parameters are the
%                corresponding row of alpha.
function x = dirichletrnd(alpha)
if sum(sum(alpha <= 0)) > 0
    error('alpha must be strictly positive')
end

x = gamrnd(alpha, 1);
s = sum(x, 2);
x = x ./ repmat(s, 1, size(alpha, 2));
end