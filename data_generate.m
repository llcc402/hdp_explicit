% Input:
%      M        a scalar. The number of documents.
%      D        a scalar. The size of the dictionary.
%      lambda   a scalar. The expected number of words in each document.
%      K        a scalar. The maximum number of topics in each document.
%      alpha    a scalar. The concentration parameter of G1 to GM.
%      beta     a scalar. The parameter of the topics.
% Output:
%      mixing   a matrix of order M * actN. The mixing measure of each document.
%      data     a matrix of order M * actN. The documents.
%      topic    a matrix of order actN * D.
function [data, mixing, topic] = data_generate(M, D, lambda, K, alpha, beta)
if nargin < 1
    M = 100;
end
if nargin < 2
    D = 500;
end
if nargin < 3
    lambda = 80;
end
if nargin < 4
    K = 10;
end
if nargin < 5
    alpha = 1;
end
if nargin < 6
    beta = 1;
end

%--------------------------------------------------------------------------
% STEP 1: Generate topic
%--------------------------------------------------------------------------
topic = dirichletrnd(ones(K, D));
% make the topics concentrate on some important words
topic = dirichletrnd(beta * topic);

%--------------------------------------------------------------------------
% STEP 2: Generate mixing measure
%--------------------------------------------------------------------------
mixing = zeros(M, K);
for i = 1:M
    ix = find(rand(1, K) > .3);
    mixing(i, ix) = dirichletrnd(alpha * ones(1, length(ix)));
end

%--------------------------------------------------------------------------
% STEP 3: Generate documents
%--------------------------------------------------------------------------

% generate lengths of the documents
len = poissrnd(lambda, 1, M);

% generate data
data = zeros(M, D);
for i = 1:M
    for j = 1:len(i)
        % sample a topic
        t = discreternd(mixing(i,:), 1);
        % sample a word
        w = discreternd(topic(t,:), 1);
        data(i,w) = data(i,w) + 1;
    end
end
end
    