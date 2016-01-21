% Input:
%      M        a scalar. The number of documents.
%      D        a scalar. The size of the dictionary.
%      lambda   a scalar. The expected number of words in each documents.
%      actN     a scalar. The maximum number of activated atoms.
%      gamma    a scalar. The concentration parameter of G0.
%      alpha    a scalar. The concentration parameter of G1 to GM.
%      beta     a scalar. The parameter of the topics.
% Output:
%      mixing   a matrix of order M * actN. The mixing measure of each document.
%      data     a matrix of order M * actN. The documents.
%      topic    a matrix of order actN * D.
function [data, mixing, topic] = data_generate(M, D, lambda, actN, gamma, alpha, beta)
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
    actN = 100;
end
if nargin < 5
    gamma = 5;
end
if nargin < 6
    alpha = 1;
end
if nargin < 7
    beta = 1;
end

%--------------------------------------------------------------------------
% STEP 1: Generate topic
%--------------------------------------------------------------------------
topic = dirichletrnd(beta * ones(actN, D));

%--------------------------------------------------------------------------
% STEP 2: Generate mixing measure
%--------------------------------------------------------------------------
mixing = zeros(M, actN);
G0 = gem(actN, gamma);
for i = 1:M
    mixing(i,:) = dpDisrnd(alpha, G0);
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
    