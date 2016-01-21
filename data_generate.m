% Input:
%      M        a scalar. The number of documents.
%      D        a scalar. The size of the dictionary.
%      lambda   a scalar. The expected number of words in each documents.
%      actN     a scalar. The maximum number of activated atoms.
%      gamma    a scalar. The concentration parameter of G0.
%      alpha    a scalar. The concentration parameter of G1 to GM.
% Output:
%      mixing   a matrix of order M * actN. The mixing measure of each document.
%      data     a matrix of order M * actN. The documents.
[data, mixing] = data_generate(M, D, lambda, actN, gamma, alpha)

