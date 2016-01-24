% Input:
%     data        a matrix. The rows are corresponding to documents and the
%                 columns are corresponding to words in the dictionary.
%     gamma       a scalar. The concentration parameter of G0. 
%     alpha       a scalar. The concentration parameter of G1, ..., GM,
%                 where M is the number of documents.
%     beta        a scalar. The parameter of the topics.    
%     actN        a scalar. The maximum number of activated atoms of the
%                 mixing measure.
%     maxIter     a scalar. The number of Gibbs iterations.
% Output:
%     mixing_post a matrix of order M * actN. 
%     topic_post  a matrix of order actN * D, where D is the size of the
%                 dictionary.
function [mixing_post, topic_post] = hdp(data, gamma, alpha, beta, actN, maxIter)
if nargin < 5
    actN = 100;
end
if nargin < 6
    maxIter = 1000;
end

%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
topic_post = zeros(actN, size(data, 2));
mixing_post = zeros(size(data, 1), actN);

% Z is the latent variables, we assume all words are sampled from topic 1
% initially
Z = ones(size(data)) .* (data > 0);

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------

for iter = 1:maxIter
    
    % sample the positions of G0
    for k = 1:actN
        counts = sum(data .* (Z == k)) + beta;
        topic_post(k,:) = dirichletrnd(counts);
    end
%     topic_post(topic_post < 1e-200) = 1e-100;
    
    % sample weights of G0
    a = histcounts(Z(:), 1:actN+1);
    b = [cumsum(a(2:end), 'reverse'), 0];
    a = a + 1;
    b = b + gamma;
    V = betarnd(a, b);
    G0_weights = V;
    V = cumprod(1 - V);
    G0_weights(2:end) = G0_weights(2:end) .* V(1:end-1);
    
    % sample G1, ..., GM
    for i = 1:actN
        counts = histcounts(Z(i,:), 1:actN+1);
        a = alpha * G0_weights + counts;
        b =[cumsum(a(2:end), 'reverse'), 0];
        
        V = betarnd(a,b);
        mixing_post(i,:) = V;
        V = cumprod(1-V);
        mixing_post(i,2:end) = mixing_post(i,2:end) .* V(1:end-1); 
    end
    
    % sample Z
    for i = 1:size(Z, 1)
        for j = 1:size(Z, 2)
            if data(i,j) ~= 0
                log_p = log(mixing_post(i,:)) + log(topic_post(:,j)') * data(i,j);
                log_p = log_p - max(log_p); % in case all probs are 0
                prob = exp(log_p);
                prob = prob / sum(prob);
                if sum(isnan(prob)) > 0
                    fprintf(num2str(log_p))
                    h = 1
                end
                prob = cumsum(prob);
                [~, ~, Z(i,j)] = histcounts(rand(1), [0, prob]);
            end
        end
    end
    fprintf('iteration %d done\n', iter)
end

end
    