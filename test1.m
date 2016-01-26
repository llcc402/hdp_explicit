clear
clc

%--------------------------------------------------------------------------
% STEP 1: Generate data
%--------------------------------------------------------------------------
N = 500; % number of words in each doc
mixing = [.5 .5 0 0 0 0; 0 .5 .5 0 0 0; 0 0 .5 .5 0 0; 0 0 0 .5 .5 0; 0 0 0 0 .5 .5];

topic = zeros(6, 100);
for i = 1:6
    ix = randperm(100, 3);
    for j = 1:3
        topic(i,ix(j)) = rand(1);
    end
    topic(i,:) = topic(i,:) / sum(topic(i,:));
end

data = zeros(5, 100);
for i = 1:5
    for j = 1:N
        % sample latent variable
        z = discreternd(mixing(i,:), 1);
        % sample word
        w = discreternd(topic(z,:), 1);
        data(i,w) = data(i,w) + 1;
    end
end

ix = find(sum(data) > 0);
data = data(:,ix);

%--------------------------------------------------------------------------
% STEP 2: Init hdp pamameters
%--------------------------------------------------------------------------
% concentration parameter for G0
gamma = .5;
% concentration parameter for G1 to G5
alpha = .1;
% hyper parameter for topics
beta = .1;

% init all in topic 1
Z = ones(size(data)) .* (data > 0);

% assume at most activate 100 topics
actN = 100;

% maximum number of iterations for gibbs sampling
maxIter = 5000 * 2;

% mixing measures
mixing_post = zeros(5, actN);

% topics
topic_post = zeros(actN, size(data, 2));

%--------------------------------------------------------------------------
% STEP 3: HDP sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    % sample G0
    counts = accumarray(Z(:)+1, data(:));
    counts = counts(2:end);
    if length(counts) < actN
        a = [counts', zeros(1, actN - length(counts))];
    elseif length(counts) == actN
        a = counts';
    end
        
    b = [cumsum(a(2:end), 'reverse'), 0];

    a = a + 1;
    b = b + gamma;
    
    V = betarnd(a,b);
    G0 = V;
    V = cumprod(1-V);
    G0(2:end) = G0(2:end) .* V(1:end-1);
    
    % sample G1 to G5
    for i = 1:5
        counts = accumarray((Z(i,:)+1)',data(i,:));
        counts = counts(2:end);
        counts = counts';
        if length(counts) < actN
            counts = [counts, zeros(1,actN - length(counts))];
        end
        
        a = alpha * G0 + counts;
        b = [cumsum(a(2:end), 'reverse'), 0];
        
        V = betarnd(a,b);
        mixing_post(i,:) = V;
        V = cumprod(1-V);
        mixing_post(i,2:end) = mixing_post(i,2:end) .* V(1:end-1);
    end
    
    % sampling topics
    for i = 1:actN
        counts = sum((Z == i) .* data) + beta;
        topic_post(i,:) = dirichletrnd(counts);
    end
    
    % sample Z
    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            if data(i,j) > 0
                log_w = log(mixing_post(i,:)) + log(topic_post(:,j)') * data(i,j);
                log_w = log_w - max(log_w);
                w = exp(log_w);
                w = w / sum(w);
                w = cumsum(w);
                [~,~,Z(i,j)] = histcounts(rand(1), [0,w]);
            end
        end
    end
end
               

