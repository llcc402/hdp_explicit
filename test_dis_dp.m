G0 = gem(100, 5);

x = discreternd(G0, 1000);

G1 = gem(100, 5); % prior
alph = .1;

count = histcounts(x, 1:101);
a = alpha * G1 + count;
b = [cumsum(a(2:end), 'reverse'), 0];

V = betarnd(repmat(a, 100, 1),repmat(b, 100, 1));
G2 = V;
V = cumprod(1-V, 2);
G2(:,2:end) = V(:,1:end-1) .* G2(:,2:end);
G2 = mean(G2);

plot(1:100, G0, 'o', 1:100, G2, '*')