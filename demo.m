% Demo for PAA codes
% copyright (c) Sohan Seth, sohan.seth@hiit.fi

options = generate_options(); 
options.verbose = false; 
options.display = true;
options.eps = 10^-8;
options.maxIter = 50000;
n = 500;

% Poisson observation
matFeatSam = [ceil(rand(1,n)*10); ceil(rand(1,n)*10)];
[~, ~, ~] = paa_Poisson(matFeatSam,4,options);

% Normal observation
matFeatSam = [rand(1,n); rand(1,n)];
[~, ~, ~] = paa_normal(matFeatSam,4,options);

% Bernoulli observation
matFeatSam = [rand(1,n) > 0.3; rand(1,n) > 0.6];
[~, ~, ~] = paa_Bernoulli(matFeatSam,4,options);

% Term-frequency observation
matFeatSam = rand(3,n); matFeatSam = matFeatSam ./ repmat(sum(matFeatSam), 3, 1);
nFeatSam = 100 * ceil( repmat(ones(1, n) * 1000, 3, 1) .* matFeatSam);
[~, ~, ~] = paa_stochastic(nFeatSam, 3, options);