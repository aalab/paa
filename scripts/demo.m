%% Demo for PAA codes
% copyright (c) Sohan Seth, sohan.seth@hiit.fi


%% Poisson observation
close all
rng default,
options = generate_options(); 
options.verbose = true; 
options.display = true;
options.eps = 10^-10;
options.maxIter = 100000;
options.matFeatLat = [];
n = 100;

% Generating data
matFeatSam = [ceil(rand(1,n)*10); ceil(rand(1,n)*10)];

% Learning archetypes
[matSamLat, matLatSam_1, ~] = paa_Poisson(matFeatSam, 4, options);
archetypes = matFeatSam * matSamLat;
axis([0 11 0 11]), box on, set(gca, 'fontsize', 20)

% Computing projections given archetypes
options.matFeatLat = archetypes;
options.display = false;
[~, matLatSam_2, ~] = paa_Poisson(matFeatSam, [], options);

fprintf('difference between projections %0.6f\n', ...
    norm(archetypes * matLatSam_2 - archetypes * matLatSam_1, 1) / numel(matFeatSam))

%% Bernoulli observation
close all
rng default,
options = generate_options(); 
options.verbose = true; 
options.display = true;
options.eps = 10^-6;
options.maxIter = 1000;
options.matFeatLat = [];
n = 500;

% Generating data
matFeatSam = [rand(1,n) > 0.8; rand(1,n) > 0.2];

% Learning archetypes
[matSamLat, matLatSam_1, ~] = paa_Bernoulli(matFeatSam, 4, options);
archetypes = matFeatSam * matSamLat;
axis([-0.1 1.1 -0.1 1.1]), box on, set(gca, 'fontsize', 20)

% Computing projections given archetypes
options.matFeatLat = archetypes;
options.display = false;
[~, matLatSam_2, ~] = paa_Bernoulli(matFeatSam, [], options);

fprintf('difference between projections %0.6f\n', ...
    norm(archetypes * matLatSam_2 - archetypes * matLatSam_1, 1) / numel(matFeatSam))

%% Multinomial observations
close all
rng default,
options = generate_options(); 
options.verbose = true; 
options.display = true;
options.eps = 10^-6;
options.maxIter = 10000;
options.matFeatLat = [];
n = 500;

% Generating data
matFeatSam = rand(3, n); matFeatSam = bsxfun(@rdivide, matFeatSam, sum(matFeatSam));
nFeatSam = mnrnd(1000, matFeatSam')';
matFeatSam = bsxfun(@rdivide, nFeatSam, sum(nFeatSam)); % empirical probabilities

% Learning archetypes
[matSamLat, matLatSam_1, ~] = paa_stochastic(nFeatSam, 3, options);
archetypes = matFeatSam * matSamLat;
axis([0 1000 0 1000]), box on, set(gca, 'fontsize', 20)

% Computing projections given archetypes
options.matFeatLat = archetypes;
options.display = false;
[~, matLatSam_2, ~] = paa_stochastic(nFeatSam, [], options);

fprintf('difference between projections %0.6f\n', ...
    norm(archetypes * matLatSam_2 - archetypes * matLatSam_1, 1) / numel(matFeatSam))

%% Normal observation
close all
rng default,
options = generate_options(); 
options.verbose = true; 
options.display = true;
options.eps = 10^-6;
options.maxIter = 20;
options.matFeatLat = [];
n = 100;

% Generating data
matFeatSam = [rand(1,n); rand(1,n)];

% Learning archetypes
[matSamLat, matLatSam_1, ~] = paa_normal(matFeatSam, 4, options);
archetypes = matFeatSam * matSamLat;
axis([0 1 0 1]), box on, set(gca, 'fontsize', 20)

% Computing projections given archetypes
options.matFeatLat = archetypes;
options.display = false;
[~, matLatSam_2, ~] = paa_normal(matFeatSam, [], options);

fprintf('difference between projections %0.6f\n', ...
    norm(archetypes * matLatSam_2 - archetypes * matLatSam_1, 1) / numel(matFeatSam))

%% Normal observation with R interface
close all
rng default,
% Generating data
matFeatSam = [rand(1,n); rand(1,n)];

% Learning archetypes
[matSamLat, matLatSam_1, ~] = classic_aa(matFeatSam, 4);
archetypes = matFeatSam * matSamLat;

% Computing projections given archetypes
matLatSam_2 = classic_aa_test(matFeatSam);

fprintf('difference between projections %0.6f\n', ...
    norm(archetypes * matLatSam_2 - archetypes * matLatSam_1, 1) / numel(matFeatSam))

publish('')