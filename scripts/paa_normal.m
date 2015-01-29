function [matSamLat, matLatSam, obj] = paa_normal(matFeatSam, nLat, options)

% PAA_NORMAL computes archetypal patterns from Normal observations
%   This code uses the nonnegative least squares solver... SLOW IMPLEMENTATION
%   [matSamLat,matLatSam] = paa_normal(matFeatSam,nLat) returns archetypal
%   loading matrix matSamLat, and archetypal factor matrix matLatSam, from
%   real valued observations stored in matrix matFeatSam, given number
%   of archetypes nLat. Each column of matFeatSam is an observation. The
%   archetypes can be computed from loading matrix, and observation matrix
%   as matFeatSam x matSamLat.
%
%   options is an optional structure specifying paramters,
%       eps, the convergence criteria, default is 10^-6, and
%       verbose, switch for textual display, default is false
%       display, switch for graphical display, default is false
%       maxIter, maximum number of iterations, default is 10000
%
%   obj is an optional vector storing the value of objective function
%
%   The code uses a PENALTY parameter to enforce stochasticity, default 200.
%
%   copyright (c) Sohan Seth, sohan.seth@hiit.fi
global alpha_0 beta_0
alpha_0 = 1; beta_0 = 1;

if nargin < 2
    error('Observation matrix and number of archetypes must be provided');
end

if nargin == 2
    options = generate_options();
end

eps = options.eps;
verbose = options.verbose;
display = options.display;
maxIter = options.maxIter;
matFeatLat = options.matFeatLat; LEARN_ARCHETYPES = true;
if ~isempty(matFeatLat)
    nLat = size(matFeatLat, 2); LEARN_ARCHETYPES = false;
end
clear options

debug = false;
if verbose, fprintf('alpha %0.2f beta %0.2f\n', alpha_0, beta_0); end

penalty = 200; % Large number to ensure stochasticity
problem.solver = 'lsqnonneg';
problem.options = optimset('maxIter',1000,'Display', 'off');

[nFeat, nSam] = size(matFeatSam);
if display
    figureMain = figure('papertype', 'a4', 'paperposition', [0 0 3.5 3.5]);
end
if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

% computeCost = @(matSamLat,matLatSam,matFeatSam,matFeatLat,reg)...
%     (norm(matFeatSam - matFeatLat*matLatSam,'fro')^2 +...
%     norm(matFeatLat - matFeatSam*matSamLat,'fro')^2);

% Initialization
matSamLat = rand(nSam,nLat); matSamLat = bsxfun(@rdivide,matSamLat,sum(matSamLat));
matLatSam = rand(nLat,nSam); matLatSam = bsxfun(@rdivide,matLatSam,sum(matLatSam));

if LEARN_ARCHETYPES
    matFeatLat = matFeatSam * matSamLat;
    obj = Inf * ones(maxIter+1,1);
    obj(1) = computeCost(matSamLat,matLatSam,matFeatSam,matFeatLat);
    
    for countIter = 1:maxIter
        if verbose, fprintf('[Iteration %d]\n',countIter), elseif ~mod(countIter,10), fprintf('.'), end
        
        problem.C = [matFeatSam; penalty*ones(1,nSam)];
        for countCol = 1:nLat
            problem.d = [matFeatLat(:,countCol); penalty];
            problem.x0 = matSamLat(:,countCol);
            matSamLat(:,countCol) = lsqnonneg(problem);
        end
        
        problem.C = [matFeatLat; penalty*ones(1,nLat)];
        for countCol = 1:nSam
            problem.d = [matFeatSam(:,countCol); penalty];
            problem.x0 = matLatSam(:,countCol);
            matLatSam(:,countCol) = lsqnonneg(problem);
        end
        
        T = computeCost(matSamLat,matLatSam,matFeatSam,matFeatLat);
        alpha_1 = alpha_0 + nSam*nFeat/2;
        beta_1 = beta_0 + norm(matFeatSam - matFeatLat*matLatSam,'fro')^2/2;
        alpha_2 = alpha_0 + nLat*nFeat/2;
        beta_2 = beta_0 + norm(matFeatLat - matFeatSam*matSamLat,'fro')^2/2;
        regS = (alpha_2 / beta_2) / (alpha_1 / beta_1);
        %regS = 1/regS;
        if debug; fprintf('eps1 %0.4f eps2 %0.4f regularizer %0.4f\n', ...
                alpha_1/beta_1, alpha_2/beta_2, regS); end
        temp = ((matLatSam * matLatSam' + regS*eye(nLat)) \ (matLatSam + regS*matSamLat'))';
        matFeatLat = matFeatSam * temp;
        
        obj(countIter+1) = computeCost(matSamLat,matLatSam,matFeatSam,matFeatLat);
        if obj(countIter+1) < T
            error('cost decreased')
        end
        if debug, fprintf('[cost %0.6f, change %0.6f%%]\n', ...
                obj(countIter+1),abs((obj(countIter+1) - obj(countIter))/obj(countIter))), end
        
        if abs((obj(countIter+1) - obj(countIter))/obj(countIter)) < eps
            break
        end
        if display && ~mod(countIter,10)
            figure(figureMain), clf, hold on,
            plot(matFeatSam(1,:),matFeatSam(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
            matFeatLat_ = matFeatSam * matSamLat;
            plot(matFeatLat_(1,:),matFeatLat_(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
            hold off; drawnow
        end
    end
else
    matLatSam = rand(nLat, nSam); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
    problem.C = [matFeatLat; penalty*ones(1,nLat)];
    for countCol = 1:nSam
        problem.d = [matFeatSam(:,countCol); penalty];
        problem.x0 = matLatSam(:,countCol);
        matLatSam(:,countCol) = lsqnonneg(problem);
    end
    obj = [];
end

if ~isempty(matFeatLat)
    obj(isinf(obj)) = [];
end
fprintf('\n')

% Variational bound
function val = entGamma(alpha, beta)
val = -log(beta) + alpha + gammaln(alpha) + (1 - alpha) * psi(alpha);
function val = logGamma(alpha, beta)
val = psi(alpha) - log(beta);
function val = computeCost(matSamLat,matLatSam,matFeatSam,matFeatLat)
global alpha_0 beta_0
alpha_1 = alpha_0 + numel(matFeatSam)/2;
beta_1 = beta_0 + norm(matFeatSam - matFeatLat*matLatSam,'fro')^2/2;
alpha_2 = alpha_0 + numel(matFeatLat)/2;
beta_2 = beta_0 + norm(matFeatLat - matFeatSam*matSamLat,'fro')^2/2;
val = - alpha_1 - alpha_2 + ...
    numel(matFeatSam) * logGamma(alpha_1, beta_1) / 2 + ...
    numel(matFeatLat) * logGamma(alpha_2, beta_2) / 2 + ...
    entGamma(alpha_1, beta_1) + entGamma(alpha_2, beta_2);