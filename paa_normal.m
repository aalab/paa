function [matSamLat,matLatSam,matFeatLat,obj] = paa_normal(matFeatSam,nLat,options)

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
clear options

penalty = 200; % Large number to ensure stochasticity
problem.solver = 'lsqnonneg';
problem.options = optimset('maxIter',1000);

[nFeat, nSam] = size(matFeatSam);
if display
    figureMain = figure; 
end
if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

computeCost = @(matSamLat,matLatSam,matFeatSam,matFeatLat,reg)...
    (norm(matFeatSam - matFeatLat*matLatSam,'fro')^2 +...
        norm(matFeatLat - matFeatSam*matSamLat,'fro')^2);
    
% Initialization
matSamLat = rand(nSam,nLat); matSamLat = bsxfun(@rdivide,matSamLat,sum(matSamLat));
matLatSam = rand(nLat,nSam); matLatSam = bsxfun(@rdivide,matLatSam,sum(matLatSam));
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
    
    regS = ( (norm(matFeatLat - matFeatSam*matSamLat,'fro')^2) / (1 + nLat*nFeat/2)) ...
        / ( (norm(matFeatSam - matFeatLat*matLatSam,'fro')^2) / (1+ nSam*nFeat/2));
    temp = ((matLatSam * matLatSam' + regS*eye(nLat)) \ (matLatSam + regS*matSamLat'))';
    matFeatLat = matFeatSam * temp;
    
    obj(countIter+1) = computeCost(matSamLat,matLatSam,matFeatSam,matFeatLat);
    if verbose, fprintf('[cost %0.6f, change %0.6f%%]\n',obj(countIter+1),(obj(countIter+1) - obj(countIter))/obj(countIter)), end
    
    if abs(obj(countIter+1) - obj(countIter))/obj(countIter) < eps
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
obj(isinf(obj)) = [];
fprintf('\n')