function [matSamLat,matLatSam,obj] = paa_Bernoulli(matFeatSam,nLat,options)

% PAA_BERNOULLI computes archetypal patterns from Bernoulli observations
%   [matSamLat,matLatSam] = paa_Bernoulli(matFeatSam,nLat) returns archetypal
%   loading matrix matSamLat, and archetypal factor matrix matLatSam, from
%   {0,1} valued observations stored in matrix matFeatSam, given number
%   of archetypes nLat. Each column of matFeatSam is an observation. The
%   archetypes can be computed from loading matrix, and observation matrix
%   as matFeatSam x matSamLat.
%
%   options is an optional structure specifying paramters,
%       eps, the convergence criteria, default is 10^-6, and
%       verbose, switch for textual display, default is false
%       display, switch for graphical display, default is false
%       maxIter, maximum number of iterations, default is 10000
%       matFeatLat, archetypes, if not empty then matFeatSam is projected on these archetypes
%
%   obj is an optional vector storing the value of objective function
%
%   copyright (c) Sohan Seth, sohan.seth@hiit.fi

if any(isnan(matFeatSam))
    error('nan in data matrix')
end

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
matFeatLat = options.matFeatLat;
if ~isempty(matFeatLat)
    nLat = size(matFeatLat, 2);
end

clear options

obj = zeros(maxIter, 1); obj(:) = Inf;

[nFeat, nSam] = size(matFeatSam);

if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

if display
    figureMain = figure('papertype', 'a4', 'paperposition', [0 0 3.5 3.5]);
end

% Negative Bernoulli likelihood
computeCost = @(matFeatLat,matLatSam,matFeatSam)...
    (sum(sum(-matFeatSam.*log(matFeatLat*matLatSam))) ...
    + sum(sum(-(1 - matFeatSam).*log(1 - matFeatLat*matLatSam))) );

% Initialization
matSamLat = rand(nSam,nLat);
matLatSam = rand(nLat,nSam);

matLatSamNorm = bsxfun(@rdivide,matLatSam,sum(matLatSam));
matSamLatNorm = bsxfun(@rdivide,matSamLat,sum(matSamLat));

% Convert 0/1 to double, and add small value to avoid division by zero
epsilon = 10^-12;
matFeatSam = double(matFeatSam(:));
matFeatSam(matFeatSam <= epsilon) = epsilon;
matFeatSam(matFeatSam >= 1 - epsilon) = 1 - epsilon;
matFeatSam = reshape(matFeatSam, nFeat, nSam);
matFeatSam_ = 1 - matFeatSam;

if isempty(matFeatLat)
    obj(1) = computeCost(matFeatSam * matSamLatNorm, matLatSamNorm, matFeatSam);
else
    obj(1) = computeCost(matFeatLat, matLatSamNorm, matFeatSam);
end
for iter = 1:maxIter
    if verbose
        if ~mod(iter, 10)
            fprintf('.')
        end
        if ~mod(iter, 100)
            fprintf(' [%d]\n',iter);
        end
    end
    
    if isempty(matFeatLat)
        % Update matLaSam
        matLatSam = matLatSam .* ( (matFeatSam * matSamLatNorm)' ...
            * (matFeatSam ./ (matFeatSam * matSamLatNorm * matLatSamNorm)) ...
            + (matFeatSam_ * matSamLatNorm)' ...
            * (matFeatSam_ ./ (matFeatSam_ * matSamLatNorm * matLatSamNorm)) ) ...
            ./ nFeat; %(repmat(sum(matFeatSam), nLat, 1) + repmat(sum(matFeatSam_), nLat, 1));
        matLatSamNorm = bsxfun(@rdivide,matLatSam,sum(matLatSam));
        
        % Update matSamLat
        normalization = diag((matFeatSam * matSamLatNorm)' ...
            * (matFeatSam ./ (matFeatSam * matSamLatNorm * matLatSamNorm) ) * matLatSamNorm') ...
            + diag((matFeatSam_ * matSamLatNorm)' ...
            * (matFeatSam_ ./ (matFeatSam_ * matSamLatNorm * matLatSamNorm) ) * matLatSamNorm');
        matSamLat = matSamLat .* ( matFeatSam' ...
            * (matFeatSam ./ (matFeatSam * matSamLatNorm * matLatSamNorm) ) * matLatSamNorm' ...
            + matFeatSam_' * (matFeatSam_ ./ (matFeatSam_ * matSamLatNorm * matLatSamNorm) ) ...
            * matLatSamNorm') ./ repmat(normalization(:)', nSam, 1);
        matSamLatNorm = bsxfun(@rdivide,matSamLat,sum(matSamLat));
        
        obj(iter+1) = computeCost(matFeatSam * matSamLatNorm, matLatSamNorm, matFeatSam);
    else
        % Update matLaSam
        matLatSam = matLatSam .* ( (matFeatLat)' ...
            * (matFeatSam ./ (matFeatLat * matLatSamNorm)) ...
            + (1 - matFeatLat)' ...
            * (matFeatSam_ ./ (1 - matFeatLat * matLatSamNorm)) ) ...
            ./ nFeat; %(repmat(sum(matFeatSam), nLat, 1) + repmat(sum(matFeatSam_), nLat, 1));
        matLatSamNorm = bsxfun(@rdivide,matLatSam,sum(matLatSam));
        
        obj(iter+1) = computeCost(matFeatLat, matLatSamNorm, matFeatSam);
    end
    
    % Convergence
    if abs(obj(iter+1) - obj(iter))/obj(iter) < eps
        if verbose
            fprintf('\nconvergence reached in %d iterations\n', iter)
        end
        break
    end
    
    % Display
    if display && ~mod(iter,10)
        matFeatLatPlot = matFeatSam * matSamLatNorm;
        figure(figureMain), clf, hold on,
        plot(matFeatSam(1,:),matFeatSam(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
        plot(matFeatLatPlot(1,:),matFeatLatPlot(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
        hold off, pause(0.01)
    end
end

matSamLat = matSamLatNorm;
matLatSam = matLatSamNorm;

if iter == maxIter && verbose
    fprintf('\nmaximum iteration reached\n')
end
fprintf('\n')
if nargout == 3
    obj(isinf(obj)) = [];
end

if ~isempty(matFeatLat)
    matSamLat = [];
end