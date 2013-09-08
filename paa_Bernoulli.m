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
%
%   obj is an optional vector storing the value of objective function
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

obj = zeros(maxIter, 1); obj(:) = Inf;

[nFeat, nSam] = size(matFeatSam);

if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

if display
    figureMain = figure;
end

% Negative Bernoulli likelihood
computeCost = @(matSamLat,matLatSam,matFeatSam)...
    (sum(sum(-matFeatSam.*log(matFeatSam*matSamLat*matLatSam))) ...
        + sum(sum(-(1 - matFeatSam).*log(1 - matFeatSam*matSamLat*matLatSam))) );
    
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

obj(1) = computeCost(matSamLatNorm,matLatSamNorm,matFeatSam);
for iter = 1:10000
    if verbose
        if ~mod(iter, 10)
            fprintf('.')
        end
        if ~mod(iter, 100)
            fprintf(' [%d]\n',iter);
        end
    end
    
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
       
    % Convergence
    obj(iter+1) = computeCost(matSamLatNorm,matLatSamNorm,matFeatSam);
    if abs(obj(iter+1) - obj(iter))/obj(iter) < eps
        if verbose
            fprintf('\nconvergence reached in %d iterations\n', iter)
        end
        break
    end
    
    % Display
    if display && ~mod(iter,10)
       matFeatLat = matFeatSam * matSamLatNorm;
       figure(figureMain), clf, hold on,
       plot(matFeatSam(1,:),matFeatSam(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
       plot(matFeatLat(1,:),matFeatLat(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
       hold off, pause(0.1)
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