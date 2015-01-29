function [matSamLat, matLatSam, obj] = paa_stochastic(nFeatSam, nLat, options)

% PAA_STOCHASTIC computes archetypal patterns from term frequency observations
%   [matSamLat,matLatSam] = paa_stochastic(nFeatSam,nLat) returns archetypal
%   loading matrix matSamLat, and archetypal factor matrix matLatSam, from
%   term frequency observations stored in matrix matFeatSam, given number
%   of archetypes nLat. Each column of matFeatSam is an observation. The
%   archetypes can be computed from loading matrix, and observation matrix
%   as matFeatSam x matSamLat where matFeatSam is column normalized nFeatSam.
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

if any(isnan(nFeatSam))
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

[nFeat, nSam] = size(nFeatSam);

if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

if display
    figureMain = figure('papertype', 'a4', 'paperposition', [0 0 3.5 3.5]);
end

% Normalize each column
matFeatSam = nFeatSam ./ repmat(sum(nFeatSam), nFeat, 1);

% Initialization
matSamLat = rand(nSam,nLat); matSamLat = bsxfun(@rdivide,matSamLat,sum(matSamLat));
matLatSam = rand(nLat,nSam); matLatSam = bsxfun(@rdivide,matLatSam,sum(matLatSam));

% Negative log likelihood
computeCost = @(matFeatLat, matLatSam, nFeatSam, epsilon)...
    (- sum(sum(nFeatSam .* log(epsilon + matFeatLat * matLatSam),2)));

epsilon = 10^-16;  % Avoiding log(0) and division by 0
if isempty(matFeatLat)
    obj(1) = computeCost(matFeatSam * matSamLat, matLatSam, nFeatSam, epsilon);
else
    obj(1) = computeCost(matFeatLat, matLatSam, nFeatSam, epsilon);
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
        % Expectation
        temp = nFeatSam ./ (matFeatSam * matSamLat * matLatSam);
        
        matSamLatNew = (epsilon + matFeatSam' * temp * matLatSam') .* matSamLat;
        matLatSamNew = (epsilon + matSamLat' * matFeatSam' * temp) .* matLatSam;
        
        % Maximization
        matSamLat = bsxfun(@rdivide,matSamLatNew,sum(matSamLatNew));
        matLatSam = bsxfun(@rdivide,matLatSamNew,sum(matLatSamNew));
        
        % Convergence
        obj(iter+1) = computeCost(matFeatSam * matSamLat, matLatSam, nFeatSam, epsilon);;
    else
        % Expectation
        temp = nFeatSam ./ (matFeatLat * matLatSam);
        
        matLatSamNew = (epsilon + matFeatLat' * temp) .* matLatSam;
        
        % Maximization
        matLatSam = bsxfun(@rdivide,matLatSamNew,sum(matLatSamNew));
        
        obj(iter+1) = computeCost(matFeatLat, matLatSam, nFeatSam, epsilon);
    end
    
    % Convergence
    if abs((obj(iter+1) - obj(iter))/obj(iter)) < eps
        if verbose
            fprintf('\nconvergence reached in %d iterations\n', iter)
        end
        break;
    end
    
    % Display
    if display && ~mod(iter,10)
        figure(figureMain), clf, hold on,
        plot(nFeatSam(1,:), nFeatSam(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
        nFeatLat = nFeatSam * matSamLat;
        plot(nFeatLat(1,:), nFeatLat(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
        hold off, pause(0.01)
    end
end

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