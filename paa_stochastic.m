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

[nFeat, nSam] = size(nFeatSam);

if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

if display
    figureMain = figure;
end

% Normalize each column
matFeatSam = nFeatSam ./ repmat(sum(nFeatSam), nFeat, 1);

% Initialization
matSamLat = rand(nSam,nLat); matSamLat = bsxfun(@rdivide,matSamLat,sum(matSamLat));
matLatSam = rand(nLat,nSam); matLatSam = bsxfun(@rdivide,matLatSam,sum(matLatSam));

computeCost = @(nFeatSam, matFeatSam, matSamLat, matLatSam, epsilon)...
    (sum(sum(nFeatSam .* log(epsilon + matFeatSam * matSamLat * matLatSam),2)));

epsilon = 10^-12;  % Avoiding log(0) and division by 0
obj(1) = computeCost(nFeatSam, matFeatSam, matSamLat, matLatSam, epsilon);
for iter = 1:maxIter
    if verbose
        if ~mod(iter, 10)
            fprintf('.')
        end
        if ~mod(iter, 100)
            fprintf(' [%d]\n',iter);
        end
    end
    
    % Expectation
    temp = nFeatSam ./ (matFeatSam * matSamLat * matLatSam);
    
    matSamLatNew = (epsilon + matFeatSam' * temp * matLatSam') .* matSamLat;
    matLatSamNew = (epsilon + matSamLat' * matFeatSam' * temp) .* matLatSam;
    
    % Maximization
    matSamLat = matSamLatNew ./ repmat(sum(matSamLatNew), nSam, 1);
    matLatSam = matLatSamNew ./ repmat(sum(matLatSamNew), nLat, 1);
    
    % Convergence
    obj(iter+1) = computeCost(nFeatSam, matFeatSam, matSamLat, matLatSam, epsilon);    
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
        hold off,
    end 
end

if iter == maxIter && verbose
    fprintf('\nmaximum iteration reached\n')    
end
fprintf('\n')
if nargout == 3
    obj(isinf(obj)) = [];
end