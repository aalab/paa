function [matSamLat, matLatSam, obj] = paa_Poisson(matFeatSam, nLat, options)

% PAA_POISSON computes archetypal patterns from Poisson observations
%   [matSamLat,matLatSam] = paa_Poisson(matFeatSam,nLat) returns archetypal
%   loading matrix matSamLat, and archetypal factor matrix matLatSam, from
%   integer valued observations stored in matrix matFeatSam, given number
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
%   The code uses a PENALTY parameter to enforce stochasticity, default 20.
%       The variables are normalized to be bounded between [0,1].
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

debug = false;

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

% Penalty and scaling
matFeatSamScale = max(matFeatSam, [], 2); %ones(nFeat, 1); % $ No normalization $ %max(matFeatSam, [], 2);
matFeatSam = matFeatSam ./ repmat(matFeatSamScale, 1, nSam);
if ~isempty(matFeatLat)
    matFeatLat = matFeatLat ./ repmat(matFeatSamScale, 1, size(matFeatLat, 2));
end
penalty = 20 * max(matFeatSam(:));
epsilon = 10^-12; matFeatSam = matFeatSam + epsilon;

% Negative Poisson likelihood
computeCost = @(matFeatLat, matLatSam, matFeatSam, penalty)...
    (sum(sum(- matFeatSam .* log(matFeatLat * matLatSam) ...
    + (matFeatLat * matLatSam ))));
% Penalized negative Poisson likelihood
% computeCost = @(matSamLat, matLatSam, matFeatSam, penalty)...
%     (sum(sum(- matFeatSam .* log(matFeatSam * matSamLat * matLatSam) ...
%     + (matFeatSam * matSamLat * matLatSam ))) ...
%     + penalty * sum(-log(sum(matLatSam)) + sum(matLatSam)) ...
%     + penalty * sum(-log(sum(matSamLat)) + sum(matSamLat)));

% Initialization
matSamLat = rand(nSam,nLat); matSamLat = bsxfun(@rdivide,matSamLat,sum(matSamLat));
matLatSam = rand(nLat,nSam); matLatSam = bsxfun(@rdivide,matLatSam,sum(matLatSam));

if isempty(matFeatLat)
    obj(1) = computeCost(matFeatSam * matSamLat, matLatSam, matFeatSam, penalty);
else
    obj(1) = computeCost(matFeatLat, matLatSam, matFeatSam, penalty);
end
for iter = 1:maxIter
    if verbose
        if ~mod(iter, 100)
            fprintf('.')
        end
        if ~mod(iter, 1000)
            fprintf(' [%d]\n',iter);
        end
    end
    
    if isempty(matFeatLat)
        % Update matLaSam
        %B = [matFeatSam; penalty * ones(1,nSam)];
        %A = [matFeatSam * matSamLat; penalty * ones(1,nLat)];
        %gradPosLatSam = repmat(sum(A,1)', 1, nSam);
        %gradNegLatSam = A' * (B ./ (A * matLatSam));
        %matLatSam = matLatSam .* gradNegLatSam ./ gradPosLatSam;
        
        gradPosLatSam = repmat((sum(matFeatSam) * matSamLat)' + penalty, 1, nSam);
        gradNegLatSam = matSamLat' * matFeatSam' * (matFeatSam ...
            ./ (matFeatSam * matSamLat * matLatSam)) ...
            + penalty * repmat(1 ./ sum(matLatSam,1), nLat, 1);
        matLatSam = matLatSam .* gradNegLatSam ./ gradPosLatSam;
        
        % Update matSamLat
        %b = [matFeatSam(:); penalty * ones(nLat, 1)];
        %A = [kron(matLatSam', matFeatSam); zeros(nLat, nSam * nLat)];
        %for countLat = 1:nLat
        %A(nFeat * nSam + countLat, ((countLat-1) * nSam + 1) : countLat * nSam) = penalty;
        %end
        %gradPosSamLat = sum(A,1)';
        %gradNegSamLat = A' * (b ./ (A * matSamLat(:)));
        %matSamLat = reshape(matSamLat(:) .* gradNegSamLat ./ gradPosSamLat, nSam, nLat);
        
        gradPosSamLat = (sum(matFeatSam',2) * sum(matLatSam')) + penalty;
        gradNegSamLat = matFeatSam' * (matFeatSam ...
            ./ (matFeatSam * matSamLat * matLatSam)) ...
            * matLatSam' + penalty * repmat(1 ./ sum(matSamLat,1), nSam, 1);
        matSamLat = matSamLat .* (gradNegSamLat ./ gradPosSamLat);
        
        obj(iter+1) = computeCost(matFeatSam * matSamLat, matLatSam, matFeatSam,penalty);
        if debug == 1; fprintf('obj %0.2f\n', obj(iter + 1)); end
        
    else
        
        gradPosLatSam = repmat((sum(matFeatLat))' + penalty, 1, nSam);
        gradNegLatSam = matFeatLat' * (matFeatSam ...
            ./ (matFeatLat * matLatSam)) ...
            + penalty * repmat(1 ./ sum(matLatSam,1), nLat, 1);
        matLatSam = matLatSam .* gradNegLatSam ./ gradPosLatSam;
        obj(iter+1) = computeCost(matFeatLat, matLatSam, matFeatSam,penalty);
    end
    
    % Convergence
    if abs((obj(iter+1) - obj(iter))/obj(iter)) < eps
        if verbose
            fprintf('\nconvergence reached in %d iterations\n', iter);
        end
        break
    end
    
    % Display
    if display && ~mod(iter,100)
        figure(figureMain), clf, hold on,
        plot(matFeatSamScale(1) * matFeatSam(1,:), matFeatSamScale(2) * matFeatSam(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
        matFeatLatPlot = diag(matFeatSamScale) * matFeatSam * matSamLat;
        plot(matFeatLatPlot(1,:),matFeatLatPlot(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
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