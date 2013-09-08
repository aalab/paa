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
%
%   obj is an optional vector storing the value of objective function
%
%   The code uses a PENALTY parameter to enforce stochasticity, default 20.
%       The variables are normalized to be bounded between [0,1].
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

% Penalty and scaling
matFeatSamScale = max(matFeatSam, [], 2); %ones(nFeat, 1); % $ No normalization $ %max(matFeatSam, [], 2);
matFeatSam = matFeatSam ./ repmat(matFeatSamScale, 1, nSam);
penalty = 20 * max(matFeatSam(:));
epsilon = 10^-12; matFeatSam = matFeatSam + epsilon;

% Penalized negative Poisson likelihood
computeCost = @(matSamLat, matLatSam, matFeatSam, penalty)...
    (sum(sum(- matFeatSam .* log(matFeatSam * matSamLat * matLatSam) ...
    + (matFeatSam * matSamLat * matLatSam ))));
% computeCost = @(matSamLat, matLatSam, matFeatSam, penalty)...
%     (sum(sum(- matFeatSam .* log(matFeatSam * matSamLat * matLatSam) ...
%     + (matFeatSam * matSamLat * matLatSam ))) ...
%     + penalty * sum(-log(sum(matLatSam)) + sum(matLatSam)) ...
%     + penalty * sum(-log(sum(matSamLat)) + sum(matSamLat)));

% Initialization
matSamLat = rand(nSam,nLat); matSamLat = bsxfun(@rdivide,matSamLat,sum(matSamLat));
matLatSam = rand(nLat,nSam); matLatSam = bsxfun(@rdivide,matLatSam,sum(matLatSam));

obj(1) = computeCost(matSamLat, matLatSam, matFeatSam, penalty);
for iter = 1:maxIter
    if verbose
        if ~mod(iter, 10)
            fprintf('.')
        end
        if ~mod(iter, 100)
            fprintf(' [%d]\n',iter);
        end
    end
    
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
    
    % Convergence
    obj(iter+1) = computeCost(matSamLat,matLatSam,matFeatSam,penalty);
    if abs(obj(iter+1) - obj(iter)/obj(iter)) < eps
        if verbose
            fprintf('\nconvergence reached in %d iterations\n', iter); keyboard
        end
        break
    end
    
    % Display
    if display && ~mod(iter,10)
        figure(figureMain), clf, hold on,
        plot(matFeatSamScale(1) * matFeatSam(1,:), matFeatSamScale(2) * matFeatSam(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
        matFeatLat = diag(matFeatSamScale) * matFeatSam * matSamLat;
        plot(matFeatLat(1,:),matFeatLat(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
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
