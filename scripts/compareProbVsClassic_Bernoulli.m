% Qualitative assessment of probabilistic archetypal analysis and standard
% archetypal analysis by Cutler and Breiman on different observation models
%
% Additional info: Figures 4, 5 and 6 in article Probabilistic Archetypal
% Analysis by Seth and Eugster in Machine Learning Journal
%
% Copyright Sohan Seth sohan.seth@hiit.fi

close all,  %clearvars -except saveFigure plotFigure maxStart;
chooseFigure = 4;
saveFigure = 0; plotFigure = 1; maxTrial = 10; % 10 for figure, 100 for NLL
maxStart = 10; % Run PAA maxStart times and take best solution

computeCostBernoulli = @(X, P)(sum(sum(-X.*log(double(P + eps)))) ...
    + sum(sum(-(1 - X).*log(1 - double(P) + eps))) );
computeCostPoisson = @(X, R)(sum(sum(- X .* log(R + eps) + R)));
computeCostStochastic = @(X, P)...
    (- sum(sum(X .* log(P + eps),2)));

%% Bernoulli data, multiple data multiple trials
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
rng default,
clear nFeatSam
d = 10;         % Dimensionality
K = 6;          % Number of archetypes
n = 100;        % Number of observations
gammaParam = 0.4; % Concentration of Dirichlet distribution for H
options = generate_options();
options.verbose = false; options.maxIter = 10000;
FONTSIZE = 6;
if plotFigure; hMain = figure('papersize',[8.5 11],'paperposition',[0 0 8 2]); colormap(flipud(gray)); end
uniqMatchesPAA = zeros(maxTrial, 1); uniqMatchesAA = zeros(maxTrial, 1);
l1List = zeros(maxTrial, 2); nllList = zeros(maxTrial, 2);
for countTrial = 1:maxTrial
    matFeatLat = rand(d, K) > 0.7;
    while any(sum(matFeatLat') == K) || any(sum(~matFeatLat') == K)
        % avoid all 1s in row for classical solution
        matFeatLat = rand(d, K) > 0.7;
    end
    matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
    matFeatSam = rand(d, n) < (matFeatLat * matLatSam);
    while any(sum(matFeatSam') == n)
        % avoid all 1s in row for classical solution
        matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
        matFeatSam = rand(d, n) < (matFeatLat * matLatSam);
    end
    
    matFeatSam(:,end) = 1; % Testing failure of classical-aa
    if plotFigure
        subplot(3,maxTrial,countTrial)
        imagesc(matFeatLat)
        if countTrial == 1
            ylabel('True','fontsize',FONTSIZE)
        end
        title(num2str(countTrial),'fontsize',FONTSIZE)
        set(gca,'fontsize',FONTSIZE,'xtick',[],'ytick',[])
    end
    
    % Run PAA and take the best solution
    matSamLatEM = cell(maxStart,1);
    matLatSamEM = cell(maxStart,1);
    objEM = cell(maxStart,1);
    for countEM = 1:maxStart
        [matSamLatEM{countEM}, matLatSamEM{countEM}, objEM{countEM}] = paa_Bernoulli(matFeatSam, K, options);
    end
    [~, ind] = min(cellfun(@(x)(x(end)), objEM));
    matSamLatEM = matSamLatEM{ind};
    matLatSamEM = matLatSamEM{ind};
    
    [~, ind] = min(pdist2(double(matFeatSam * matSamLatEM > 0.5)', double(matFeatLat)' ,'jaccard'));
    TTT = matFeatSam * matSamLatEM(:,ind);
    indNonUniq = (sum((bsxfun(@eq,ind,ind'))) ~= 1);
    TTT(:, indNonUniq) = 0; % Remove non-unique matches
    uniqMatchesPAA(countTrial) = K - sum(indNonUniq);
    if length(unique(ind)) ~= K
        disp('Index are not unique');
    end
    
    if plotFigure
        subplot(3,maxTrial,maxTrial + countTrial)
        imagesc(TTT)% > 0.5)
        text(find(sum((bsxfun(@eq,ind,ind'))) ~= 1),0*double(find(sum((bsxfun(@eq,ind,ind'))) ~= 1) > 0),'o','fontsize',FONTSIZE/2)
        if countTrial == 1
            ylabel('PAA','fontsize',FONTSIZE)
        end
        set(gca,'fontsize',FONTSIZE,'xtick',[],'ytick',[])
    end
    
    [matSamLatClassic, matLatSamClassic, obj] = classic_aa(matFeatSam, K);
    [~, ind] = min(pdist2(double(matFeatSam * matSamLatClassic > 0.5)', double(matFeatLat)' ,'jaccard'));
    TTT = matFeatSam * matSamLatClassic(:,ind);
    indNonUniq = (sum((bsxfun(@eq,ind,ind'))) ~= 1);
    TTT(:, indNonUniq) = 0; % Remove non-unique matches
    uniqMatchesAA(countTrial) = K - sum(indNonUniq);
    if length(unique(ind)) ~= K
        disp('Index are not unique');
    end
    
    if plotFigure
        subplot(3,maxTrial,2*maxTrial + countTrial)
        imagesc(TTT)% > 0.5)
        text(find(sum((bsxfun(@eq,ind,ind'))) ~= 1),0*double(find(sum((bsxfun(@eq,ind,ind'))) ~= 1) > 0),'o','fontsize',FONTSIZE/2)
        if countTrial == 1
            ylabel('AA','fontsize',FONTSIZE)
        end
        set(gca,'fontsize',FONTSIZE,'xtick',[],'ytick',[])
    end
    
    if ~plotFigure
        % Testing dataset, compute matLatSam given archetypes
        matFeatLatEM = matFeatSam * matSamLatEM;
        matFeatLatClassic = matFeatSam * matSamLatClassic;
        
        matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
        matFeatSam = rand(d, n) < (matFeatLat * matLatSam);
        
        clear matLatSamEM objEM
        options.matFeatLat = matFeatLatEM;
        for countEM = 1:maxStart
            [~, matLatSamEM{countEM}, objEM{countEM}] = paa_Bernoulli(matFeatSam, K, options);
        end
        [~, ind] = min(cellfun(@(x)(x(end)), objEM));
        matLatSamEM = matLatSamEM{ind};
        nllList(countTrial, 1) = computeCostBernoulli(matFeatSam, matFeatLatEM * matLatSamEM);
        l1List(countTrial, 1) = norm(matFeatSam - matFeatLatEM * matLatSamEM);
        options.matFeatLat = [];
        
        matLatSamClassic = classic_aa_test(matFeatSam);
        nllList(countTrial, 2) = computeCostBernoulli(matFeatSam, matFeatLatClassic * matLatSamClassic);
        l1List(countTrial, 2) = norm(matFeatSam - matFeatLatClassic * matLatSamClassic, 1);
    end
end
if saveFigure == 1
    saveas(gcf,'../Paper-NIPS-ICML/compareProbVsStandardBinary.eps','epsc')
    !epstopdf  ../Paper-NIPS-ICML/compareProbVsStandardBinary.eps
end
save compareProbVsDefault_Bernoulli nllList l1List