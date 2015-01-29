% Qualitative assessment of probabilistic archetypal analysis and standard
% archetypal analysis by Cutler and Breiman on different observation models
%
% Additional info: Figures 4, 5 and 6 in article Probabilistic Archetypal
% Analysis by Seth and Eugster in Machine Learning Journal
%
% Copyright Sohan Seth sohan.seth@hiit.fi

close all,  %clearvars -except saveFigure plotFigure maxStart;
chooseFigure = 5;
saveFigure = 0; plotFigure = 1; maxTrial = 10; % 10 for figure, 100 for NLL
maxStart = 10; % Run PAA maxStart times and take best solution

computeCostBernoulli = @(X, P)(sum(sum(-X.*log(double(P + eps)))) ...
    + sum(sum(-(1 - X).*log(1 - double(P) + eps))) );
computeCostPoisson = @(X, R)(sum(sum(- X .* log(R + eps) + R)));
computeCostStochastic = @(X, P)...
    (- sum(sum(X .* log(P + eps),2)));

%% Poisson data
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
rng default,
clear nFeatSam
d = 12;         % Dimensionality
K = 6;          % Number of archetypes
n = 100;        % Number of observations (500)
options.eps = 10^-8;
gammaParam = 0.4; % Concentration of Dirichlet distribution for H
options = generate_options();
options.verbose = false; options.maxIter = 10000;
FONTSIZE = 6;
if plotFigure, hMain = figure('papersize',[8.5 11],'paperposition',[0 0 8 2]); colormap(flipud(gray)); end
uniqMatchesPAA = zeros(maxTrial, 1); uniqMatchesAA = zeros(maxTrial, 1);
l1List = zeros(maxTrial, 2); nllList = zeros(maxTrial, 2);
for countTrial = 1:maxTrial
    matFeatLat = 0.0*ones(d, K);
    matFeatLat(:,K) = ceil(10*rand(d,1));
    for countK = 2:K-1
        matFeatLat((countK-1)*2 + 1:countK* 2,countK) = [ceil(rand*10), ceil(rand*10)]; %[3 3];
    end
    
    matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
    matFeatSam = poissrnd(matFeatLat * matLatSam);
    
    if plotFigure
        subplot(3,maxTrial,countTrial)
        imagesc(matFeatLat)
        if countTrial == 1
            ylabel('True','fontsize',FONTSIZE)
        end
        title(num2str(countTrial),'fontsize',FONTSIZE)
        set(gca,'fontsize',FONTSIZE,'xtick',[],'ytick',[])
    end
    
    
    % Run PAA and take best solution
    matSamLatEM = cell(maxStart,1);
    matLatSamEM = cell(maxStart,1);
    objEM = cell(maxStart,1);
    for countEM = 1:maxStart
        [matSamLatEM{countEM}, matLatSamEM{countEM}, objEM{countEM}] = paa_Poisson(matFeatSam, K, options);
    end
    [~, ind] = min(cellfun(@(x)(x(end)), objEM));
    matSamLatEM = matSamLatEM{ind};
    matLatSamEM = matLatSamEM{ind};
    
    [~, ind] = min(pdist2((matFeatSam * matSamLatEM)', matFeatLat' ,'cityblock'));
    TTT = matFeatSam * matSamLatEM(:,ind);
    indNonUniq = (sum((bsxfun(@eq,ind,ind'))) ~= 1);
    TTT(:, indNonUniq) = 0; % Remove non-unique matches
    uniqMatchesPAA(countTrial) = K - sum(indNonUniq);
    if length(unique(ind)) ~= K
        disp('Index are not unique');
    end
    
    if plotFigure
        subplot(3,maxTrial,maxTrial + countTrial)
        imagesc(TTT)
        text(find(sum((bsxfun(@eq,ind,ind'))) ~= 1),0*double(find(sum((bsxfun(@eq,ind,ind'))) ~= 1) > 0),'o','fontsize',FONTSIZE/2)
        if countTrial == 1
            ylabel('PAA','fontsize',FONTSIZE)
        end
        set(gca,'fontsize',FONTSIZE,'xtick',[],'ytick',[])
    end
    
    [matSamLatClassic, matLatSamClassic, obj] = classic_aa(matFeatSam, K);
    [~, ind] = min(pdist2((matFeatSam * matSamLatClassic)', matFeatLat' ,'cityblock'));
    TTT = matFeatSam * matSamLatClassic(:,ind);
    indNonUniq = (sum((bsxfun(@eq,ind,ind'))) ~= 1);
    TTT(:, indNonUniq) = 0; % Remove non-unique matches
    uniqMatchesAA(countTrial) = K - sum(indNonUniq);
    if length(unique(ind)) ~= K
        disp('Index are not unique');
    end
    
    if plotFigure
        subplot(3,maxTrial,2*maxTrial + countTrial)
        imagesc(TTT)
        text(find(sum((bsxfun(@eq,ind,ind'))) ~= 1),0*double(find(sum((bsxfun(@eq,ind,ind'))) ~= 1) > 0),'o','fontsize',FONTSIZE/2)
        if countTrial == 1
            ylabel('AA','fontsize',FONTSIZE)
        end
        set(gca,'fontsize',FONTSIZE,'xtick',[],'ytick',[])
        drawnow
    end
    
    if ~plotFigure
        % Testing dataset, compute matLatSam given archetypes
        matFeatLatEM = matFeatSam * matSamLatEM;
        matFeatLatClassic = matFeatSam * matSamLatClassic;
        
        matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
        matFeatSam = poissrnd(matFeatLat * matLatSam);
        
        clear matLatSamEM objEM
        options.matFeatLat = matFeatLatEM;
        for countEM = 1:maxStart
            [~, matLatSamEM{countEM}, objEM{countEM}] = paa_Poisson(matFeatSam, K, options);
        end
        [~, ind] = min(cellfun(@(x)(x(end)), objEM));
        matLatSamEM = matLatSamEM{ind};
        nllList(countTrial, 1) = computeCostPoisson(matFeatSam, matFeatLatEM * matLatSamEM);
        l1List(countTrial, 1) = norm(matFeatSam - matFeatLatEM * matLatSamEM);
        options.matFeatLat = [];
        
        matLatSamClassic = classic_aa_test(matFeatSam);
        nllList(countTrial, 2) = computeCostPoisson(matFeatSam, matFeatLatClassic * matLatSamClassic);
        l1List(countTrial, 2) = norm(matFeatSam - matFeatLatClassic * matLatSamClassic, 1);
    end
end
if saveFigure == 1
    saveas(gcf,'compareProbVsStandardPoisson.eps','epsc')
    !epstopdf  compareProbVsStandardPoisson.eps
end
save compareProbVsDefault_Poisson nllList l1List