% Qualitative assessment of probabilistic archetypal analysis and standard
% archetypal analysis by Cutler and Breiman on different observation models
%
% Additional info: Figures 4, 5 and 6 in article Probabilistic Archetypal
% Analysis by Seth and Eugster in Machine Learning Journal
%
% Copyright Sohan Seth sohan.seth@hiit.fi

close all,  %clearvars -except saveFigure plotFigure maxStart;
chooseFigure = 4; % 4 or 5 or 6;
saveFigure = 0; plotFigure = 0; maxTrial = 100; % 10 for figure
maxStart = 10; % Run PAA maxStart times and take best solution

computeCostBernoulli = @(X, P)(sum(sum(-X.*log(double(P + eps)))) ...
    + sum(sum(-(1 - X).*log(1 - double(P) + eps))) );
computeCostPoisson = @(X, R)(sum(sum(- X .* log(R + eps) + R)));
computeCostStochastic = @(X, P)...
    (- sum(sum(X .* log(P + eps),2)));

switch chooseFigure
    case 4 % Bernoulli data, multiple data multiple trials
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
        if saveFigure == 1
            saveas(gcf,'../Paper-NIPS-ICML/compareProbVsStandardBinary.eps','epsc')
            !epstopdf  ../Paper-NIPS-ICML/compareProbVsStandardBinary.eps
        end
        save compareProbVsDefault_Bernoulli nllList l1List
        
    case 5 % Poisson data
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
        if saveFigure == 1
            saveas(gcf,'compareProbVsStandardPoisson.eps','epsc')
            !epstopdf  compareProbVsStandardPoisson.eps
        end
        save compareProbVsDefault_Poisson nllList l1List
        
    case 6 % Multinomial data
        rng default, close all,
        MARKERSIZE = 4; FONTSIZE = 6;
        if plotFigure; hMain = figure('papersize',[8.5 11],'paperposition',[0 0 8 2]); end;
        setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
        
        d = 3;          % Dimensionality
        K = 5;          % Number of archetypes
        n = 100;        % Number of observations (500)
        colorList = summer(maxTrial);
        gammaParam = 0.5; % Concentration of Dirichlet distribution for H
        if d ~= 3
            fprintf('Only valid for d = 3.\n')
        end
        matFeatLat = convexCircle(K);
        matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
        matFeatSam = matFeatLat * matLatSam;
        nFeatSam = mnrnd(ceil(1000*rand(n,1) + 1000), matFeatSam')';
        
        options = generate_options();
        options.maxIter = 1000;
        options.verbose = false;
        options.eps = 10^-6;
        options.display = 0;
        
        matFeatSam = bsxfun(@rdivide, nFeatSam, sum(nFeatSam));
        l1List = zeros(maxTrial, 2); nllList = zeros(maxTrial, 2);
        for countTrial = 1:maxTrial
            matSamLatEM = cell(maxStart,1);
            matLatSamEM = cell(maxStart,1);
            objEM = cell(maxStart,1);
            for count = 1:maxStart
                [matSamLatEM{count}, matLatSamEM{count}, objEM{count}] = paa_stochastic(nFeatSam, K, options);
            end
            [~, ind] = min(cellfun(@(x)(x(end)), objEM));
            matSamLatEM = matSamLatEM{ind};
            matLatSamEM = matLatSamEM{ind};
            
            if plotFigure
                figure(hMain)
                subplot(131),
                if countTrial == 1
                    set(gca,'fontsize',FONTSIZE), box on
                    hold on
                    plot(matFeatSam(1,:),matFeatSam(2,:),'o','markersize',MARKERSIZE,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','b')
                end
                if countTrial == maxTrial
                    %ind = convhull(matFeatSam(1:2,:)');
                    %hold on, plot(matFeatSam(1,ind), matFeatSam(2,ind), '-', 'color', [0.5 0.5 0.5]); hold off
                    title('(a) EM solution in the latent space')
                    grid on
                end
                matFeatLatEM = matFeatSam * matSamLatEM;
                plot(matFeatLatEM(1,:),matFeatLatEM(2,:),'ro','markersize',MARKERSIZE+2,'markerfacecolor',colorList(countTrial,:),'markeredgecolor',[0.5 0.5 0.5])
                ind = convhull(matFeatLatEM(1:2,:)');
                plot(matFeatLatEM(1,ind), matFeatLatEM(2,ind), '-', 'color', [0.5 0.5 0.5]);
            end
            
            [matSamLatClassic, matLatSamClassic, objClassic] = classic_aa(nFeatSam, K);
            
            if plotFigure
                figure(hMain)
                subplot(132),
                if countTrial == 1
                    set(gca,'fontsize',FONTSIZE), box on
                    hold on
                    plot(nFeatSam(1,:),nFeatSam(2,:),'o','markersize',MARKERSIZE,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','b')
                end
                if countTrial == maxTrial
                    title('(b) Standard solution in the observation space')
                    grid on
                end
                matFeatLatClassic = nFeatSam * matSamLatClassic;
                plot(matFeatLatClassic(1,:),matFeatLatClassic(2,:),'ro','markersize',MARKERSIZE+2,'markerfacecolor',colorList(countTrial,:),'markeredgecolor',[0.5 0.5 0.5])
                ind = convhull(matFeatLatClassic(1:2,:)');
                plot(matFeatLatClassic(1,ind), matFeatLatClassic(2,ind), '-', 'color', [0.5 0.5 0.5]);
            end
            
            if plotFigure
                figure(hMain)
                subplot(133),
                if countTrial == 1
                    set(gca,'fontsize',FONTSIZE), box on
                    hold on
                    plot(matFeatSam(1,:),matFeatSam(2,:),'o','markersize',MARKERSIZE,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','b')
                end
                if countTrial == maxTrial
                    title('(b'') Standard solution projected to latent space')
                    grid on
                end
                matFeatLatClassic = matFeatSam * matSamLatClassic;
                hold on, plot(matFeatLatClassic(1,:),mamatFeatLatClassictFeatLat(2,:),'ro','markersize',MARKERSIZE+2,'markerfacecolor',colorList(countTrial,:),'markeredgecolor',[0.5 0.5 0.5])
                ind = convhull(matFeatLatClassic(1:2,:)');
                plot(matFeatLatClassic(1,ind), matFeatLatClassic(2,ind), '-', 'color', [0.5 0.5 0.5]);
            end
            
            % Testing dataset, compute matLatSam given archetypes
            matFeatLatEM = matFeatSam * matSamLatEM;
            matFeatLatClassic = matFeatSam * matSamLatClassic;
            
            matFeatSam = matFeatLat * matLatSam;
            nFeatSam = mnrnd(ceil(1000*rand(n,1) + 1000), matFeatSam')';
            matFeatSam = bsxfun(@rdivide, nFeatSam, sum(nFeatSam)); % empirical
            
            clear matLatSamEM objEM
            options.matFeatLat = matFeatLatEM;
            for countEM = 1:maxStart
                [~, matLatSamEM{countEM}, objEM{countEM}] = paa_stochastic(nFeatSam, K, options);
            end
            [~, ind] = min(cellfun(@(x)(x(end)), objEM));
            matLatSamEM = matLatSamEM{ind};
            nllList(countTrial, 1) = computeCostStochastic(nFeatSam, matFeatLatEM * matLatSamEM);
            l1List(countTrial, 1) = norm(matFeatSam - matFeatLatEM * matLatSamEM);
            options.matFeatLat = [];
            
            matLatSamClassic = classic_aa_test(nFeatSam);
            nllList(countTrial, 2) = computeCostStochastic(nFeatSam, matFeatLatClassic * matLatSamClassic);
            l1List(countTrial, 2) = norm(matFeatSam - matFeatLatClassic * matLatSamClassic, 1);
        end
        if plotFigure
            legend('samples','archetypes','location','NorthEast')
        end
        if saveFigure == 1
            saveas(gcf,'compareProbVsStandard.eps','epsc')
            !epstopdf  compareProbVsStandard.eps
        end
        save compareProbVsDefault_stochastic nllList l1List
end

%% Bernoulli data, one data multiple trials
% rng default, clear all, close all
% clear nFeatSam
% d = 10;         % Dimensionality
% K = 5;          % Number of archetypes
% n = 500;        % Number of observations
% maxTrial = 10;  % Number of trials
% gammaParam = 0.4; % Concentration of Dirichlet distribution for H
% options = generate_options();
% options.verbose = 1; options.maxIter = 10000;
% matFeatLat = rand(d, K) > 0.7;
% [d, K] = size(matFeatLat); n = 100;
% matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
% matFeatSam = rand(d, n) < (matFeatLat * matLatSam);
% hMain = figure('papersize',[8.5 11],'paperposition',[0 0 8 2]);
% subplot(2,maxTrial+1,1)
% imagesc(matFeatLat), colormap gray
% title('True')
% for countTrial = 1:maxTrial
%     subplot(2,maxTrial+1,countTrial+1)
%     [matSamLatEM, matLatSamEM, objEM] = paa_Bernoulli(matFeatSam, K, options);
%     [~, ind] = min(pdist2(double(matFeatSam * matSamLatEM > 0.5)', double(matFeatLat)' ,'jaccard'));
%     if length(unique(ind)) ~= K
%         disp('Index are not unique');
%         imagesc((matFeatSam * matSamLatEM))% > 0.5)
%         title([num2str(countTrial),' unmatched'])
%     else
%         imagesc((matFeatSam * matSamLatEM(:,ind)))% > 0.5)
%         title([num2str(countTrial),' matched'])
%     end
%
%     subplot(2,maxTrial+1,maxTrial + 1 + countTrial+1)
%     [matSamLatClassic, matLatSamClassic, obj] = classic_aa(matFeatSam, K);
%     [~, ind] = min(pdist2(double(matFeatSam * matSamLatClassic > 0.5)', double(matFeatLat)' ,'jaccard'));
%     if length(unique(ind)) ~= K
%         disp('Index are not unique');
%         imagesc((matFeatSam * matSamLatClassic))% > 0.5)
%         title([num2str(countTrial),' unmatched'])
%     else
%         imagesc((matFeatSam * matSamLatClassic(:,ind)))% > 0.5)
%         title([num2str(countTrial),' matched'])
%     end
% end
%myunique(bin2dec(num2str(matFeatSam')))
%myunique(bin2dec(num2str(matFeatLat')))