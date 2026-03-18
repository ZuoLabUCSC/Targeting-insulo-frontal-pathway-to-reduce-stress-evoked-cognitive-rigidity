function sessionStats = PopulationDistanceDirect(sigArray,sessionName,distMethod,varargin)

%% parse inputs
winID = 2:5; % default optostimulation window, from post trial start to pre dig (2:5)
%
if ~isempty(varargin)
    idx = find(strncmpi(varargin,'window',3));
    if ~isempty(idx)
        winID = varargin{idx+1};       
    end
end

%% extract behaivoral conditions and calcium signals in selected sessions
% select sessions to analyze
nSession = length(sessionName);
trialKeepId = ismember(sigArray.TrialLabel.Session,sessionName);

% behavior conditions in selected sessions
sessionLabel = sigArray.TrialLabel.Session(trialKeepId);
CurrRwd = sigArray.TrialLabel.CurrRwd(trialKeepId);

% Ca signals in selected sessions
caA = sigArray.caA; % cell x time x window x trial

% normalize caA by cell 
caN = caA(:,:); % cell x (time x window x trial)
caN = normalize(caN,2,"zscore");
caN = reshape(caN,size(caA));

% select windows within a trial and get Ca signals in selected windows
caS = caN(:,:,winID,trialKeepId); % use only opto windows, all frame
%caS = mean(caS,2); % additional analysis: average within window

% initialize session summary stats
sessionStats = struct;
winNameS = sigArray.winName(winID);
sessionStats.winNameS = winNameS;

%% grand trial average
trialGrandMean = mean(caN,4);
trialGrandMean = trialGrandMean(:,:); % cell x (time x window)
sessionStats.trialGrandMean = trialGrandMean;

%% trial-by-trial correlation distance of cell activation patterns 
% activation pattern by trial
trialPattern = permute(caS,[4 2 3 1]); % trial x time x window x cell
trialPattern = trialPattern(:,:); % trial x (time x window x cell)

%% session summary stats
for i = 1:nSession
    
    idxS = strcmp(sessionLabel,sessionName{i}); % session index
    nTrialS = sum(idxS);
    trialPattern_i = trialPattern(idxS,:);

    %% trial-by-trial correlation distance of cell activation patterns 
    trialPD_i = pdist(trialPattern_i,distMethod); % pairwise distance
    corrMat_i = squareform(1-trialPD_i); % correlation matrix

    % mean correlation
    corrMat_i(eye(nTrialS)==1) = NaN; % set diagnal to NaN for ease of analysis later
    sessionStats.MeanCorr(i) = mean(corrMat_i(:),'omitnan'); 

    % number of errors per session
    idxC = CurrRwd(idxS);% all correct
    idxI = ~idxC; % all incorrect
    sessionStats.errorPerSession(i) = sum(idxI);

    % sequential updates dependent on prior trial outcome
    seqDist = 1-diag(corrMat_i,1); % sequential correlation distance
    idxC1 = idxC;
    idxC1(end) = []; % exclude last trial (correct), no update data, change from (end)
    updateI = mean(seqDist(idxI));
    updateC = mean(seqDist(idxC1)); 
    updateDif = updateI - updateC;
    sessionStats.updateI(i) = updateI;
    sessionStats.updateC(i) = updateC;
    sessionStats.updateDif(i) = updateDif;
    sessionStats.seqDist{i} = seqDist;

    % convergence of activity patterns to expert trials
    m = 3; % define expert trials as last m trials
    expTrial = nTrialS + 1 - (1:m); % define expert trials 
    expertPattern = mean(trialPattern_i(expTrial,:),1); % average expert trials
    novTrial = 1:nTrialS-m; % novice trials
    novicePattern = trialPattern_i(novTrial,:); % novice patterns
    corrToExpert = 1-pdist2(novicePattern,expertPattern,distMethod);% similarity to expert trials
    beta = robustfit(novTrial,corrToExpert); % regression fit for convergence rate
    
    sessionStats.corrMat{i} = corrMat_i;
    sessionStats.corrToExpert{i} = corrToExpert;
    sessionStats.convergenceRate(i) = beta(2);
    sessionStats.convergenceOffset(i) = beta(1);  
    sessionStats.trialPerSession(i) = nTrialS;

    sessionStats.expertPattern{i} = expertPattern;
    sessionStats.meanCellActPerTrial{i} = mean(trialPattern_i,2);
    sessionStats.cellActPattern{i} = trialPattern_i;

    sessionStats.trialPerSession(i) = nTrialS;
end

%% plot

%% pairwise correlation between single-trial activity patterns
if any(strncmpi (varargin,'similarity',3))
    %%
    for i = 1:nSession
        ax = nexttile;hold on;
        tmp = sessionStats.corrMat{i};
        colorRange = [-0.35 0.35];
        for k = 1:size(tmp,1)
            rectangle('Position',[k-0.5, k-0.5,1,1],...
                'FaceColor',[0.75 0.75 0.75],'EdgeColor','none');
        end
        imagesc(tmp,'AlphaData', ~isnan(tmp),colorRange);
        colormap(ax,green_white_magenta);
        axis image ij
        
        title ([sessionName{i},' Correlation of trial activity patterns']); 
        xlabel('Trial -->'); ylabel('<-- Trial');axis square
        h = colorbar; h.Label.String = 'Correlation coefficient';
    end
end       

%% Convergence to expert trials 
if any(strncmpi (varargin,'convergence',3))
    %%
     for i = 1:nSession 
        idxS = strcmp(sessionLabel,sessionName{i});
        nTrial = sessionStats.trialPerSession(i);
        novTrial = 1:nTrial-m; % novice trials
        CorrToExpert = sessionStats.corrToExpert{i};
        polyCoef = [sessionStats.convergenceRate(i),...
            sessionStats.convergenceOffset(i)];

        ax = nexttile; hold on;
        RwdIdx = CurrRwd(idxS);
        h = zeros(2,1);
        for j = 1:2
            idx = intersect(novTrial,find(RwdIdx==j-1));
            xt = novTrial(idx);
            yt = CorrToExpert(idx);
            h(j) = scatter(xt,yt,100,'o','filled','MarkerEdgeColor','none');
        end
        hL= refline(polyCoef);
        hL.LineWidth = 1;
        ylim([-0.2 0.4]); xlim([1 30]);
        title([sessionName{i}, ' trial convergence rate'], ...
            ['slope = ', num2str(sessionStats.convergenceRate(i),'%.2g')]);
        xlabel('Trial');ylabel('Correlation to the average of last 3 trials')
        legend(h,{'Incorrect','Correct'},'Location','bestoutside');
     end
end
%% IC Difference  plot 
if any(strncmpi(varargin,'ICdifference',3))
    %%
    for i = 1:nSession
        ax = nexttile;hold on;axis square;
        idxS = strcmp(sessionLabel,sessionName{i});
        RwdIdx = CurrRwd(idxS);
        seqDist = sessionStats.seqDist{i};

        colorL = {'b','r'};
        gMean = NaN(2,1);
        for k = 1:2
            pidx = RwdIdx(1:end-1) == k-1;
            gMean(k) = mean(mean(seqDist(pidx)));
            scatter(k,seqDist(pidx),50,colorL{k}, 'o','LineWidth',1);
            scatter(k,gMean(k),50,colorL{k},'*','Linewidth',1)
        end

        xlim([0 3]);xticks([1 2]);xticklabels({'I','C'});
        xlabel([]);ylabel('Distance to next trial')
        ylim([0.6 1.1]);
        legend(h,{'Incorrect','Correct'},'location','bestoutside')
        title([sessionName{i},' difference in update rate '],...,
            ['I - C =', num2str(gMean(1)-gMean(2),'%.2g')]);
    end
end    