function sigArray = periEventSigArray(sig,behEvtTbl,fps)
% This functions organizes ca signals into a periEvent array
%   Inputs:
%       sig: calcium signal matrix: nFrame x nCell
%       behEvtTbl: behavior event table: nEvent x nFeature
%       fps: calcium signal frame rate per second
%   Outputs:
%       sigArray: a structure with subfields for signals and metadata
%           .caA: calcium signals arranged in a 4d array
%           .Events: structure for behavior variables and event times organized by trial
%           .dimension = {'neuron','framePerWindow','eventWindow','trial'};

%% normalize ca activity by cell
[~,nNeuron] = size(sig);

% mask movie start caused transients in ca signal
maskWin = round(2*fps); % 2 sec mask
sigM = mean(sig,2); 
movieStartIdx = find(sigM == 0); % movie starts have 0 intensity, this also include some Ca intensity == 0 frames which cause problem
transientFrame = movieStartIdx + (0:maskWin);
transientFrame = transientFrame(:);

% zscore ca sig
sigN = sig;
sigN(transientFrame,:) = NaN;
sigN = normalize(sigN,1,'zscore');

%% peri-event time windows
tau = 2; %2 sec
tWin = round(tau*fps); %  2s windoow

%% event marker index in behEvtTbl
trialStartIdx = contains(behEvtTbl.Behavior,'TrialStart');
approachIdx = contains(behEvtTbl.Behavior,'Approach');
digIdx = contains(behEvtTbl.Behavior,'Dig') & contains(behEvtTbl.Status,'START');
outcomeIdx = (contains(behEvtTbl.Behavior,'Eat') & contains(behEvtTbl.Status,'START'))|(contains(behEvtTbl.Behavior,'Leave'));

%% calculate behavior events in calcium movie frames

%% behavior variables to decode by trials: current and prior reward
uniqueTrials = unique(behEvtTbl.Trial);
nTrial = length(uniqueTrials);

% structure for behavior variables and event times organized by trial
Events = struct;
for i = 1:nTrial
    trialIdx = behEvtTbl.Trial == uniqueTrials(i);
    Events.Frame.start(i,1) = behEvtTbl.caFrame(trialIdx & trialStartIdx);
    Events.Frame.approach(i,1) = max(behEvtTbl.caFrame(trialIdx & approachIdx)); % the last apprachFrame
    Events.Frame.dig(i,1) = behEvtTbl.caFrame(trialIdx & digIdx); 
    Events.Frame.outcome(i,1) = behEvtTbl.caFrame(trialIdx & outcomeIdx); 
    Events.behVar.outcome(i,1) = any(contains(behEvtTbl.Behavior(trialIdx),'Eat'));
end
Events.behVar.priorRwd = [false;Events.behVar.outcome(1:end-1)];

%% calcium activity organized by nNeuron x tWin x nWin x nTrial
markerName = fieldnames(Events.Frame);
nMarker = length(markerName);
winName = {'preSta','posSta','preApp','posApp','preDec','posDec','preOut','posOut'};
nWin = length(winName);

caA = zeros(nNeuron,tWin,nWin,nTrial); % ca signal array
for i = 1:nTrial    
    for j = 1:nMarker
        caFrame = Events.Frame.(markerName{j})(i);
        caA(:,:,2*j-1,i) = sigN(caFrame-tWin : caFrame-1,:)';
        caA(:,:,2*j,i) = sigN(caFrame : caFrame+tWin-1,:)';
    end
end

% interpolate caA to 10Hz
xq = linspace(1,tWin,round(tau*10)); 
tWinIn = length(xq);
caL = reshape(permute(caA,[2 1 3 4]),tWin,nNeuron*nWin*nTrial);
caL = interp1(1:tWin,caL,xq); %
caA = permute(reshape(caL,[tWinIn nNeuron nWin nTrial]),[2 1 3 4]); % [nNeuron,tWin,nWin,nTrial]

% output structure
sigArray.caA = caA;
sigArray.Events = Events;
sigArray.dimension = {'neuron','framePerWindow','eventWindow','trial'};
sigArray.winName = winName;
sigArray.winDurationSec = tau;
sigArray.fps = 10;
sigArray.ca = sigN;
end