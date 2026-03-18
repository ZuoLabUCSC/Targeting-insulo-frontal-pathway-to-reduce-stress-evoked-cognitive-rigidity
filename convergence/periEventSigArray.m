function sigArray = periEventSigArray(dn,behEvtTbls)
% This functions organizes ca signals into an periEvent array
%   Inputs:
%      dn: ca signal structure
%       .ca: calcium signal matrix: nFrame x nCell
%       .time: timestamps of each frame
%       .sttFrameT: table for the start time of each session of interest
%      behEvtTbls: structure with session name as field
%           behavior event table: nEvent x nFeature

%   Outputs:
%       sigArray: a structure with subfields for signals and metadata
%           .caA: calcium signals arranged in a 4d array
%           .Events: structure for behavior variables and event times organized by trial
%           .dimension = {'neuron','framePerWindow','eventWindow','trial'};

%% organize multi-session signals
sessionNames = fieldnames(behEvtTbls);
nSession = length(sessionNames);

%% extract ca signals corresponding to sessions
% find the row (frame) index for the start of each session
sessionStartIdx = zeros(nSession,1);
for k = 1:nSession
    rowIdx = strcmp(dn.sttFrameT.session,sessionNames{k});
    sessionStartIdx(k) = find(dn.time==dn.sttFrameT.sttT(rowIdx));
end

sig = dn.ca(sessionStartIdx(1):end,:);
time = dn.time(sessionStartIdx(1):end); % in sec
fps = 1./diff(time(1:2)); % calcium frame rate ~10 Hz

%% excluding start transients
[~,nNeuron] = size(sig);

% mask recording start caused transients in ca signal
maskWin = round(3*fps); % 3 sec mask
sigM = mean(sig,2); 
movieStartIdx = find(sigM == 0); % movie starts have 0 intensity
transientFrame = movieStartIdx + (0:maskWin); % nStart x nWin
transientFrame = transientFrame(:); % linearize

sigN = sig;
sigN(transientFrame,:) = NaN;
sigN(sigN<0) = 0; % negative numbers are small and related to start transient in CNMFe

% fill NaN in sigN by interpolation
nanId = isnan(sigN(:,1));
sigN = interp1(time(~nanId),sigN(~nanId,:),time);

%% peri-event time windows
tau = 2; %2 sec
tWin = round(tau*fps); 

%% behavior variables to decode by trials: current and prior reward
% event markers in behEvtTbls
markerName = {'TrialStart','Approach','Dig','Outcome'};
nMarker = length(markerName);

% structures for behavior event times and labels, organized by trial
EvtFrame = struct;
TrialLabel = struct;
counter = 1;
for k = 1:nSession
    % single session behavior event table
    behEvtTbl = behEvtTbls.(sessionNames{k});
    uniqueTrials = unique(behEvtTbl.Trial); 
    nTrial = length(uniqueTrials);
    
    % event index
    evtIdx = [];
    evtIdx.TrialStart = contains(behEvtTbl.Behavior,'TrialStart');
    evtIdx.Approach = contains(behEvtTbl.Behavior,'Approach');
    evtIdx.Dig = contains(behEvtTbl.Behavior,'Dig') & contains(behEvtTbl.Status,'START');
    evtIdx.Outcome = contains(behEvtTbl.Behavior,'Dig') & contains(behEvtTbl.Status,'STOP');
    for i = 1:nTrial
        trialIdx = behEvtTbl.Trial == uniqueTrials(i);        
        if any(trialIdx & evtIdx.Dig) % exclude no dig (decision) trial
            
            % calculate event frames in ca movies
            for j = 1:nMarker
                markerIdx =  trialIdx & evtIdx.(markerName{j}); 
                markerIdx = find(markerIdx,1,'last'); % keep last approach
                markerFrame = behEvtTbl.caFrame(markerIdx);                     
                EvtFrame.(markerName{j})(counter,1) = markerFrame;                  
            end  

            TrialLabel.CurrRwd(counter,1) = any(contains(behEvtTbl.Behavior(trialIdx),'Eat'));
            TrialLabel.Trial(counter,1) = i;
            TrialLabel.Session(counter,1) = sessionNames(k);
            counter = counter+1;
        end
    end     
end
TrialLabel.PriorRwd = [true;TrialLabel.CurrRwd(1:end-1)]; % first trial, true for prior reward

%% calcium activity organized by nNeuron x tWin x nWin x nTrial
winName = {'preSta','posSta','preApp','posApp','preDec','posDec','preOut','posOut'};
nWin = length(winName);
totalTrial = length(TrialLabel.Trial);
caA = zeros(nNeuron,tWin,nWin,totalTrial); % ca signal array

for i = 1:totalTrial   
    for j = 1:nMarker
        caFrame = EvtFrame.(markerName{j})(i);
        caA(:,:,2*j-1,i) = sigN(caFrame-tWin : caFrame-1,:)';
        caA(:,:,2*j,i) = sigN(caFrame : caFrame+tWin-1,:)';
    end
end

% interpolate caA to 10Hz
xq = linspace(1,tWin,round(tau*10)); 
tWinIn = length(xq);
caL = reshape(permute(caA,[2 1 3 4]),tWin,nNeuron*nWin*totalTrial);
caL = interp1(1:tWin,caL,xq); 
caA = permute(reshape(caL,[tWinIn nNeuron nWin totalTrial]),[2 1 3 4]); % [nNeuron,tWin,nWin,nTrial]

%% output structure
sigArray.caA = caA;
sigArray.dimension = {'neuron','framePerWindow','eventWindow','trial'};
sigArray.TrialLabel = TrialLabel;
sigArray.winName = winName;
sigArray.winDurationSec = tau;
sigArray.fpsOrigin = fps;
sigArray.fps = 10;

end