function [ROCstats] = ROCanalysis(sigArray,iter)
% This function is used to calculate the auROC of behavioral events
% observations: averaged ca in 2s window, pre and post
% label: lagical with pre 0 post 1
% output: auROC
% 
% Shaorong 2024-5-29
%% Get auROC from data
caA = sigArray.caA;
[nNeuron,~,nWin,nTrial] = size(caA); % 2nd dimension is frame, averaged below
lbl = logical([zeros(nTrial,1);ones(nTrial,1)]);
auROCarray = nan(nNeuron,nWin/2);

for j=1:nWin/2
    caPre = squeeze(mean(caA(:,:,2*j-1,:),2,"omitnan")); % neuron X trial with mean of Ca in a window
    caPost = squeeze(mean(caA(:,:,2*j,:),2,"omitnan"));
    for neuron=1:nNeuron
        caC = [caPre(neuron,:)';caPost(neuron,:)'];
        input = table(caC,lbl); % observation + label
        input.Properties.VariableNames = {'Ca' 'Label'};
        [~,~,~,auROC] = perfcurve(input.Label,input.Ca,true);
        auROCarray(neuron,j) = auROC;       
    end
end
%% Permutation
auROCperm = nan(iter,1);
fps = sigArray.fps;
winT = sigArray.winDurationSec;
ca = sigArray.ca;
nFrame = size(ca,1);
for i=1:iter
    permF = randperm(nFrame-winT*fps*2, nTrial)+winT*fps;
    permPre = NaN(nNeuron,nTrial);
    permPost = NaN(nNeuron,nTrial);
    for j=1:nTrial
        permPre(:,j) = mean(ca(permF(j)-winT*fps : permF(j)-1, :),"omitnan")';
        permPost(:,j) = mean(ca(permF(j)+1 : permF(j)+winT*fps, :),"omitnan")';
    end   
    for neuron = 1:nNeuron
        permC = [permPre(neuron,:)';permPost(neuron,:)'];
        inputPerm = table(permC,lbl);
        inputPerm.Properties.VariableNames = {'permCa' 'Label'};
        [~,~,~,aucPerm] = perfcurve(inputPerm.Label,inputPerm.permCa,true);
        auROCperm(i) = aucPerm;
    end
end
p975 = prctile(auROCperm,97.5);
p25 = prctile(auROCperm,2.5);
%% Identify acitve or inhibit cell
act = auROCarray > p975;
inh = auROCarray < p25;
non = ~act & ~inh;
rep.act = act;
rep.inh = inh;
rep.non = non;
ROCstats.auROCarray = auROCarray;
ROCstats.auROCperm = auROCperm;
ROCstats.response = rep;