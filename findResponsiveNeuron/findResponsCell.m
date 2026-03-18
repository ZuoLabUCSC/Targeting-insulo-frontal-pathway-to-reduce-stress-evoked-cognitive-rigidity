%% find response cell for each event

% This code is to identify behavioral event (trial start, approach, decision, outcome) responsive cells

% input: sig (calcium data), behEvtTbl (behavior event table), fps: calcium frame rate 10Hz
% output: logic matrix for activated, inhibited and non-responsvie cell in each event

% Shaorong 2024-5-30

%% load input data
load('signal.mat');
load('behEvtTbl.mat')
fps = 10;    

%% extract peri-event signal array
sigArray = periEventSigArray(sig,behEvtTbl,fps);

%% use auROC to identify cell response
iter = 1000;
ROCstats = ROCanalysis(sigArray,iter);
