%% Master_script_AST

clearvars; close all;clc;

%% load data
load('dn.mat') % calcium signal
load('behEvtTbls.mat')

%% extract peri-event signal array
sigArray = periEventSigArray(dn,behEvtTbls);

%% calculate convergence and visualization
sessionName = {'Rev','EDS'}; 
distMethod = 'correlation';
%distMethod = 'spearman'; % Choose between correlation or spearman

figure
sessionStats = PopulationDistanceDirect(sigArray,sessionName,distMethod,...
    'similarity','convergence','update','ICdiff');



