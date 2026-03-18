findResponseCell

To identify neurons that responsive to specific behavioral events (trial start, approach, decision, outcome). 

This script calls two functions:
periEventSigArray: organizes calcium signals into a periEvent Array
ROCstats: calculate the area under the ROC curve (auROC) of specific behavioral event, and identify is the neuron is activated, inhibited or not responsive to the event


Estimated running time: ~mins

Input (demo data):
need to be saved in the same folder with the codes. 

1. sig: calcium data after preprocessing, organized by frame X neuron, see example signal.mat
2. behEvtTbl: behavior table generated after annotation, see example behEvtTbl.mat

Instructions for use: open matlab R2022a and navigate to the folder with codes and demo data, open findResponseCell.m and click 'run'.

Expected output:

ROCstats: a structure array with 3 fields (auROCarray, auROCperm, response)
 auROCarray - stores the auROC results, organized by neuron X event
 auROCperm - stores the auROC result of 1000 permutation of shuffling time, these data are used to generate null distribution
 response - is a structure array stores 3 logical arrays of activated neurons, inhibited neurons and non-responsive neurons in each event, each organized by neuron X event

   