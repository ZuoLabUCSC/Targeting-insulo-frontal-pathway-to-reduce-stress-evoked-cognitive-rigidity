Master_script_AST

Calculate the trial-by-trial convergence rate and update rate of attentional set-shifting task (AST) by using miniscope calcium data and behavior annotation data.

This script calls two functions:
perEventSigArray: organizes calcium signals into an periEvent array
PopulationDistanceDirect: calculate the trial-by-trial convergence and update

System requirements: 
Desktop with Windows 10 pro, Matlab R2022a, need Matlab Statistics and Machine Learning Toolbox. Does not require any non-standard hardware, does not need to be installed.

Estimated running time: ~seconds

Input(demo data): 
need to be saved in the same folder with the codes. 

1 - dn.mat, preprocessed miniscope calcium imaging data stored as a structure array with 5 fields: frame, cell, ca, time and sttFrameT
2 - behEvtTbls.mat, a structure array of behavioral annotations formatted in table array.

Instructions for use: data need to be formatted following the demo data structure, open matlab 2022a and navigate to the folder with codes and demo data, open Master_script_AST.m and click 'run'.

Expected output:
figure with 6 panels: 
1 - Rev session correlation matrices of pairwise comparisons of neuron population activity pattern; 
2 - EDS session correlation matrices of pairwise comparisons of neuron population activity pattern;
3 - Rev session convergence rate: the correlation between each trial's population activity pattern and the final pattern (average of the last three trials), the slope of the regression line relating the correlation coefficients to the trial numver indicates the rate of convergence; 
4 - EDS session convergence rate: the correlation between each trial's population activity pattern and the final pattern (average of the last three trials), the slope of the regression line relating the correlation coefficients to the trial numver indicates the rate of convergence; 
5 - Rev session update rate: sequential update quantified as the distance between the population activity pattern (1-correlation coefficient), and grouped by following incorrect vs correct trials.
6 - EDS session update rate: sequential update quantified as the distance between the population activity pattern (1-correlation coefficient), and grouped by following incorrect vs correct trials.
