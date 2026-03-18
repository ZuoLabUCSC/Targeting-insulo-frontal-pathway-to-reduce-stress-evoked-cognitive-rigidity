%% photometryPriProcess

% This function reads the csv file and cut before and after movie recording  
% fit data to correct the drifting effect: two term exponential
% correct movement flactuations by 405 ch
% get zscore of 470

% Inputs: 
% ca_file: data of 405 and 470, 2nd column is the data
% event_file: data of event file recorded by fiber photometry recording system
% with TTL signal recorded in Ch1
% cafps: frame rate of ca signal, 40 Hz
% evtfps: frame rate of event file, 80 Hz

% Output: caTbl, A structure array of corrected calcium and zscored calcium

% Shaorong Ma 2023-5-18
%% data source - input 3 files
cafps = 40; 
evtfps = 80;
pathname = pwd;
ca405 = dir(fullfile(pathname,'*405.csv'));
ca470 = dir(fullfile(pathname,'*470.csv'));
caEvt = dir(fullfile(pathname,'*event*.csv'));
%% read event file find video start and stop
caful405 = readmatrix(fullfile(ca405.folder,ca405.name));
caful470 = readmatrix(fullfile(ca470.folder,ca470.name));
caF = min(length(caful405), length(caful470));
caful = [caful405(1:caF,2) caful470(1:caF,2)];
evt = readmatrix(fullfile(caEvt.folder,caEvt.name));
evt = evt(1:end, 1);
k = find(evt);
bestt = k(1);
bestp = k(end);
factor = evtfps / cafps;
castt = round ( bestt / factor);
castp = round ( bestp / factor);
%% Cut Calcium data by frame
cacut = caful(castt:castp,:);
%% Data fiting
% two-term exponential
figure;
ch = {'405' '470'};
r = cafps;
m =size(cacut,1);
time = linspace(0, m/r, m)';

for i=1:length(ch)
    f = fit(time,cacut(:,i),'exp2');
    curfit = f.a * exp(f.b*time) + f.c * exp(f.d*time); % f(x) = a*exp(b*x) + c*exp(d*x)
    F_detrend(:,i) = cacut(:,i) - curfit;
    nexttile;
    plot(time,cacut(:,i),'r'); 
    hold on; 
    plot(time,F_detrend(:,i),'b'); 
    plot(time,curfit,'g')
    legend('Ori','Fitted','CurFit');
    axis tight;
    xlabel ('Time (s)'); ylabel('Ca');
    title([ch{i} ' fit']);
end
%% correct 470 with 405
ca_corrected = F_detrend(:,2)-F_detrend(:,1);
zCa = zscore(ca_corrected);
figure;
subplot(4,1,1)
plot(time,F_detrend(:,1),'b');
hold on
legend('405'); axis tight; xlabel ('Time(s)'); ylabel('Ca');
subplot(4,1,2)
plot(time,F_detrend(:,2),'g');
legend('470'); axis tight; xlabel ('Time(s)'); ylabel('Ca');
subplot(4,1,3)
plot(time,ca_corrected,'r');
legend('Corrected'); axis tight; xlabel ('Time(s)'); ylabel('Ca');
subplot(4,1,4)
plot(time,zCa,'k')
legend('Zscore'); axis tight; xlabel ('Time(s)'); ylabel('Ca');
sgtitle('Zscored Ca correted by 405');
%% output
caTbl.ca = ca_corrected;
caTbl.zCa = zCa;