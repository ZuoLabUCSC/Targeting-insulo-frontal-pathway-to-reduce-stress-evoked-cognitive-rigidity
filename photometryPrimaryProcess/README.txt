photometryPriProcess

Primary processing of fiber photometry recorded calcium data. Read fiber photometry recorded csv file, cut before and after behavioral video recording, correct the drifting effect caused by photo bleaching, correct movement fluctuations by 405 ch, and get zscored 470 ch data.

System requirements: 
Desktop with Windows 10 pro, Matlab R2022a, does not require any non-standard hardware, does not need to be installed.

Estimated running time: ~seconds

Input(demo data): 
need to be saved in the same folder with the codes. 

1 - 405 data: uv channel
2 - 470 data: GCaMP channel
3 - event data: recorded TTL signal driven by behavioral video recording

Instructions for use: open matlab 2022a and navigate to the folder with codes and demo data, open photometryPriProcess.m and click 'run'.

Expected output:

1 - figure with 405 and 470 ch original signals, fitted curve and fitted signals;
2 - figure with fitted 405 and 470 ch, 470 ch after corrected by 405, and zscored 470 ch after correction
3 - caTbl: a structure array with correct 470 ch (ca) and zcored 470 ch (zCa)
