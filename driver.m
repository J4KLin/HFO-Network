%% Run full analysis for one patient

addpath('helpers');

patientID = 28;


umid = strcat("UMHS-00",string(patientID));

hfodatafolderpath = fullfile(pwd,'\Data\HFOEEG');
bkgdatafilepath = fullfile(pwd,'\Data\BKGEEG',strcat(umid,'.mat'));
patientfilepath = fullfile(pwd,'\Data\PatientData',strcat(umid,'.mat'));

%% Run HFO EEG Cross Correlations
runCCHFO_wKDEConfidence(umid,hfodatafolderpath,patientfilepath);

%% Run Background EEG Cross Correlations
runCCBKG_wKDEConfidence(bkgdatafilepath,'rmsdata');
        
%% Compute Functional Connectivity
runFCN_wKDEZTest(umid,hfodatafolderpath,bkgdatafilepath,patientfilepath);

%% Compute Connectivity Matricies and Centralities
runCentralities(patientfilepath,bkgdatafilepath);

%% Split centrality rank into percentile groups by electrode types for all features
getPercentileRanks(patientfilepath);

%% Calculate Critical Resection Percentage (CReP) for all features
getCReP(patientfilepath);
