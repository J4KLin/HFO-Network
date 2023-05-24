%%%
%(For HFO centered data) Compute the maximum cross correlation between 
%every pair of electrodes' EEG data and kernel density of the distribution 
%of lags of the maximum cross correlation lags
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% umid:                 patient identifier (ie. UMHS-0028)
% datafolderpath:       folder directory of patients HFO EEG
% savefilepath:         savepath
% datatype:             type of EEG data to process (default: rmsdata);
%                       set filtdata for 80-500 Hz filtered EEG data
%%%
function runCCHFO_wKDEConfidence(umid,datafolderpath,savefilepath,datatype)
%%
% umid = 'UMHS-0028'; 
% datatype = 'rmsdata';
% datafolderpath = fullfile(pwd,'..','Data','HFOEEG');
% savefilepath = fullfile(pwd,'..','Data','PatientData');
%%%

%%
load(fullfile(datafolderpath,strcat(umid,'-elec1.mat')),'PatientInfo');

elecct = PatientInfo.Electrodes.nChan;

%%%%Make sure data folder exist and has all electrode files%%%%%%%%%%%%%%%%
if(~isfolder(datafolderpath)); error(strcat(umid, " DATA DOES NOT EXIST")); end

datafolderdir = dir(datafolderpath);
elecRegexp =strcat('^',umid,'-elec\d*.mat$');
elecfiletokens = regexp({datafolderdir.name},elecRegexp);
elecfilect = sum([elecfiletokens{:}]);

if(elecfilect ~= elecct)
    error(strcat("ACTUAL ", string(elecfilect), " EXPECTED ", string(elecct),...
        "elec.mat FILES for ", umid));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%ANALYSIS WINDOW PARAMETERS
MINSAMPREQ = 5;                 %ignore electrodes with less than minimal sample requirements
mainbound = [-100 100];         %(ms) boundaries for cross correlation
zerolagradius = 1;              %(ms) radius around zero lag CC to ignore

alpha0 = .05;
binwidth = .3;      %ms
edges = (binwidth/2:binwidth:100);
edges = [-edges(end:-1:1),edges];

%PROCESSING PARAMS
if(~exist('datatype','var')); datatype = 'rmsdata'; end
    
%SETTING VARS
if(~exist('datainfo','var')); datainfo = struct(); end
datainfo.allHFOi = {};                    %all actual HFO sample indicies used per electrode
datainfo.hfosamplects = zeros(elecct);    %number of samples used per electrode

%baseline normalization info
datainfo.(datatype).maxcclags = {};
datainfo.(datatype).maxccamps = {};

datainfo.(datatype).alpha0 = alpha0;
datainfo.(datatype).edges = edges;

datainfo.(datatype).lagkdemodelsCB = {};
datainfo.(datatype).lagkdepdfCB = {};
datainfo.(datatype).lagkdesupCB = {};

%% Run analysis
for eleci = 1:elecct
    fprintf('%s: Elec_ %i\n', umid, eleci);
    load(fullfile(char(datafolderpath),char(strcat(umid,"-elec",string(eleci),".mat"))),...
        'FS','data','elecROI','hfoI');
    
    curhfoi = 1:length(hfoI);
    
    datainfo.allHFOi{eleci,1} = hfoI(curhfoi);
    cursampct = length(curhfoi);
    datainfo.hfosamplects(eleci,1) = cursampct;

    [~,dlengths] = cellfun(@size,data.(datatype));
    if(~all(dlengths == dlengths(1))); error('DATALENGTH NOT ALL THE SAME'); end
    midpointidx = floor(dlengths(1)/2)+1;
    mainidxoffset = (mainbound.*FS./1000)+midpointidx;

    %Run CC only if there is enough samples
    if(cursampct >= MINSAMPREQ)
        for elecj = 1:elecct
            if(eleci ~= elecj)
                sprintf('%s Elec: %d vs Elec: %d\n',PatientInfo.umid,eleci,elecj)
                curmaxampsij = nan(cursampct,1);
                curmaxlagsij = nan(cursampct,1);
        
                for i = 1:cursampct
                    sampi = curhfoi(i);
                    curdata = data.(datatype){sampi}([eleci,elecj],:);
                    
                    if(isempty(find(isnan(curdata),1)))     %Valid data pair only
                        [CCamp, CClag] = xcorr(curdata(1,mainidxoffset(1):mainidxoffset(2)),...
                            curdata(2,mainidxoffset(1):mainidxoffset(2)),'coeff');
                        [~,maxidx] = max(abs(CCamp));        %max CC
                        curmaxampsij(sampi) = CCamp(maxidx);
                        curmaxlagsij(sampi) = CClag(maxidx);                
                    end
                end
                curmaxlagsij = curmaxlagsij .* (1000/FS);
                
                datainfo.(datatype).maxccamps{eleci,elecj} = curmaxampsij;
                datainfo.(datatype).maxcclags{eleci,elecj} = curmaxlagsij;

                %Remove zero lag radius
                curmaxlagsij(abs(curmaxlagsij) <= zerolagradius) = nan;
                curmaxlagsij = curmaxlagsij(~isnan(curmaxlagsij));

                if(~isempty(curmaxlagsij))
                    %KDE Confidence Band
                    [kdepdfCB,supCB,kdemodelCB] = kde_wBTCBunbiased(curmaxlagsij,alpha0,edges);
                    datainfo.(datatype).lagkdemodelsCB{eleci,elecj} = kdemodelCB;
                    datainfo.(datatype).lagkdepdfCB{eleci,elecj} = kdepdfCB;
                    datainfo.(datatype).lagkdesupCB{eleci,elecj} = supCB;
                end

            end
        end
    
    end
end

%% Save Results
save(fullfile(char(savefilepath)),'PatientInfo','elecct','datainfo')














