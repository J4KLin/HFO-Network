%%%
%(For background data) Compute the maximum cross correlation between 
%every pair of electrodes' EEG data and kernel density of the distribution 
%of lags of the maximum cross correlation lags
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% datafolderpath:       folder directory of patients background EEG
% datatype:             type of EEG data to process (default: rmsdata);
%                       set filtdata for 80-500 Hz filtered EEG data
%%%
function runCCBKG_wKDEConfidence(datafilepath,datatype)
%%
% umid = 'UMHS-0028'; 
% datafilepath = fullfile(pwd,'\Data\BKGEEG',strcat(umid,'.mat'));
% datatype = 'rmsdata';

%%
load(datafilepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #JUST TO GET ELECCT
elecct = PatientInfo.Electrodes.nChan;
    
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

%baseline normalization info
datainfo.(datatype).maxcclags = {};
datainfo.(datatype).maxccamps = {};

datainfo.(datatype).alpha0 = alpha0;
datainfo.(datatype).edges = edges;

datainfo.(datatype).lagkdemodelsCB = {};
datainfo.(datatype).lagkdepdfCB = {};
datainfo.(datatype).lagkdesupCB = {};
    
%% Run analysis
[~,dlengths] = cellfun(@size,data.(datatype));
if(~all(dlengths == dlengths(1))); error('DATALENGTH NOT ALL THE SAME'); end
midpointidx = floor(dlengths(1)/2)+1;
mainidxoffset = (mainbound.*FS./1000)+midpointidx;

cursampct = length(data.(datatype));
for eleci = 1:elecct
    for elecj = eleci+1:elecct
        sprintf('%s Elec: %d vs Elec: %d\n',PatientInfo.umid,eleci,elecj)
        curmaxampsij = nan(cursampct,1);
        curmaxlagsij = nan(cursampct,1);
            
        for sampi = 1:cursampct
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
        datainfo.(datatype).maxccamps{elecj,eleci} = curmaxampsij;
        datainfo.(datatype).maxcclags{eleci,elecj} = curmaxlagsij;
        datainfo.(datatype).maxcclags{elecj,eleci} = -curmaxlagsij;
        
        %Remove zero lag radius
        curmaxlagsij(abs(curmaxlagsij) <= zerolagradius) = nan;
        curmaxlagsij = curmaxlagsij(~isnan(curmaxlagsij));
        
        if(~isempty(curmaxlagsij))
            %KDE Confidence Band
            [kdepdfCB,supCB,kdemodelCB] = kde_wBTCBunbiased(curmaxlagsij,alpha0,edges);
            datainfo.(datatype).lagkdemodelsCB{eleci,elecj} = kdemodelCB;
            datainfo.(datatype).lagkdepdfCB{eleci,elecj} = kdepdfCB;
            datainfo.(datatype).lagkdepdfCB{elecj,eleci} = kdepdfCB(end:-1:1);
            datainfo.(datatype).lagkdesupCB{eleci,elecj} = supCB;
            datainfo.(datatype).lagkdesupCB{elecj,eleci} = supCB(end:-1:1);
        end
    end
end

%% Save Results
save(datafilepath,'datainfo','-append','-v7.3')








            