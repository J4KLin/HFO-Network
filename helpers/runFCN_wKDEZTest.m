%%%
%Compute functional connectivity matrix from cross correlation data
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% umid:                 patient identifier (ie. UMHS-0028)
% datafolderpath:       folder directory of patients HFO correlation data
% bkgdatafilepath:      file path of patient background correlation data
% analysisfilepath:     file path of patient analysis data
%
%%%
function runFCN_wKDEZTest(umid,datafolderpath, bkgdatafilepath,analysisfilepath)
%Addition to previous runCC_FCM@2: ignore CC results near 0 lag
% umid = "UMHS-0022";    
% datafolderpath = fullfile(pwd,'..','Data','HFOEEG');
% analysisfilepath = fullfile(pwd,'..','Data','PatientData',strcat(umid,'.mat'));
% bkgdatafilepath = fullfile(pwd,'..','Data','BKGEEG',strcat(umid,'.mat'));
%%%

%%
%Load Background Analysis Data
load(bkgdatafilepath,'datainfo');
bkgdatainfo = datainfo;
edges = bkgdatainfo.rmsdata.edges;
allBKGpdfCB = bkgdatainfo.rmsdata.lagkdepdfCB;
allBKGsupCB = bkgdatainfo.rmsdata.lagkdesupCB;
clear bkgdatainfo
                
%Load Real Analysis Data
load(analysisfilepath,'PatientInfo','elecct','datainfo');

%ANALYSIS WINDOW PARAMETERS
MINSAMPREQ = 5;                 %ignore electrodes with less than minimal sample requirements
zeroradius = 1;

%PROCESSING PARAMS
datatype = 'rmsdata';
alpha0 = .05;

%SETTING VARS
CM = struct();
    
CM.(datatype).ccFCM = nan(elecct,elecct);      %elec x elec of unnormalized FCM
CM.(datatype).ccLCM = nan(elecct,elecct);      %Lag connectivity matrix

%baseline normalization info
datainfo.(datatype).usedmaxlagsamps = {};

%bootstrap info
datainfo.(datatype).BSCC = {};
datainfo.(datatype).BSmeans = nan(elecct,elecct);
datainfo.(datatype).BSstds = nan(elecct,elecct);

%bootstrap info  
allRealpdfCB = datainfo.(datatype).lagkdepdfCB;
allRealsupCB = datainfo.(datatype).lagkdesupCB;   

%% Run analysis
for eleci = 1:elecct
    fprintf('%s: Elec_ %i\n', umid, eleci);
    filename = strcat(umid,"-",char(strcat("elec",string(eleci),".mat")));
    load(fullfile(char(datafolderpath),filename),...
        'FS','data','hfoI');
    
    curhfoi = 1:length(hfoI);
    
    datainfo.allHFOi{eleci,1} = hfoI(curhfoi);
    cursampct = length(curhfoi);
    datainfo.hfosamplects(eleci,1) = cursampct;

    %Run CC only if there is enough samples
    if(cursampct >= MINSAMPREQ)
        curusedmaxlagsamps = zeros(elecct,length(curhfoi));
        for elecj = 1:elecct
            curmaxlagmat = datainfo.(datatype).maxcclags{eleci,elecj};
            curmaxampmat = datainfo.(datatype).maxccamps{eleci,elecj};
            if(eleci ~= elecj)
                curRealpdfCB = allRealpdfCB{eleci,elecj};
                curRealsupCB = allRealsupCB{eleci,elecj};        

                curkdepdfCB = allBKGpdfCB{eleci,elecj};
                curkdesupCB = allBKGsupCB{eleci,elecj};

                %Z test of confidence bands for signifance
                curZs = [];
                for i = 1:length(curRealpdfCB)
                    curZs(i,1) = ztest2s(curRealpdfCB(i),curRealsupCB,...
                        curkdepdfCB(i),curkdesupCB,alpha0); 
                end
                curissiglags = curZs >= norminv(1-alpha0);
                %get sig lag bounds

                xdiff = diff([0; curissiglags(:); 0]);
                leftboundidx = find(xdiff == 1);
                rightboundidx = find(xdiff == -1)-1;
                leftboundlags = edges(leftboundidx + 1);
                rightboundlags = edges(rightboundidx + 1);
                
                for sidx = 1:length(curhfoi)
                    curusedmaxlagsamps(elecj,sidx) = ...
                        any((curmaxlagmat(sidx) >= leftboundlags) & ...
                        (curmaxlagmat(sidx) <= rightboundlags)) & ...
                        (abs(curmaxlagmat(sidx)) > zeroradius);
                end 
                
                cursampcts = sum(curusedmaxlagsamps(elecj,:));
                
                if(cursampcts >= MINSAMPREQ)
                    meancclagatpromlag = nansum(curmaxlagmat .* curusedmaxlagsamps(elecj,:)') ./ ...
                        sum(curusedmaxlagsamps(elecj,:));     %actual mean lag times within most prominent lag bin
                    meanccampatpromlag = nansum(curmaxampmat .* curusedmaxlagsamps(elecj,:)')./ ...
                        sum(curusedmaxlagsamps(elecj,:));   %mean CCamp at most prominent lag bin
                    CM.(datatype).ccFCM(eleci,elecj) = meanccampatpromlag;
                    CM.(datatype).ccLCM(eleci,elecj) = meancclagatpromlag;
                end
                    
            end
        end        
        
        datainfo.(datatype).usedmaxlagsamps{eleci} = curusedmaxlagsamps;
        
        %Max CC of bootrapped data (electrode x samples)
        BSCC = singleElecBS_FFT(data.(datatype)(curhfoi),eleci,FS,...
            datainfo.(datatype).usedmaxlagsamps{eleci});
        datainfo.(datatype).BSCC{eleci,1} = BSCC;
        datainfo.(datatype).BSmeans(eleci,:) = nanmean(BSCC,2);
        datainfo.(datatype).BSstds(eleci,:) = nanstd(BSCC,[],2);
    end
end

% Save Results
save(analysisfilepath,'CM','datainfo','-append')

end
function zs = ztest2s(mu1,sup1,mu2,sup2,alpha)

z = norminv(1-alpha);

SE1 = sup1/z;
SE2 = sup2/z;

zs = (mu1-mu2)/sqrt(SE1^2 + SE2^2);

end