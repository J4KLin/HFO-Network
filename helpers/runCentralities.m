%%%
%Compute the functional connectivity matricies (FCN,uLAN,fLAN) and their
%respective centralities (FCN-EIG,uLAN-EIG,fLAN-OUT)
%
%Jack Lin
%MATLAB R2022b
%5/1/2023
%
% INPUT
%
% analyasisfilepath:    file path of patient HFO cross correlation data
% bkgdatafilepath:      file path of patient background cross correlation
%                       data
%
%%%

function runCentralities(analysisfilepath,bkgdatafilepath)
%%%
% umid = "UMHS-0028"; 
% analysisfilepath = fullfile(pwd,'..','Data','PatientData',strcat(umid,'.mat'));
% bkgdatafilepath = fullfile(pwd,'..','Data','BKGEEG',strcat(umid,'.mat'));   
%%%
%%
load(bkgdatafilepath,'datainfo');
bkgmaxcclags = datainfo.rmsdata.maxcclags;
bkgmaxccamps = datainfo.rmsdata.maxccamps;
zeroradius = 1;
clear datainfo
load(analysisfilepath);

%% Get uFCN-EIGenvector centrality
nmat = nan(elecct);     %# of selected sampled for each eleci vs elecj CC
for eleci = 1:elecct
   nmat(eleci,:) = sum(datainfo.rmsdata.usedmaxlagsamps{eleci},2);
end

BKGBSnmat = nan(elecct);
BKGBS_mean = nan(elecct);
BKGBS_std = nan(elecct);

%Get mean and std of maxCCamp of BKG data
for eleci = 1:elecct
    for elecj = 1:elecct
        curn = nmat(eleci,elecj);
        if(eleci ~= elecj && curn > 0)
           %Get BKG maxamp data ignoring zero lag radius samples
           curBKGdata = bkgmaxccamps{eleci,elecj}(abs(bkgmaxcclags{eleci,elecj}) > zeroradius);
           BKGBSnmat(eleci,elecj) = length(curBKGdata);
           
           BKGBS_mean(eleci,elecj) = mean(curBKGdata);
           BKGBS_std(eleci,elecj) = std(curBKGdata);         
        end
    end
end

FCN_sim = datainfo.rmsdata.BSmeans;
stdFCNsim = datainfo.rmsdata.BSstds;

FCN_obs = CM.rmsdata.ccFCM;
FCN_obs(nmat <= 5) = nan;   %remove ei vs ej with sample less than 5
FCNZ_obs = ((FCN_obs - FCN_sim)) ./ stdFCNsim;
FCNZ_obs(isinf(FCNZ_obs)) = nan;
q2575 = quantile(FCNZ_obs,[.25 .75],'all');
FCNZ_obsQR = FCNZ_obs;
FCNZ_obsQR(FCNZ_obs > (q2575(2) + 3*diff(q2575))) = nan;
FCNZ_obsQR(FCNZ_obs < (q2575(1) - 3*diff(q2575))) = nan;

FCN_EIG = get_undirCM_cent(FCNZ_obs,0);

FCNinfo = struct();
FCNinfo.nmat = nmat;
FCNinfo.FCN_sim_mean = FCN_sim;
FCNinfo.FCN_sim_std = stdFCNsim;
FCNinfo.BKGBS_mean = BKGBS_mean;
FCNinfo.BKGBS_std = BKGBS_std;
FCNinfo.FCN_obs = FCN_obs;
FCNinfo.FCNZ_obs = FCNZ_obs;
FCNinfo.FCN_EIG = FCN_EIG;
%% Create LAN Matricies and get fLAN-OUTcloseness and uLAN-EIGenvector centralities

%Create Lag Asymmetry matrix (Only include Lags within 10ms radius)
lagradius = 10;
innerbound = 1;

%Directed Asyms (Z method)
[asymZmat,asymmat,asymbino_std,nmatZ] = ...
    fullPairwiseLagAsymmetryBino(datainfo.rmsdata.maxcclags,lagradius,innerbound);

%Get Centralities
fLAN_OUT = get_fwdasym_cent(asymZmat,0);
uLAN_EIG = get_undirCM_cent(asymZmat,0);

LANinfo = struct();
LANinfo.nmat = nmatZ;
LANinfo.asymbino_std = asymbino_std;
LANinfo.LAN_obs = asymmat;
LANinfo.LANZ_obs = asymZmat;
LANinfo.fLAN_OUT = fLAN_OUT;
LANinfo.uLAN_EIG = uLAN_EIG;

%% Combine Centralities into a table
valtabVarNames = {'elec','elec_type','val_type','val','norm_val','rank','norm_rank'};
valTable = cell2table(cell(0,length(valtabVarNames)),'VariableNames',valtabVarNames);
valTable = convertvars(valTable,{'elec_type','val_type'},'categorical');
    
infotabVarNames = {'val_type','val_algorithm','elec_remainct','elec_remain'};
valInfoTable = cell2table(cell(0,length(infotabVarNames)),'VariableNames',infotabVarNames);
valInfoTable = convertvars(valInfoTable,{'val_type','val_algorithm'},'categorical');


elabels = PatientInfo.Electrodes.elabels;
elecVec = 1:elecct;

%fLAN_OUT
curvals = fLAN_OUT;
valInfoTable = [valInfoTable; table("fLAN_OUT","outcloseness",...
        sum(~isnan(curvals)),{elecVec(~isnan(curvals))},'VariableNames',infotabVarNames)];
valTable.elec(end+1:end+elecct) = elecVec;
valTable.elec_type(end-elecct+1:end) = elabels;
valTable.val_type(end-elecct+1:end) = "fLAN_OUT";
curvals = curvals .* (sum(~isnan(curvals))-1);
valTable.val(end-elecct+1:end) = curvals;
valTable.norm_val(end-elecct+1:end) = curvals ./ max(curvals);
currnk = nanranking(curvals,10);
valTable.rank(end-elecct+1:end) = currnk;
normrank = currnk - 1;
normrank = normrank ./ max(normrank);
normrank = abs(normrank-1);
valTable.norm_rank(end-elecct+1:end) = normrank;
%uLAN_EIG
curvals = uLAN_EIG;
valInfoTable = [valInfoTable; table("uLAN_EIG","eigenvector",...
        sum(~isnan(curvals)),{elecVec(~isnan(curvals))},'VariableNames',infotabVarNames)];
valTable.elec(end+1:end+elecct) = elecVec;
valTable.elec_type(end-elecct+1:end) = elabels;
valTable.val_type(end-elecct+1:end) = "uLAN_EIG";
valTable.val(end-elecct+1:end) = curvals;
valTable.norm_val(end-elecct+1:end) = curvals ./ max(curvals);
currnk = nanranking(curvals,10);
valTable.rank(end-elecct+1:end) = currnk;
normrank = currnk - 1;
normrank = normrank ./ max(normrank);
normrank = abs(normrank-1);
valTable.norm_rank(end-elecct+1:end) = normrank;
%FCN_EIG
curvals = FCN_EIG;
valInfoTable = [valInfoTable; table("FCN_EIG","eigenvector",...
        sum(~isnan(curvals)),{elecVec(~isnan(curvals))},'VariableNames',infotabVarNames)];
valTable.elec(end+1:end+elecct) = elecVec;
valTable.elec_type(end-elecct+1:end) = elabels;
valTable.val_type(end-elecct+1:end) = "FCN_EIG";
valTable.val(end-elecct+1:end) = curvals;
valTable.norm_val(end-elecct+1:end) = curvals ./ max(curvals);
currnk = nanranking(curvals,10);
valTable.rank(end-elecct+1:end) = currnk;
normrank = currnk - 1;
normrank = normrank ./ max(normrank);
normrank = abs(normrank-1);
valTable.norm_rank(end-elecct+1:end) = normrank;
%HFORATE
% curvals = histcounts(PatientInfo.HFO.stdChanI,elecct);
% valInfoTable = [valInfoTable; table("HFORATE","count",...
%         sum(~isnan(curvals)),{elecVec(~isnan(curvals))},'VariableNames',infotabVarNames)];
% valTable.elec(end+1:end+elecct) = elecVec;
% valTable.elec_type(end-elecct+1:end) = elabels;
% valTable.val_type(end-elecct+1:end) = "HFORATE";
% valTable.val(end-elecct+1:end) = curvals;
% valTable.norm_val(end-elecct+1:end) = curvals ./ max(curvals);
% currnk = nanranking(curvals,10);
% valTable.rank(end-elecct+1:end) = currnk;
% normrank = currnk - 1;
% normrank = normrank ./ max(normrank);
% normrank = abs(normrank-1);
% valTable.norm_rank(end-elecct+1:end) = normrank;

%%
save(analysisfilepath,'LANinfo','FCNinfo','valTable','valInfoTable','-append');

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [asymZmat,asymmat,SDmat,nmat] = fullPairwiseLagAsymmetryBino(allLagCellArr,boundradius,innerbound)

elecct = length(allLagCellArr);
asymmat = nan(elecct);
nmat = nan(elecct);
SDmat = nan(elecct);

for eleci = 1:elecct
    for elecj = 1:elecct
        if(eleci ~= elecj)
            elecijLags = allLagCellArr{eleci,elecj};
            elecijLags(abs(elecijLags) > boundradius | abs(elecijLags) <= innerbound) = nan;  %set all 0 lags and lags outside bound to NAN

            elecijposct = nansum(elecijLags > 0);
            elecijnegct = nansum(elecijLags < 0);
            elecijN = elecijposct + elecijnegct;

            nmat(eleci,elecj) = elecijN;
            asymmat(eleci,elecj) = (elecijposct - elecijnegct)';
            SDmat(eleci,elecj) = 2 .* sqrt(elecijN .* (elecijposct./elecijN) .* (elecijnegct./elecijN));
        end
    end
end
asymZmat = asymmat ./ SDmat;
end

function [asymZmat,asymmat,SDmat,nmat] = fullPairwiseLagAsymmetryBino0(allLagCellArr,boundradius,innerbound)
%PARAMS
% boundradius = 10;     %Only count lags within some radius to 0 (ms)
%
elecct = length(allLagCellArr);

asymmat = nan(elecct);

nmat = nan(elecct);

SDmat = nan(elecct);

for eleci = 1:elecct
    eleciLags = allLagCellArr{eleci};
    eleciLags(abs(eleciLags) > boundradius | abs(eleciLags) <= innerbound) = nan;  %set all 0 lags and lags outside bound to NAN
    
    eleciposct = nansum(eleciLags > 0,2);
    elecinegct = nansum(eleciLags < 0,2);
    eleciN = eleciposct + elecinegct;
    
    nmat(eleci,:) = eleciN;
    
    asymmat(eleci,:) = (eleciposct - elecinegct)';
    
    SDmat(eleci,:) = 2 .* sqrt(eleciN .* (eleciposct./eleciN) .* (elecinegct./eleciN));

end

asymZmat = asymmat ./ SDmat;
end

function processedmat = iqr_outier_removal(asymmat,doinv)

processedmat = asymmat;
if(exist('doinv','var'))
    if(doinv)
        processedmat = 1./(abs(processedmat));
        processedmat(isinf(processedmat)) = nan;
    end
end

%remove 3x interquantile range above
q2575 = quantile(processedmat,[.25 .75],'all');
processedmat(processedmat > (q2575(2) + 3*diff(q2575))) = nan;

end

function undirectedgraph = forcebidirection2(adjmat,thresh)

stackmat = adjmat;
stackmat(:,:,2) = adjmat';
undirectedgraph = (nanmax(abs(stackmat),[],3));
undirectedgraph(isnan(undirectedgraph)) = 0;

if(~exist("thresh",'var')); thresh = 0; end

undirectedgraph(undirectedgraph < thresh) = 0;
undirectedgraph(undirectedgraph == 0) = nan;
end

function [centvals,undirasymmat,elec_remain] = get_undirCM_cent(asymZmat,threshold,eidxonly)
centralitytype = 'eigenvector';
weighttype = 'Importance';

undirasymmat = forcebidirection2(asymZmat,threshold);
undirasymmat = iqr_outier_removal(undirasymmat);

%Ignore all channels that are not specificed in eidxonly
if(exist('eidxonly','var'))
    nanmat = nan(size(asymZmat));
    nanmat(eidxonly,eidxonly) = undirasymmat(eidxonly,eidxonly);
    undirasymmat = nanmat;
end

nanmat = ~isnan(undirasymmat);
elec_remain = find((sum(nanmat,1)' + sum(nanmat,2)) ~= 0);

subasymmat = undirasymmat(elec_remain,elec_remain);  %subset of asymmat with elec_remain

graphobj = coord3DGrapher(subasymmat,subasymmat,false);
subcentvals = centrality(graphobj,centralitytype,weighttype,graphobj.Edges.Weight);

centvals = nan(size(asymZmat,1),1);
centvals(elec_remain) = subcentvals;
end

function [centvals,fwdasymmat,elec_remain] = get_fwdasym_cent(asymZmat,threshold,eidxonly)
centralitytype = 'outcloseness';
weighttype = 'Cost';

fwdasymmat = asymZmat;
fwdasymmat(fwdasymmat >= -threshold) = nan;
fwdasymmat = iqr_outier_removal(fwdasymmat,true);   %Remove outlier & INVERT values

%Ignore all channels that are not specificed in eidxonly
if(exist('eidxonly','var'))
    nanmat = nan(size(asymZmat));
    nanmat(eidxonly,eidxonly) = fwdasymmat(eidxonly,eidxonly);
    fwdasymmat = nanmat;
end

nanmat = ~isnan(fwdasymmat);
elec_remain = find((sum(nanmat,1)' + sum(nanmat,2)) ~= 0);

subasymmat = fwdasymmat(elec_remain,elec_remain);  %subset of asymmat with elec_remain

graphobj = coord3DGrapher(subasymmat,subasymmat);
subcentvals = centrality(graphobj,centralitytype,weighttype,graphobj.Edges.Weight);

centvals = nan(size(asymZmat,1),1);
centvals(elec_remain) = subcentvals;
end