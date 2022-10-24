%%%
%Get top X percentile electrodes based on each feature type (ie.
%centrality, hforate, etc.) for each electrode type (ie. resected volume
%RV, seizure onset zone SOZ, etc.)
%
%Jack Lin
%MATLAB R2021a
%10/24/2022 
%%%
%%
umid = "UMHS-0028";
analysisfilepath = fullfile(pwd,'Data','PatientData',strcat(umid,'.mat'));
%%

load(analysisfilepath,'PatientInfo','valTable','valInfoTable'); 

einfo = PatientInfo.Electrodes;
valtypes = unique(valTable.val_type);
valtabVarNames = valTable.Properties.VariableNames;

%New Centrality by Percentile Table
groupPrtile = [0,.5,.9,1];   %Percentile threshold for each group (SOZ/RV/REST)
prtileTable = cell2table(cell(0,length(valtabVarNames)+1),'VariableNames',...
    [valtabVarNames,{'percentile'}]);
prtileTable = convertvars(prtileTable,{'elec_type','val_type'},'categorical');

for cidx = 1:length(valtypes)
    curvaltype = valtypes(cidx);
    cureremains = valInfoTable.elec_remain{valInfoTable.val_type == curvaltype};

    %Translate percentile # to acutal numbers within each group (SOZ/RV/REST)
    SOZremains = cureremains(einfo.SOZ(cureremains));
    RVremains = cureremains(einfo.RV(cureremains));
    RESTremains = cureremains(einfo.REST(cureremains));
    notRVremains = cureremains(~einfo.RV(cureremains));
    BOTHremains = cureremains(einfo.SOZ(cureremains) & einfo.RV(cureremains));

    SOZonlyremains = cureremains(einfo.SOZ(cureremains) & ~einfo.RV(cureremains));
    RVonlyremains = cureremains(~einfo.SOZ(cureremains) & einfo.RV(cureremains));

    SOZk = [ceil(length(SOZremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(SOZremains)); SOZk(end) = 0; end
    RVk = [ceil(length(RVremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(RVremains)); RVk(end) = 0; end
    RESTk = [ceil(length(RESTremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(RESTremains)); RESTk(end) = 0; end
    notRVk = [ceil(length(notRVremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(notRVremains)); notRVk(end) = 0; end
    BOTHk = [ceil(length(BOTHremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(BOTHremains)); BOTHk(end) = 0; end

    SOZonlyk = [ceil(length(SOZonlyremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(SOZonlyremains)); SOZonlyk(end) = 0; end
    RVonlyk = [ceil(length(RVonlyremains) .* (1-groupPrtile(1:end-1)))';1];
    if(isempty(RVonlyremains)); RVonlyk(end) = 0; end

    iscurctypeanderemain = (valTable.val_type == curvaltype &...
        ismember(valTable.elec,cureremains));

    curSOZtable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'SOZ','BOTH'}),:);
    curSOZtable.elec_type(:) = categorical("SOZ");
    curRVtable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'RV','BOTH'}),:);
    curRVtable.elec_type(:) = categorical("RV");
    curRESTtable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'REST'}),:);
    curNotRVtable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'REST','SOZ'}),:);
    curNotRVtable.elec_type(:) = categorical("notRV");
    curBOTHtable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'BOTH'}),:);

    curSOZonlytable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'SOZ'}),:);
    curSOZonlytable.elec_type(:) = categorical("SOZonly");
    curRVonlytable = valTable(iscurctypeanderemain & ...
        ismember(valTable.elec_type,{'RV'}),:);
    curRVonlytable.elec_type(:) = categorical("RVonly");

    for gpidx = 1:length(groupPrtile)
        [~,subSOZi] = mink(curSOZtable.rank(:),SOZk(gpidx));
        prtileTable(end+1:end+SOZk(gpidx),valtabVarNames) = curSOZtable(subSOZi,:);
        prtileTable.percentile(end-SOZk(gpidx)+1:end) = groupPrtile(gpidx);

        [~,subRVi] = mink(curRVtable.rank(:),RVk(gpidx));
        prtileTable(end+1:end+RVk(gpidx),valtabVarNames) = curRVtable(subRVi,:);
        prtileTable.percentile(end-RVk(gpidx)+1:end) = groupPrtile(gpidx);

        [~,subRESTi] = mink(curRESTtable.rank(:),RESTk(gpidx));
        prtileTable(end+1:end+RESTk(gpidx),valtabVarNames) = curRESTtable(subRESTi,:);
        prtileTable.percentile(end-RESTk(gpidx)+1:end) = groupPrtile(gpidx);

        [~,subnotRVi] = mink(curNotRVtable.rank(:),notRVk(gpidx));
        prtileTable(end+1:end+notRVk(gpidx),valtabVarNames) = curNotRVtable(subnotRVi,:);
        prtileTable.percentile(end-notRVk(gpidx)+1:end) = groupPrtile(gpidx);

        [~,subBOTHi] = mink(curBOTHtable.rank(:),BOTHk(gpidx));
        prtileTable(end+1:end+BOTHk(gpidx),valtabVarNames) = curBOTHtable(subBOTHi,:);
        prtileTable.percentile(end-BOTHk(gpidx)+1:end) = groupPrtile(gpidx);

        [~,subSOZonlyi] = mink(curSOZonlytable.rank(:),SOZonlyk(gpidx));
        prtileTable(end+1:end+SOZonlyk(gpidx),valtabVarNames) = curSOZonlytable(subSOZonlyi,:);
        prtileTable.percentile(end-SOZonlyk(gpidx)+1:end) = groupPrtile(gpidx);

        [~,subRVonlyi] = mink(curRVonlytable.rank(:),RVonlyk(gpidx));
        prtileTable(end+1:end+RVonlyk(gpidx),valtabVarNames) = curRVonlytable(subRVonlyi,:);
        prtileTable.percentile(end-RVonlyk(gpidx)+1:end) = groupPrtile(gpidx);
    end
end
