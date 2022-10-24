%%%
%Compute Critical Resection Percentage from ranked feature values (ie.
%centrality, hforate, etc.)
%
%Jack Lin
%MATLAB R2021a
%10/24/2022 
%%%
%%
umid = "UMHS-0022";
analysisfilepath = fullfile(pwd,'Data','PatientData',strcat(umid,'.mat'));
%%

load(analysisfilepath,'PatientInfo','valTable','elecct'); 

elabels = PatientInfo.Electrodes.elabels;
valtypes = unique(valTable.val_type);

%Participation by percentile thresholding
groupPrtile = [1,.9:-.1:.5,0];
sumTabVarNames = {'val_type','percentile','SOZct','RVct','BOTHct',...
    'RESTct','SOZper','RVper','BOTHper','RESTper'};
crepTable = cell2table(cell(0,length(sumTabVarNames)),'VariableNames',sumTabVarNames);
crepTable = convertvars(crepTable,{'val_type'},'categorical');

for cidx = 1:length(valtypes)
    curvaltype = valtypes(cidx);
    curvals = valTable.val(valTable.val_type == curvaltype);
    
    assert(elecct == length(curvals));

    [~,sortidx] = sort(curvals,'descend','MissingPlacement',...
        'last');       %Nan values will be delegated LEAST important position
    
    sortedelabels = elabels(sortidx);
    sSOZlogic = (ismember(sortedelabels,{'SOZ','BOTH'}));
    sRVlogic = (ismember(sortedelabels,{'RV','BOTH'}));
    sBOTHlogic = (sortedelabels == 'BOTH');
    sRESTlogic = (sortedelabels == 'REST');

    SOZcumct = cumsum(sSOZlogic);
    RVcumct = cumsum(sRVlogic);
    BOTHcumct = cumsum(sBOTHlogic);
    RESTcumct = cumsum(sRESTlogic);
    
    prtileK = max([1],[ceil(length(elabels) .* (1-groupPrtile))])';
    gpn = length(groupPrtile);

    SOZatK = SOZcumct(prtileK);
    RVatK = RVcumct(prtileK);
    BOTHatK = BOTHcumct(prtileK);
    RESTatK = RESTcumct(prtileK);

    
    crepTable.val_type(end+1:end+gpn) = curvaltype;
    crepTable.percentile(end-gpn+1:end) = groupPrtile;
    crepTable.SOZct(end-gpn+1:end) = SOZatK;
    crepTable.SOZper(end-gpn+1:end) = SOZatK./prtileK;
    crepTable.RVct(end-gpn+1:end) = RVatK;
    crepTable.RVper(end-gpn+1:end) = RVatK./prtileK;
    crepTable.BOTHct(end-gpn+1:end) = BOTHatK;
    crepTable.BOTHper(end-gpn+1:end) = BOTHatK./prtileK;
    crepTable.RESTct(end-gpn+1:end) = RESTatK;
    crepTable.RESTper(end-gpn+1:end) = RESTatK./prtileK;
    
end

