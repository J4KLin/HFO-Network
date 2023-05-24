%%%
%Compute kernel density with confidence band of the distribution of data
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% data:         s x 1 double of lag times
% alpha0:       alpha of confidence
% edges:        n x 1 double of edges to compute the density
%
% OUTPUT
%
% datakdepdf_db:n x 1 debiased kernel density
% sup:          n x 1 confidence band supports
% datafit:      n x 1 kernel density
% edges:        n x 1 double of edges for the density
%%%
function [datakdepdf_db,sup,datafit,edges] = kde_wBTCBunbiased(data,alpha0,edges)
%%

data = data(~isnan(data));

datalen = length(data);
n_BT = 1000;
%%
der2 = @(x,y) diff(diff(y)./diff(x))./diff(x(2:length(x)));
debias = @(x,x2,bwidth) x(2:end-1) - (bwidth^2 * .5 .* x2); %debias KDE pdf

datafit = fitdist(data(:),'kernel','Kernel','normal');
datakdepdf = pdf(datafit,edges);
useband = datafit.BandWidth;

datakdepdf_der2 = der2(edges,datakdepdf);   %2nd order der. estimate of KDE
datakdepdf_db = debias(datakdepdf,datakdepdf_der2,useband);  %debiased KDE pdf
        
%%
BTdata = nan(n_BT,datalen);
BTkdepdf = nan(n_BT,length(edges));
BTkdepdf_db = nan(n_BT,length(edges)-2);     %pdf of KDE of Bootstrapped data

parfor bidx = 1:n_BT
   curdataBT = datasample(data,datalen);        %create a sample bootstrap
   BTdata(bidx,:) = curdataBT;
   curBTfit = fitdist(curdataBT(:),'kernel','Kernel','normal','Width',useband); %get KDE of bootstraped data
   curBTkdepdf = pdf(curBTfit,edges);
   BTkdepdf(bidx,:) = curBTkdepdf;
   
   curBTkdepdf_der2 = der2(edges,curBTkdepdf);
   BTkdepdf_db(bidx,:) = debias(curBTkdepdf,curBTkdepdf_der2,useband);
end
kdepdfDIFF_db = abs(BTkdepdf_db - datakdepdf_db);        %difference in debiased-pdfs of BT and data
maxkdepdfDIFF_db = max(kdepdfDIFF_db,[],2);

sup = prctile(maxkdepdfDIFF_db,(1-alpha0)*100);
