%%%
%Compute the mean maximum cross correlation between EEG of electrode i with
%the randomly scrambled EEGs of every other electrode j
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% data:         s x 1 cell (samples) of n x m (n electrodes) x (m datapoints)
% elecRoi:      electrode of interest to perform cross correlation against
% FS:           Hz sampling rate
% usesample:    1 x n (electrode) cell of n x s boolean (n electrodes) x (s
%               samples) of samples to use
%
% OUTPUT
%
% BSCC:         n x 1 (electrode) cell of n x s of mean of maximum cross
%               correlations
%%%

function [BSCC] = singleElecBS_FFT(data,elecRoi,FS,usesample)
rng(1);
Shuffle(1,'seed');

mainbound = [-100 100];         %time for CC
midpointidx = floor(size(data{1,1},2)/2)+1;
mainidxoffset = (mainbound.*FS./1000)+midpointidx;
datalength = diff(mainidxoffset) + 1;

elecct= size(usesample,1);
iterct = 100;

BSCC = nan(elecct,size(data,1));

paddedzerovec = zeros(1,2*datalength -1);
paddedzeromat = zeros(iterct,2*datalength-1);
istartidx = datalength;

for sampi = 1:length(data)
    elecusingcursamp = find(usesample(:,sampi));
    if(~isempty(elecusingcursamp))
        elecROIOGrms = data{sampi}(elecRoi,mainidxoffset(1):mainidxoffset(2));   %RMS data for ElecROI
        elecROIAutoCorr = elecROIOGrms * elecROIOGrms';
        
        curn = length(elecROIOGrms);    
        pelecROIrms = paddedzerovec;
        elecROIBSFFT = paddedzeromat;
        for i = 1:iterct
            pelecROIrms(1,istartidx:end) = elecROIOGrms(randperm(curn));
            elecROIBSFFT(i,:) = fft(pelecROIrms);
        end
        
        elecijAutoFactor = sqrt(elecROIAutoCorr * ...
            sum(data{sampi}(elecusingcursamp,mainidxoffset(1):mainidxoffset(2)).^2,2));
        
        %elecj        
        for elecj = 1:length(elecusingcursamp)
            
            pelecjrms = paddedzerovec;
            pelecjrms(1:datalength) = data{sampi}(elecusingcursamp(elecj),mainidxoffset(1):mainidxoffset(2));
            elecjBSFFT = fft(pelecjrms);
            
            elecijpreCC = elecROIBSFFT .* conj(elecjBSFFT);
%             elecijCC = ifft(elecijpreCC');      %data x sample (flipped)
            elecijCC = paddedzeromat';
            
            for bi = 1:iterct
                elecijCC(:,bi) = ifft(elecijpreCC(bi,:));
            end
            
            negsigns = sign(max(elecijCC)-min(elecijCC));
            elecjBSCC = max(abs((1/elecijAutoFactor(elecj))* elecijCC )).*negsigns;

            BSCC(elecusingcursamp(elecj),sampi) = nanmean(elecjBSCC);
        end
    end
end

end