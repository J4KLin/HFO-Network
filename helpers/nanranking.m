%%%
%Ranking vector of numbers from highest (1) to lowest (n)
%IGNORES all nan elements in the vector
% 
%ex:    nanranking([1 5 5 10]) = [4 2 2 1]
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% data:         n x 1 double
%
% optional:
% bins:         scalar double of bins to sort the data into
%
% varargin:
% sort:   sorting directing (default 'descend'); set to 'ascend'
%
% OUTPUT
%
% rnk:          n x 1 rank from 1 to n
% bin_rnk:      bin x 1 rank from 1 to bin
%%%
function [rnk,bin_rnk] = nanranking(data, varargin)

defSorting = 'descend'; %default ranking highest as 1

p = inputParser;
isNumVec = @(x) isnumeric(x) && isvector(x);
binConstraint = @(x) isscalar(x) && isnumeric(x);
addRequired(p,'data',isNumVec);
addOptional(p,'bins',length(data),binConstraint);
addParameter(p,'sort',defSorting);
parse(p,data,varargin{:});

%Code
nonnanidx = find(~isnan(p.Results.data));
nndata = p.Results.data(nonnanidx);

[~,nnrnk] = ismember(nndata,sort(nndata,p.Results.sort));
rnk = nan(size(p.Results.data));
rnk(nonnanidx) = nnrnk;

nnbin_rnk = discretize(nnrnk,linspace(1,max(nnrnk),p.Results.bins+1));
bin_rnk = nan(size(p.Results.data));
bin_rnk(nonnanidx) = nnbin_rnk;
end