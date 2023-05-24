%%%
%Compute graph object from adjacency matrix
%
%Jack Lin
%MATLAB R2022b
%5/1/23
%
% INPUT
%
% adjmat:       n x n adjacency matrix
% weightmat:    n x n weight matrix
% directed:     boolean if adjacency matrix is directed (asymmetric)
%
% OUTPUT
%
% graphobj:     graph object representation of adjacency matrix
%%%

function graphobj = coord3DGrapher(adjmat,weightmat,directed)

if(~exist('directed','var')); directed = true; end

if(~directed)
    adjmat = triu(adjmat,1);
    weightmat = triu(weightmat,1);
end

adjmat(isnan(adjmat)) = 0;
weightmat(isnan(weightmat)) = 0;

elecct = size(adjmat,1);
elecidxperm = repmat(1:elecct,elecct,1);
fromidx = elecidxperm';
fromidx = fromidx(logical(adjmat));     %from nodes
toidx = elecidxperm(logical(adjmat));   %to nodes

if(directed)
    graphobj = digraph(fromidx,toidx,weightmat(logical(adjmat)));
else
    graphobj = graph(fromidx,toidx,weightmat(logical(adjmat)));
end

if(size(graphobj.Nodes,1) ~= elecct)
    graphobj = addnode(graphobj,elecct-size(graphobj.Nodes,1));
end