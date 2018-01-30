function [mappedX mapping]=PMFcompute_mapping(X,opts)
%% Dimensionality reduction methods
% Usage: [mappedX mapping]=PMFcompute_mapping(X,opts);
% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology
SUPPORTED_DR={'tSNE','ISOMAP','LLE','PCA'};
defopts=struct('method','tSNE','no_dims',3,'NNeighbors',12,'init_dims',30);
if ~exist('opts','var')
    opts=struct();
end
[drmethod no_dims k init_dims]=scanparam(defopts,opts);
%% define default method
if ~find(strcmpi(drmethod,SUPPORTED_DR)==1)
    drmethod='PCA';
end

X=double(X);
no_dims=min(no_dims,size(X,2));

switch lower(drmethod)
    case 'tsne'
        init_dims=min(init_dims,size(X,2));
        mappedX=tsne(X,[],no_dims,init_dims);
    case 'lle'
        mappedX=lle(X,no_dims,k);
    case 'isomap'
        mappedX=isomap(X,no_dims,k);
    case 'pca'        
        mappedX=pca(X,no_dims);
end
mapping=mappedX\X;        

end