function [w h]=PMFnnmf(v,opts)
% call nnmf provided by MATLAB
if ~exist('nnmf.m','file')
    error('Does not find nnmf in your MATLAB.');
end

defopts=struct('NumOfComp',[],'algorithm','als','Display','off','MaxIter',100,'TolFun',1e-4,'TolX',1e-4,...
    'UseParallel','never','UseSubstreams','never','replicates',1);
if ~exist('opts','var')
    opts=struct;
end
[k,alg,Display,MaxIter,TolFun,TolX,UseParallel,UseSubstreams,replicates]=scanparam(defopts,opts);
opts=statset('Display',Display,'MaxIter',MaxIter,'TolFun',TolFun,'TolX',TolX,'UseParallel',UseParallel,...
    'UseSubstreams',UseSubstreams);
[w h]=nnmf(v,k,'algorithm',alg,'replicates',replicates,'options',opts);