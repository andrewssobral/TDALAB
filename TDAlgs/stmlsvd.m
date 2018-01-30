function [Y fit] = stmlsvd(T,opts)
%% [Y fit]=stmlsvd(T,opts)
%    T is a tensor
%    opts.NumOfComp: multilinear rank of components (mlr)
%        .Alpha: order of recover (p)
%    Y is a ttensor
%    fit is optional and saves the fittingness
%% function [A] = stmlsvd(T,mlr,p)
% STMLSVD  Sequentially truncated multilinear singular value decomposition.

    % A = stmlsvd(T,mlr,p) computes the ttensor A, as rank-mlr ST-MLSVD 
    % approximation corresponding to the order p of the tensor T.

    % A = stmlsvd(T,mlr) computes the ttensor A, as rank-mlr ST-MLSVD
    % approximation of T, corresponding to the order of the strongest
    % reduction of the size of the partially truncated core tensor.

    % NOTE: only supports dense unstructured tensors efficiently. Structured 
    % tensors may loose their structure after the first projection step.

%% Adjust for TDALAB
defopts=struct('NumOfComp',0,'p',[]);
if ~exist('opts','var')
    opts = struct;
end
[mlr,p]=scanparam(defopts,opts);
if isempty(p)
    p=1:numel(size(p));
end
    
    
% Sanity check.
mlr = min( size(T), mlr );
    
% Determine the order of the strongest reduction.
if nargin < 3
    [tempuseless,p] = sort( size(T) ./ mlr, 'descend' );
end

% Factor matrices and core tensor.
d = ndims(T);
U = cell(1,d);
S = T;

% ST-MLSVD algorithm.
for k = 1 : d
    m = p(k);
    
    % Compute dominant subspace.
    U{m} = nvecs(S,m,mlr(m));
    
    % Project onto dominant subspace.
    S = ttm(S, U{m}, m, 't');
end

Y = ttensor(S,U);
%% compute the fitting error
if nargout==2
    normT=norm(T);
    fit=1-norm(tensor(Y)-tensor(T))/normT;
end
end
