function [Ycap varexpl]=call_HONMF(X,opts)
%% call interface for HONMF for the purpose of TDALAB
%   
%
% Input:
%   X                 Array of non negative data
%   options
%     .NumOfComp      1 x length(FACT) vector containing number of factors
%                     for each modality, i.e. NumOfComp(i) is number of factors in
%                     FACT{i}.
%     .costfcn        Cost function to optimize
%                       'ls': Least squares (default)
%                       'kl': Kullback Leibler
%     .constFACT      1 x length(FACT) vector indicating which factors to update,
%                     i.e. constFACT(i)=0: updates FACT{i}  (default)
%                          constFACT(i)=1: FACT{i} is not updated.
%     .constCore      constCore=0: updates Core  (default)
%                     constCore=1: Core is not updated.
%     .lambda         Sparsity weight on Core and Factors,
%                     (sparsity function used is the 1-norm)
%                       lambda(1):      Core
%                       lambda(i+1):    FACT{i}
%     .maxiter        Maximum number of iterations (default 100)
%     .accel          Wild driver accelleration parameter (default 1.3)
%     .displaylevel   Level of display: [off | final | iter]
%% HONMF: Higher Order NON-NEGATIVE MATRIX FACTORIZATION
%
% Authors:
%   Morten Mørup
%   Technical University of Denmark,
%   Institute for Matematical Modelling
%
% Reference:
%   M. Mørup et al. Algorithms for Sparse Higher Order Non-negative Tensor Factorization. Technical University of Denmark, 2006.
%   M. M?rup, L.K. Hansen, and S.M. Arnfred, ¡°Algorithms for sparse
%   nonnegative tucker decompositions,¡± Neural Computation, vol. 20, no.
%   8, pp. 2112-2131, 2008.
% Usage:
%       [Core, FACT] = HONMF(X, d, varargin)
%
% Example:
%       opts.lambda=[1 0 0 1];
%       opts.maxiter=1000;
%       X=rand(10,20,30);
%       [Core, FACT] = HONMF(X, [3 3 3], opts)
%
% Input:
%   X                 Array of non negative data
%   d                 1 x length(FACT) vector containing number of factors
%                     for each modality, i.e. d(i) is number of factors in
%                     FACT{i}.
%   options
%     .costfcn        Cost function to optimize
%                       'ls': Least squares (default)
%                       'kl': Kullback Leibler
%     .FACT           Initial FACT
%     .Core           Initial Core
%     .W              Indicator tensor of same size as X, 0 where data is missing
%                     1 where present
%     .constFact      1 x length(FACT) vector indicating which factors to update,
%                     i.e. constFACT(i)=0: updates FACT{i}  (default)
%                          constFACT(i)=1: FACT{i} is not updated.
%     .constCore      constCore=0: updates Core  (default)
%                     constCore=1: Core is not updated.
%     .lambda         Sparsity weight on Core and Factors,
%                     (sparsity function used is the 1-norm)
%                       lambda(1):      Core
%                       lambda(i+1):    FACT{i}
%     .maxiter        Maximum number of iterations (default 100)
%     .conv_criteria  Function exits when cost/delta_cost exceeds this
%     .accel          Wild driver accelleration parameter (default 1.3)
%     .displaylevel   Level of display: [off | final | iter]
%
% Output:
%   Core                 d(1) x d(2) x ... x d(N) array
%   FACT                 cell array of factors for each modality, i.e.
%                        FACT{i} is a size(V,i) x d(i) matrix
%   varexpl              Explained variation
%
% Copyright (C) Morten Mørup and Technical University of Denmark,
% September 2006
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%   Edit history:
%   30 Oct 2008:
%           varexplained changed from relative to
%           SST=norm(X-mean(X),'fro')^2 to SST= norm(X,'fro')^2
%           error corrected in initialization
X=double(X);
sz=size(X);
NumOfMode=numel(sz);
defopts=struct('NumOfComp',[],'costfcn','ls','constFACT',zeros(1,NumOfMode),'constCore',0,'lambda',zeros(1,NumOfMode+1),...
    'maxiter',100,'accel',1.3,'displaylevel','off','conv_criteria',1e-6);
if ~exist('opts','var')
    opts=struct;
end
opts=scanparam(defopts,opts);

[Core, FACT, varexpl] = HONMF(X, opts.NumOfComp, opts);

Ycap=ttensor(tensor(Core),FACT);

% .costfcn        Cost function to optimize
%                       'ls': Least squares (default)
%                       'kl': Kullback Leibler
%     .FACT           Initial FACT
%     .Core           Initial Core
%     .W              Indicator tensor of same size as X, 0 where data is missing
%                     1 where present
%     .constFact      1 x length(FACT) vector indicating which factors to update,
%                     i.e. constFACT(i)=0: updates FACT{i}  (default)
%                          constFACT(i)=1: FACT{i} is not updated.
%     .constCore      constCore=0: updates Core  (default)
%                     constCore=1: Core is not updated.
%     .lambda         Sparsity weight on Core and Factors,
%                     (sparsity function used is the 1-norm)
%                       lambda(1):      Core
%                       lambda(i+1):    FACT{i}
%     .maxiter        Maximum number of iterations (default 100)
%     .conv_criteria  Function exits when cost/delta_cost exceeds this
%     .accel          Wild driver accelleration parameter (default 1.3)
%     .displaylevel   Level of display: [off | final | iter]