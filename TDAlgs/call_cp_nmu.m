function [Y out]=call_cp_nmu(X,opts)
% function [Y out]=call_cp_nmu(X,opts)
%call CP_ALS in tensor toolbox 2.4 in TDALAB
% Outputs: 
%           Y: ktensor
%           A: loading factors
%
%CP_NMU Compute nonnegative CP with mutiplicative updates.
%
%   P = CP_NMU(X,R) computes an estimate of the best rank-R PARAFAC
%   model of a tensor X with nonnegative constraints on the factors.
%   This version uses the Lee & Seung multiplicative updates from
%   their NMF algorithm.  The input X can be a tensor, sptensor,
%   ktensor, or ttensor. The result P is a ktensor.
%
%   P = CP_NMU(X,R,OPTS) specify options:
%   OPTS.tol: Tolerance on difference in fit {1.0e-4}
%   OPTS.maxiters: Maximum number of iterations {50}
%   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
%   OPTS.init: Initial guess [{'random'}|'nvecs'|cell array]
%   OPTS.printitn: Print fit every n iterations {1}
%
%   [P,U0] = CP_NMU(...) also returns the initial guess.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   P = cp_nmu(X,2);
%   P = cp_nmu(X,2,struct('dimorder',[3 2 1]));
%   P = cp_nmu(X,2,struct('dimorder',[3 2 1],'init','nvecs'));
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   P = cp_nmu(X,2,struct('dimorder',[3 2 1],'init',{U0}));
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda.
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: cp_nmu.m,v 1.5 2010/03/19 23:46:32 tgkolda Exp $

N=ndims(X);
if ~exist('opts','var')
    opts = struct;
end
defoptions = struct('NumOfComp',[],'tol',1e-12,'maxiters',500,...
    'init','random','printitn',10,'dimorder',1:N);
if ~exist('opts','var')
    opts = struct;
end
[R,tol,maxiters,init,printitn,dimorder] = scanparam(defoptions,opts);
if strcmpi(init,'load')
    load('TDinit.mat','Ainit');
    init=Ainit;
    clear Ainit;
    if isempty(init)
        fprintf('[TDALAB] CP_init.mat is empty. Please check the ''init'' parameter.\n');
        return;
    end
end

opts=struct('tol',tol,'maxiters',maxiters,'init',init,'printitn',...
    printitn,'dimorder',dimorder);

[Y]=cp_nmu(X,R,opts);