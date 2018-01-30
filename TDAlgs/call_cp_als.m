function [Y]=call_cp_als(X,opts)
% call CP_ALS in tensor toolbox 2.5 in TDALAB
% Outputs: 
%           Y: ktensor
%           A: loading factors
%
%CP_ALS Compute a CP decomposition of any type of tensor.
%
%   P = CP_ALS(X,R) computes an estimate of the best rank-R
%   CP model of a tensor X using an alternating least-squares
%   algorithm.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor. The result P is a ktensor.
%
%   P = CP_ALS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations {1}
%
%   [P,U0] = CP_ALS(...) also returns the initial guess.
%
%   [P,U0,out] = CP_ALS(...) also returns additional output that contains
%   the input parameters.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   P = cp_als(X,2);
%   P = cp_als(X,2,'dimorder',[3 2 1]);
%   P = cp_als(X,2,'dimorder',[3 2 1],'init','nvecs');
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   [P,U0,out] = cp_als(X,2,'dimorder',[3 2 1],'init',U0);
%   P = cp_als(X,2,out.params); %<-- Same params as previous run
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
% $Id: cp_als.m,v 1.7 2010/03/19 23:46:32 tgkolda Exp $



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
R=R(1);
if strcmpi(init,'load')
    load('TDinit.mat','Ainit');
    init=Ainit;
    clear Ainit;
    if isempty(init)
        fprintf('[TDALAB] TD_init.mat is empty. Please check the setting of ''init''.\n');
        return;
    end
end

[Y]=cp_als(X,R,'init',init,'maxiters',maxiters,'tol',tol,'dimorder',dimorder,...
    'printitn',printitn);

Y=ktensor(Y);
% normX=norm(X);
% normresidual=sqrt( normX^2 + norm(Y)^2 - 2 * innerprod(X,Y) );
% out = 1 - (normresidual / normX);