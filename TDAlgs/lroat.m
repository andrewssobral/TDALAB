function [T, fit] = lroat(A, opts)
%function [T, fit] = lroat(A,r, opts)
%
% LROAT: Low Rank Orthogonal Approximation of Tensors. This
% function computes a rank-r Kruskal tensor T (with orthogonal mode
% vectors) that best approximates the input tensor A.
%
% This code requires the tensor toolbox developed by Brett Bader
% and Tamara Kolda at the Sandia National Lab. The toolbox can be
% downloaded from http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox/
%
% Input:
%   [A] Tensor. Datatype: tensor, ktensor or ttensor.
%   [r] Rank.
%   [opts] (Optional)
%    opts.Tol: Tolerance for test of convergence. {1e-4}
%    opts.MaxIter: Maximum number of iterations. {100}
%    opts.Init: Initial guess. [{'nvecs'}|'random'|cell array]
%
% Output:
%   [T] The approximated tensor. Datatype: ktensor.
%       The coefficients are nonnegative and in decreasing order.
%   [normT] The norm of T for all iterations.
%   [difU] The differences of U matrices between iterations.
%
% Example:
%   A = tenrand([20 16 10 32]);
%   r = 5;
%   opts = struct('Tol',1e-5, 'MaxIter',1000, 'init','nvecs');
%   [T, normT, difU] = lroat(A, r, opts);
%
% See also PARAFAC_ALS, TUCKER_ALS.

% Copyright 2009, Jie Chen.
% This program is free software; you can redistribute and/or
% modify it for NON-COMMERCIAL purposes. This program is
% distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY, including that of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE.
% $Date: 2009/01/20 19:57:32$


%------------------------------------------------------------
%  preparation
%------------------------------------------------------------
% parameters
%   [opts] (Optional)
%    opts.Tol: Tolerance for test of convergence. {1e-4}
%    opts.MaxIter: Maximum number of iterations. {100}
%    opts.Init: Initial guess. [{'nvecs'}|'random'|cell array]

defoptions = struct('NumOfComp',0,'Tol',1e-4,'MaxIter',100,'init','nvecs');
if ~exist('opts','var')
    opts = struct;
end
[r opts.Tol opts.MaxIter opts.init]= scanparam(defoptions,opts);

d = size(A);                            % dimensions
N = length(d);                          % order
if r > min(d)                           % rank
  error('r larger than min{dimensions}');
end

% force the datatype of A to be `tensor'
A = full(A);
normA = norm(A);

% allocate memory
sigma = zeros(r,1);                     % sigma
U = cell(1,N);                          % Un matrices
Uold = cell(1,N);
V = cell(1,N);                          % Vn matrices
for n = 1:N
  V{n} = zeros(d(n),r);
end
Ucp = cell(1,r);                        % a copy of the Un matrices
for i = 1:r                             % in a different data structure
  Ucp{i} = cell(1,N);
  for n = 1:N
    Ucp{i}{n} = zeros(d(n),1);
  end
end
normT = zeros(opts.MaxIter,1);          % normT
difU = zeros(opts.MaxIter,N);           % difU

% initial guess
if strcmp(opts.init,'load')
    load('TDinit.mat','Ainit'); % cell-type 'init' is loaded from CP_init.mat
    U = Ainit;
    clear Ainit;
elseif strcmp(opts.init, 'random')
  for n = 1:N
    U{n} = qr(randn(d(n),r), 0);
  end
elseif strcmp(opts.init, 'nvecs')
  % Can set U{n} = nvecs(A, n, r);
  % But nvecs() might not return vectors in correct order
  for n = 1:N
    An = double(tenmat(A,n));
    if r == d(n)
      [U{n}, D] = eig(An*An');
    else
      [U{n}, D] = eigs(An*An', r, 'LM', struct('disp',0));
    end
    [D, idx] = sort(diag(D), 'descend');
    U{n} = U{n}(:,idx);
  end
else
  error('The initialization method is not supported.');
end

% initialize Ucp
for n = 1:N
  for i = 1:r
    Ucp{i}{n} = U{n}(:,i);
  end
end


%------------------------------------------------------------
%  main iteration
%------------------------------------------------------------

for j = 1:opts.MaxIter
  iter = j;

  % run through all n
  for n = 1:N

    % compute Vn and sigma
    % ttv is indeed slow...
    for i = 1:r
      V{n}(:,i) = double(tenmat(ttv(A, Ucp{i}, -n), 1));
      sigma(i) = U{n}(:,i)'*V{n}(:,i);
    end

    % polar decomposition and update Un
    % there are more sophisticated ways to do polar-decomp
    [Up, Dp, Vp] = svd(V{n}*diag(sigma), 0);
    Un = Up*Vp';
    difU(j,n) = norm(U{n}-Un);
    U{n} = Un;

    % update Ucp
    for i = 1:r
      Ucp{i}{n} = U{n}(:,i);
    end

  end

  % normT
  for i = 1:r
    sigma(i) = ttv(A, Ucp{i});
  end
  normT(j) = norm(sigma);

  % convergence test
  if j > 1 && norm(U{1}-Uold{1}) < opts.Tol
    break;
  elseif j == opts.MaxIter
    fprintf(1,'Maximum number of iterations (%d) reached.\n',opts.MaxIter);
    fprintf(1,'Convergence has not been observed.\n');
  end
  
  % update Uold
  for n = 1:N
    Uold{n} = U{n};
  end

end


%------------------------------------------------------------
%  output
%------------------------------------------------------------
% revert sign
for i = 1:r
  if (sigma(i) < 0)
    U{1}(:,i) = -U{1}(:,i);
    sigma(i) = -sigma(i);
  end
end

% sort
[sigma, idx] = sort(sigma,'descend');
for n = 1:N
  U{n} = U{n}(:,idx);
end

% output
T = ktensor(sigma,U);
if nargout==2
    fit=1-norm(tensor(T)-tensor(A))/norm(A);
end
