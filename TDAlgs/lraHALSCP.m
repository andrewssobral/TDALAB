function [ Ycap ] = lraHALSCP( Y,opts )
%% HALS (Non-Negative) CP Decomposition
%   Usage: Ycap=lraHALSCP(Y,opts);
%   Input: Y can be ktensor/ttensor/tensor/double
%  Output: Ycap is a nonnegative ktensor
%
%  Parameters:
%   opts.NumOfComp:       rank of Ycap
%   opts.lra_model:       [CP]|Tucker|none -- low-rank approximation model
%   opts.lra_rank:        rank of low-rank approximation of Y
%   opts.lra_iter:        max. of low-rank approximation iterations
%   opts.maxiter:         max. of global iterations    [default:50]
%   opts.inner_maxiter:   max. of inner iterations. [default:20]
%   opts.tol:             stop threshold. [default:1e-6]
%   opts.constraint:      none|nonnegative|sparse
%   opts.initmode:        [CP]|random|sampling  initialization mode
%   opts.cons_param:      parameters for constraints. Currently only works for
%                       sparse. For sparsity, it is a vector of nonnegative numbers between 0 and 1 with
%                       size of N. 
%
% Cite these papers:
%  [1] A. Cichocki and Phan A-H. Fast local algorithms for large scale Nonnegative Matrix and Tensor Factorizations .
%      IEICE Transaction on Fundamentals, (2009)
%  [2] Guoxu Zhou; Cichocki, A.; Shengli Xie; , "Fast Nonnegative Matrix/Tensor Factorization Based on Low-Rank Approximation," 
%       Signal Processing, IEEE Transactions on , vol.60, no.6, pp.2928-2940, June 2012
%       URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166354&isnumber=6198804
%
% Coded by Guoxu Zhou
% E-Mail: zhouguoxu@gmail.com
