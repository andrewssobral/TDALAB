function [ Ycap ] = FastNTFAPG( Y,opts )
%% Low-rank approximation based Non-Negative CP Decomposition
%   Usage: Ycap=lraNCP(Y,opts);
%   Input: Y can be ktensor/ttensor/tensor/double
%  Output: Ycap is a nonnegative ktensor
%
%  Parameters:
%   opts.NumOfComp:  rank of Ycap
%   opts.lra_model:  [CP]|Tucker|none -- low-rank approximation model
%   opts.lra_rank:   rank of low-rank approximation of Y
%   opts.lra_iter:   max. of low-rank approximation iterations
%   opts.maxiter:    max. of global iterations    [default:50]
%   opts.inner_maxiter: max. of inner iterations. [default:20]
%   opts.tol:        stop threshold. [default:1e-6]
%   opts.alpha,beta  1-alpha*exp(-beta*it) is used as the step
%
%
%
% Cite this work currently:
%  Guoxu Zhou; Cichocki, A.; Shengli Xie; , "Fast Nonnegative Matrix/Tensor Factorization Based on Low-Rank Approximation,"
%  IEEE Transactions on Signal Processing, vol.60, no.6, pp.2928-2940, June 2012
%  doi: 10.1109/TSP.2012.2190410
%  URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166354&isnumber=6198804
%
% By Guoxu Zhou
% E-Mail: zhouguoxu@gmail.com
