function [Ycap]=call_bcdLM(X,opts)
%BCDLM_ALSLS Block Component Decomposition (BCD) in rank-(L,M,.) terms.
%   [A,B,C]=bcdLM_alsls(X,R,L,M) computes the BCD-(L,M,.) of a third-order tensor.
%   The alternating least squares (ALS) algorithm is used, possibly coupled
%   with different choices of line search to speed up convergence.
%
% DESCRIPTION:
% This function computes the Block Component Decomposition (BCD) of the tensor X
% of size (I,J,K) in a sum of R terms of rank-(L,M,.): 
% X = (C1 x1 A1 x2 B1)  + (C2 x1 A2 x2 B2) + ... + (CR x1 AR x2 BR) where xn is the mode-n product
% A = {A1, A2, ..., AR} is a cell of R matrices Ar, each of size (IxL),
% B = {B1, B2, ..., BR} is a cell of R matrices Br, each of size (JxM),
% C = {C1, C2, ..., CR} is a cell of R tensors Cr, each of size (L,M,K)
%
% INPUTS: 
% - X: tensor of size (IxJxK)  
% - R: number of rank-(L,M,.) terms
% - L: number of columns of the matrices Ar, r=1...R
% - M: number of columns of the matrices Br, r=1...R
%
% OPTIONAL INPUTS:
% - lsearch  (default 'elsr') : The type of line search to use. 
%           The different types are:
%            'none'  (standard ALS Algorithm, i.e., no Line Search is performed)
%             'lsh'  (Line Search proposed by Harshman in, i.e., 
%                     STEP-size set to 1.25) 
%             'lsb'  (Line Search proposed by Bro, i.e. STEP-size set to 
%                     n^1/3, wher n is the iteration index)
%             'elsr' (Exact Line Search with optimal real-valued step)
%             'elsc' (Exact Line Search with optimal complex-valued step)
% - comp     (default='on')     : Compression flag
%             'on' to activate or 'off' to unactivate the dimensionality 
%              reduction pre-processing step
% - Tol1   (default=1e-6)  : threshold value to stop the algorithm in compressed 
%           space if comp=='on' or original space if comp=='off'
% - MaxIt1 (default=5000)  : max number of iterations to stop the algorithm 
%           in compressed space if comp=='on' or original space if comp=='off'
% - Tol2   (default=1e-4)  : threshold value to stop the final refinement stage 
%           in the original space (after decompression), only used if comp=='on'
% - MaxIt2 (default=50)    : max number of iterations to stop the final 
%           refinement stage  (after decompression), only used if comp=='on'
% - Ninit  (default=3)     : number of random starting points used. 
% - A,B,C     : Initial loadings in cell format.
%               If one or more given, loadings used to initialize the algorithm; 
%               in this case, only this initialization is tried (Ninit=1)
%
% OUTPUTS:
% - Estimates: 
%   A : cell of R blocks of size (IxL), 
%   B : cell of R blocks of size (JxM)
%   C : cell of R tensors of size (LxMxK)
% - phi                      : Final value of cost function
% - it1 (it1<= MaxIt1)       : number of iterations in compressed space 
%                              if comp=='on' or original space if comp=='off'. 
% - it2 (it2<= MaxIt2)       : number of iterations in the final refinement 
%                              stage. If comp=='off', it2 is set to 0.
% - phi_vec                  : holds the successive values of phi 
%                              (useful for plot, to see behavior of algorithm)
%  NOTE: all outputs are the ones obtained for the BEST initialization only, the
%  latter being selected as the one that yields the smallest final value of phi.
%-------------------------------------------------------------------------------

%  Reference
% [1] L. De Lathauwer, "Decompositions of a Higher-Order Tensor in Block Terms -
%     Part I: Lemmas for Partitioned Matrices", SIAM J. Matrix Anal. Appl. 
%     (SIMAX), Special Issue on Tensor Decompositions and Applications, Vol. 30,
%     Issue 3, pp. 1022-1032, 2008.
% [2] L. De Lathauwer, "Decompositions of a Higher-Order Tensor in Block Terms -
%     Part II: Definitions and Uniqueness", SIAM J. Matrix Anal. Appl. (SIMAX),
%     Special Issue on Tensor Decompositions and Applications, Vol. 30, Issue 3,
%     pp. 1033-1066, 2008.
% [3] L. De Lathauwer and D. Nion, "Decompositions of a Higher-Order Tensor in 
%     Block Terms�Part III: Alternating Least Squares Algorithms", SIAM J. 
%     Matrix Anal. Appl. (SIMAX), Special Issue on Tensor Decompositions and 
%     Applications, Vol. 30, Issue 3, pp. 1067-1083, 2008.
% [4]  D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%      Complex-Valued Tensor Decompositions. Application in DS-CDMA", Signal 
%      Processing, vol. 88, issue 3, pp. 749-755, March 2008.
% [5]  D. Nion and L. De Lathauwer, "A Block-Component Model Based Blind DS-CDMA
%      Receiver", IEEE Trans. on Signal Processing, Vol. 56, No. 11, 
%      pp. 5567-5579, Nov. 2008.
%-------------------------------------------------------------------------------
% Author: Dimitri Nion   (Feedback: dimitri.nion@gmail.com)
% @Copyright May 2010
% All rights reserved. This M-file and the code in it belongs to the holder 
% of the copyright. For non commercial use only.
%-------------------------------------------------------------------------------

%-----
defopts=struct('R',[],'L',[],'M',[],...
    'lsearch','elsr',...
    'comp','on',...
    'Tol1',1e-6,...
    'MaxIt1',5000,...
    'Tol2',1e-4,...
    'MaxIt2',50,...
    'Ninit',3);
if ~exist('opts','var')
    opts=struct;
end
[R,L,M,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit]=scanparam(defopts,opts);

if isempty(R)||isempty(L)||isempty(M)
    fprintf(2,'R,L,M must be specified.');
    Ycap=[];
    return;
end

Tol1=min(Tol1,Tol2);
Tol2=max(Tol1,Tol2);
MaxIt1=max(MaxIt1,MaxIt2);
MaxIt2=min(MaxIt1,MaxIt2);

X=double(X);
[A,B,C,phi,it1,it2,phi_vec]=bcdLM_alsls(X,R,L,M,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit);
Ycap.A=A;
Ycap.B=B;
Ycap.C=C;
Ycap.phi=phi;
Ycap.it1=it1;
Ycap.it2=it2;
Ycap.phi_vec=phi_vec;
Ycap.type='bcdLM';
end