function [Ycap]=call_bcdLrLr1(X,opts)
%BCDLL1_ALSLS Block Component Decomposition (BCD) in rank-(Lr,Lr,1) terms.
%   [A,B,C]=bcdLrLr1_alsls(X,R,L) computes the BCD-LrLr1 of a third-order tensor.
%   The alternating least squares (ALS) algorithm is used, possibly coupled
%   with different choices of line search to speed up convergence.
%
% DESCRIPTION:
% This function computes the Block Component Decomposition (BCD) of the tensor X
% of size (I,J,K) in a sum of R terms of rank-(Lr,Lr,1): 
% X = (A1 x2 B1 x3 c1)  + (A2 x2 B2 x3 c2) + ... + (AR x2 BR x3 cR) where xn is 
% the mode-n tensor-matrix product
% A = [A1 A2 ... AR] where each block Ar is of size (IxLr),
% B = [B1 B2 ... BR] where each block Br is of size (JxLr),
% C = [c1 c2 ... cR] is of size (K,R)
%
% INPUTS: 
% - X: tensor of size (IxJxK)  
% - L_vec =[L_1 L_2 ... L_R] :  Vector of R values that gives the way A and B 
% are partitioned. (N=sum(L_vec), N is the number of columns of A and B)
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
% - Ninit  (default=3)     : number of starting points used. 
%           If dimensions allow it, the first initialization is based on GEVD 
%           (Generalized EVD) of 2 slices and the others are random. 
%           Otherwise, all initializations are random.
% - A,B,C     : Initial matrices. 
%               If one or more given, matrices used to initialize the algorithm; 
%               in this case, only this initialization is tried (Ninit=1)
%
% OUTPUTS:
% - A (IxN), B(JxN), C(KxR): estimates of the loading matrices
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
% UPDATE RULES: 
% Let X1 (IKxJ), X2 (JIxK) and X3 (KJxI) be 3 matrix representations of X 
% (with left index varying more slowly)
% Then the model in matrix forms is:
% X1=[kron(A1,c1) kron(A2,c2)]*B.'    % useful to update B
% X2=[vec(A1*B1.') vec(A2*B2.')]*C.'  % useful to update C
% X3=[kron(c1,B1) kron(c2,B2)]*A.'    % useful to update A
%-------------------------------------------------------------------------------
% STRUCTURE OF THE ALGORITHM
% (A) If comp=='on'
%    (A1) First, a dimensionality reduction is performed on X 
%         (if dimensions allow it), in a pre-processing step. This is achieved 
%         by truncating the Multilinear-Singular-Value-Decomposition (MLSVD).
%    (A2) Then, the model is fitted by ALS coupled with Line Search on the small 
%         core tensor resulting from compression, and we come back to the 
%         original space after convergence.
%    (A3) Finally, a few iterations (at most MaxIt2) are performed in the 
%        original space, in a final refinement stage.
% (B) If comp=='off', then the model is fitted on the tensor X given as input
%-------------------------------------------------------------------------------
% NOTE ON INITIALIZATION
%    In this function, Ninit initializations are tested.
% (a)- If the dimensions allow it, i.e., if I>=N and J>=N, the first 
%      initialization is built by exploiting the tensor itself. 
%      The third dimension K, is first reduced to 2 by Tucker compression, after
%      which a Generalized EVD technique is applied on the 2 IxJ slices of the 
%      matrix pencil. The others (Ninit-1) initializations are all random.
% (b)- Otherwise all Ninit initializations are random.
% (c)- If A and/or B and/or C are provided as input arguments to enforce the use
%      of these (this) matrice(s) to initialize the algorithm, then Ninit is 
%      set to 1, and the provided initialization is the only one used.

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
% [4] L. De Lathauwer and A. de Baynast, "Blind Deconvolution of DS-CDMA Signals
%     by Means of Decomposition in Rank-(1,L,L) Terms", IEEE Transactions on 
%     Signal Processing, vol. 56, no. 4, Apr. 2008, pp. 1562-1571.
% [5]  D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%      Complex-Valued Tensor Decompositions. Application in DS-CDMA", Signal 
%      Processing, vol. 88, issue 3, pp. 749-755, March 2008.

%-------------------------------------------------------------------------------
% Author: Dimitri Nion   (Feedback: dimitri.nion@gmail.com)
% @Copyright May 2010
% All rights reserved. This M-file and the code in it belongs to the holder 
% of the copyright. For non commercial use only.
%-------------------------------------------------------------------------------

defopts=struct('L_vec',[],...
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
[L_vec,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit]=scanparam(defopts,opts);

if isempty(L_vec)
    fprintf(2,'L_vec must be specified.');
    Ycap=[];
    return;
end

Tol1=min(Tol1,Tol2);
Tol2=max(Tol1,Tol2);
MaxIt1=max(MaxIt1,MaxIt2);
MaxIt2=min(MaxIt1,MaxIt2);

X=double(X);
[A,B,C,phi,it1,it2,phi_vec]=bcdLrLr1_alsls(X,L_vec,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit);
Ycap.A=A;
Ycap.B=B;
Ycap.C=C;
Ycap.phi=phi;
Ycap.it1=it1;
Ycap.it2=it2;
Ycap.phi_vec=phi_vec;
Ycap.type='LrLr1';
end