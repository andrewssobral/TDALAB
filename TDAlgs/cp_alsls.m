function [Y,phi_vec]=cp_alsls(X,opts)
%% function [Y,phi_vec]=cp345_alsls(X,opts,A)
%  This code is based on the distribution at
%             http://perso-etis.ensea.fr/~nion/Codes/Tensor_Decompositions.html
%  An uniform call interface for 3,4,5-way tensors is added for the purpose of TDALAB.
%
%   Below is the original information.
%   [Y,A,phi_vec]=cp345_alsls(X,opts,A) computes a CANDECOMP/PARAFAC decomposition
%   of a 3/4/5-order tensor X in R rank-one terms, stored in the ktensor Y, and cell A
%   consists of the components.
%   The alternating least squares (ALS) algorithm is used, possibly coupled
%   with different choices of line search to speed up convergence.
%
%   If you make use of this code, please cite the following paper:
%      D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%      Complex-Valued Tensor Decompositions. Application in DS-CDMA", Signal 
%      Processing, vol. 88, issue 3, pp. 749-755, March 2008.
%
% MODEL in matrix format:
% Stack the loading vectors in loading matrices
%   A=[a1,a2,...,aR] (IxR)
%   B=[b1,b2,...,bR] (JxR)
%   C=[c1,c2,...,cR] (KxR)
% Consider the three matrix representations of the tensor X
%   X1 (IKxJ) : obtained by stacking all KxJ up-down slices of X one above each other
%   X2 (JIxK) : obtained by stacking all IxK left-right slices of X one above each other
%   X3 (KJxI) : obtained by stacking all JxI front-bottom slices of X one above each other
% Then, the CP model can be written under the three following matrix forms:
%   X1 = (A(kr)C) * B.' + N1 , where N1 is the IKxJ matrix representation of N
%   X2 = (B(kr)A) * C.' + N2 , where N2 is the JIxK matrix representation of N
%   X3 = (C(kr)B) * A.' + N3 , where N3 is the KJxI matrix representation of N
% where kr is Khatri-Rao product (column-wise Kronecker product).
%
% COST FUNCTION:
% This algorithm minimizes the following cost function
%       phi = || X - sum{r=1,..,R} ar o br o cr ||^2,
% which is equivalent to
%       phi = ||X1 - (A(kr)C) * B.'||^2  
%           = ||X2 - (B(kr)A) * C.'||^2 
%           = ||X3 - (C(kr)B) * A.'||^2,
% where A, B and C are estimates of the loading matrices, 
% and || || denotes the Frobenius Norm
%
% ALGORITHM:
% The ALS algorithm consists of the minimization of phi in an alternating way, 
% i.e., at each iteration, minimize phi w.r.t. A, given B and C
% fixed, then update B, given A and C fixed and finally update C, given A and B. 
% See  [1,2,3,4,5]
% To speed up concergence of ALS, the unknowns can be linearly interpolated at 
% each iteration, and the interpolated matrices are used as inputs of the current
% ALS update. Several Line Search techniques are possible [2,6,7,9]
%-------------------------------------------------------------------------------
% INPUTS: 
% - X: tensor of size (IxJxK)
% - R: number of rank-1 terms in the CP3 decomposition
%-------------------------------------------------------------------------------
% OPTIONAL INPUTS:
% - lsearch  (default 'elsr') : The type of line search to use. 
%           The different types are:
%            'none'  (standard ALS Algorithm, i.e., no Line Search is performed)
%             'lsh'  (Line Search proposed by Harshman in [5], i.e., 
%                     STEP-size set to 1.25) 
%             'lsb'  (Line Search proposed by Bro in [2], i.e. STEP-size set to 
%                     n^1/3, wher n is the iteration index)
%             'elsr' (Exact Line Search with optimal real-valued step [6,7])
%             'elsc' (Exact Line Search with optimal complex-valued step [6,7])
% - comp     (default='on')     : Compression flag
%             'on' to activate or 'off' to unactivate the dimensionality 
%              reduction pre-processing step
% - min(Tol)    (default=1e-6)  : threshold value to stop the algorithm in compressed 
%               space if comp=='on' or original space if comp=='off'
% - max(MaxIt)  (default=5000)  : max number of iterations to stop the algorithm 
%               in compressed space if comp=='on' or original space if comp=='off'
% - max(Tol)        (default=1e-4)  : threshold value to stop the final refinement stage 
%               in the original space (after decompression), only used if comp=='on'
%                  -- Invalid for 4,5th-order tensor
% - min(MaxIt) (default=50)    : max number of iterations to stop the final 
%               refinement stage  (after decompression), only used if comp=='on'
%                  -- Invalid for 4,5th-order tensor
% - Ninit      (default=3)     : number of starting points used. 
%              If dimensions allow it, the first initialization is based on DTLD 
%              (Direct Trilinear Decomposition) and the others are random. 
%              Otherwise, all initializations are random.
% - A     :    Initial cell for factors. 
%               If one or more given, matrices used to initialize the algorithm; 
%               in this case, only this initialization is tried (Ninit=1)
%-------------------------------------------------------------------------------
% OUTPUTS:
% - A{1}(IxR), A{2}(JxR), A{3}(KxR)  : estimates of the loading matrices
% - Y                        : Y=ktensor(ones(R,1),A);
% - phi_vec                  : holds the successive values of phi 
%                              (useful for plot, to see behavior of algorithm)
%  NOTE: all outputs are the ones obtained for the BEST initialization only, the
%  latter being selected as the one that yields the smallest final value of phi.
%-------------------------------------------------------------------------------
% STRUCTURE OF THE ALGORITHM
% (A) If comp=='on'
%   (A1) First a dimensionality reduction is performed on X (if dimensions 
%        allow it), in a pre-processing step. This is achieved by truncating 
%        the Multilinear-Singular-Value-Decomposition (MLSVD) of X.
%   (A2) Then, the CP3 model is fitted by ALS (possibly coupled with 
%        Line Search) on the small core tensor resulting from compression, 
%        and we come back to the original space after convergence.
%   (A3) Finally, a few iterations (at most MaxIt2) are performed in the 
%        original space, in a final refinement stage. 
% (B) If comp=='off', then the CP3 model is fitted by ALS 
%     (possibly coupled with Line Search) on the tensor X given as input
%-------------------------------------------------------------------------------
% NOTE ON INITIALIZATION
%    In this function, Ninit initializations are tested.
% (a)- If the dimensions allow it, namely if X has 2 dimensions higher than R, 
%      say I>=R and J>=R, the first initialization is built by exploiting the 
%      tensor itself. The third dimension, say K, is first reduced to 2 by SVD,
%      after which a Generalized EVD technique is applied on the 2 IxJ slices of
%      the matrix pencil (DTLD for Direct Trilinear Dcomposition [8]). 
%      The others (Ninit-1) initializations are all random.
% (b)- Otherwise all Ninit initializations are random
% (c)- If A and/or B and/or C are provided as input arguments to enforce the use
%      of these (this) matrice(s) as starting point(s), then, Ninit is set to 1, 
%      and the provided initialization is the only one used.
%-------------------------------------------------------------------------------
%        REFERENCES
%-------------------------------------------------------------------------------
% [1] R. Bro, "parafac: tutorial and applications", chemom. intell. lab. syst., 
%     vol. 38, pp 149-171, 1997
% [2] R. Bro, "Multi-Way Analysis in the Food Industry: Models, Algorithms, and 
%     Applications", PhD., University of Amsterdam, 1998
% [3] A. Smilde, R. Bro and P. Geladi, "Multi-Way Analysis. Applications in the 
%     Chemical Scinces", John Wiley and Sons, 2004
% [4] N.D. Sidiropoulos, G.B. Giannakis and R. Bro, "Blind PARAFAC Receivers for
%     DS-CDMA Systems", IEEE Trans. Sig. Proc., vol 48, pp 810-823, 2000
% [5] R.A. Harshman, "Foundations of the PARAFAC procedure: Model and 
%     Conditions foran explanatory Multi-mode Factor Analysis", UCLA Working 
%     Papers in Phonetics, vol.16, pp 1-84, 1970
% [6]  M. Rajih, P. Comon and R.A. Harshman, "Enhanced Line Search: A Novel 
%      Method to Accelerate PARAFAC", SIAM Journal Matrix Anal. and Appl.
%      (SIMAX), Volume 30 , Issue 3, pp. 1128-1147, Sept 2008 
% [7]  D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%      Complex-Valued Tensor Decompositions. Application in DS-CDMA", Signal 
%      Processing, vol. 88, issue 3, pp. 749-755, March 2008.
% [8]  E. Leurgans, R.T. Ross, and R.B. Abel, "A decomposition for three-way 
%      arrays",  SIAM J. Matrix Anal. Appl., 14 (1993), pp.~1064--1083.
% [9]  G. Tomasi and R. Bro:  "A comparison of algorithms for fitting the 
%      PARAFAC model.  Computational Statistics & Data Analysis", 
%      50(7): 1700-1734 (2006)
%-------------------------------------------------------------------------------
% Author: Dimitri Nion   (Feedback: dimitri.nion@gmail.com)
% @Copyright May 2010
% All rights reserved. This M-file and the code in it belongs to the holder 
% of the copyright. For non commercial use only.
%-------------------------------------------------------------------------------

%---------------------------------------------------------------------
% Check inputs and define default parameters
%----------------------------------------------------------------
%-----
defoptions = struct('NumOfComp',0,'Tol',[1e-6 1e-4],'MaxIter',[50 5000],'lsearch','lsb',...
    'Ninit',3,'comp','off','Ainit','random');

if ~exist('opts','var')
    opts = struct;
end
[R Tol MaxIt lsearch Ninit comp Ainit]= scanparam(defoptions,opts);

X=double(X);
n=length(size(X));
switch n
    case 3
        fun='cp3_alsls';
    case 4
        fun='cp4_alsls';
    case 5
        fun='cp5_alsls';
    otherwise
        error('Only works for 3,4,5-th order tensors.');
        Y=[];A=[];phi_vec=[];
        return;
end



alsfun=str2func(fun);
if n>3
    MaxIt=max(MaxIt);Tol=min(Tol);
    [Y]=alsfun(X,R,lsearch,Tol,MaxIt,Ninit);
else
    A=[];
    if strcmpi(Ainit,'load')
        load('TDinit.mat','Ainit');
        A=Ainit;
    end

    Tol1=min(Tol);Tol2=max(Tol);
    MaxIt1=max(MaxIt);MaxIt2=min(MaxIt);
    if ~isempty(A)
        [A,B,C,phi,it1,it2,phi_vec]=alsfun(X,R,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A{1},A{2},A{3});
    else
        [A,B,C,phi,it1,it2,phi_vec]=alsfun(X,R,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit);
    end
    Y=ktensor(ones(R,1),A,B,C);
end

end
