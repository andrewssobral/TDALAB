function [A,B,C,phi,it1,it2,phi_vec]=cp3_alsls(X,R,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A,B,C)
%CP3_ALSLS CANDECOMP/PARAFAC decomposition of a third-order tensor (CP3).
%   [A,B,C] = cp3_alsls(X,R) computes a CANDECOMP/PARAFAC decomposition
%   of a third-order tensor X in R rank-one terms, stored in the factor 
%   matrices A, B and C, belonging to the first, second and third mode
%   respectively.
%   The alternating least squares (ALS) algorithm is used, possibly coupled
%   with different choices of line search to speed up convergence.
%
%   If you make use of this code, please cite the following paper:
%      D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%      Complex-Valued Tensor Decompositions. Application in DS-CDMA", Signal 
%      Processing, vol. 88, issue 3, pp. 749-755, March 2008.
%
% CP3 MODEL in tensor format:
%
%       _______                                                                         _______
%   K  /      /|         c1/             c2/                       cR/               K /      /|
%     /______/ |          /________   +   /________   +   ...   +   /________    +    /______/ |
%    |  J   |  |   =      |   b1          |   b2                    |   bR           |  J   |  |
%  I |      | /         a1|             a2|                       aR|              I |      | /
%    |______|/            |               |                         |                |______|/
%       X                                                                               N
%
% Given the IxJxK third-order tensor (three-way array) X, cp3_alsls computes its
% CP3 decomposition with R terms, i.e., X is decomposed as a sum of R rank-1 tensors
%               X = sum{r=1,..,R} ar o br o cr + N , 
% where :
% - ar, br and cr, r=1,...,R, are the Ix1, Jx1 and Kx1 loading vectors, respectively,
%   o denotes the outer product, i.e., the element indexed by (i,j,k) of the IxJxK 
%   rank-1 tensor (ar o br o cr) is ar(i)*br(j)*cr(k)
% - N is the IxJxK tensor of residuals 
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
% - Tol1   (default=1e-6)  : threshold value to stop the algorithm in compressed 
%           space if comp=='on' or original space if comp=='off'
% - MaxIt1 (default=5000)  : max number of iterations to stop the algorithm 
%           in compressed space if comp=='on' or original space if comp=='off'
% - Tol2   (default=1e-4)  : threshold value to stop the final refinement stage 
%           in the original space (after decompression), only used if comp=='on'
% - MaxIt2 (default=50)    : max number of iterations to stop the final 
%           refinement stage  (after decompression), only used if comp=='on'
% - Ninit  (default=3)     : number of starting points used. 
%           If dimensions allow it, the first initialization is based on DTLD 
%           (Direct Trilinear Decomposition) and the others are random. 
%           Otherwise, all initializations are random.
% - A,B,C     : Initial matrices. 
%               If one or more given, matrices used to initialize the algorithm; 
%               in this case, only this initialization is tried (Ninit=1)
%-------------------------------------------------------------------------------
% OUTPUTS:
% - A (IxR), B(JxR), C(KxR)  : estimates of the loading matrices
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
%---------------------------------------------------------------------
if ndims(X)~=3;error('The input argument X has to be a third-order tensor');end
if nargin<2;error('At least 2 input arguments are required');end
if exist('lsearch')~=1 || isempty(lsearch);         lsearch='elsr';          end
if exist('comp')~=1   || isempty(comp);             comp='on';              end
if exist('Tol1')~=1 || isempty(Tol1);               Tol1=1e-6;              end
if exist('MaxIt1')~=1 || isempty(MaxIt1);           MaxIt1=5e3;             end
if exist('Tol2')~=1 || isempty(Tol2);               Tol2=1e-4;              end
if exist('MaxIt2')~=1 || isempty(MaxIt2);           MaxIt2=5e1;             end
if exist('Ninit')~=1 || isempty(Ninit);             Ninit=3;                end

%-----------------------------------------------------------------------
%  Deal with initialization matrices
%-----------------------------------------------------------------------
init_vec=[0 0 0];  % [0 0 0] if no input matrices is given to initialize
init_given=0;
if exist('A')==1 && isempty(A)~=1
    if (size(A,1)==size(X,1)) && (size(A,2)==R) 
    init_vec(1)=1; 
    A_init=A;         
    else
    error('The input matrix A does not have the right size')
    end
end
if exist('B')==1 && isempty(B)~=1
    if (size(B,1)==size(X,2)) && (size(B,2)==R)
    init_vec(2)=1;
    B_init=B;
    else
    error('The input matrix B does not have the right size')
    end
end
if exist('C')==1 && isempty(C)~=1
    if (size(C,1)==size(X,3)) && (size(C,2)==R)
    init_vec(3)=1;
    C_init=C;
    else
    error('The input matrix C does not have the right size')
    end
end
if sum(init_vec)~=0  % at least one matrix is given
    init_given=1;
    Ninit=1;          % --> only the initialization provided will be tried
end

% Build matrices Pa,Ja,Pb,Jb useful only for lsearch='elsc'
if strcmp(lsearch,'elsc')==1
    Pa2=fliplr(pascal(4));
    Pa=zeros(4);
    for n=1:4
        Pa(n,n:end)=diag(Pa2,n-1).';
    end     
    Ja=repmat([1i^3,1i^2,1i,1],4,1);

    Pb2=fliplr(pascal(4));
    Pb=zeros(4);
    for n=1:4
        Pb(n,n:end)=diag(Pb2,n-1).';
    end 
    Jb=toeplitz([1,0,0,0],[1,1i,1i^2,1i^3]);
else
    Pa=[];Pb=[];Ja=[];Jb=[];
end
    
% Dimension of the problem + initialization of variables
[I,J,K]=size(X);        % Size of the problem
phi_best=inf;          % Useful to select the best initialization among Ninit


% STEP1: DIMENSION REDUCTION by truncating the MLSVD
if strcmp(comp,'on')==1
   X_input=X;                 % Store initial input tensor (we will need it in the final refinement stage)
   size_core=[I,J,K];
   if I>R+1; size_core(1)=R+1; end             % reduce first dimension
   if J>R+1; size_core(2)=R+1; end
   if K>R+1; size_core(3)=R+1; end
   % Call HOSVD, also named MLSVD (Multilinear Singular Value Decomposition)
   [U1, U2, U3, X] = mlsvd3(X,size_core);    
   if init_given==1 % then compress the matrices given to initialize accordingly
      if init_vec(1)==1;A_init=U1'*A_init;end
      if init_vec(2)==1;B_init=U2'*B_init;end
      if init_vec(3)==1;C_init=U3'*C_init;end
   end
end
 
% Matrix Unfoldings of X
X1=tens2mat(X,1);  % X1 is IKxJ (or smaller dimensions if compression was done)
X2=tens2mat(X,2);  % X2 is JIxK (idem)
X3=tens2mat(X,3);  % X3 is KJxI (idem)
    
%  STEP 2: Fit the CP3 model with several initializations
for ninit=1:Ninit   % Try Ninit different initializations
             
    % INITIALIZATION
    % Case 1: none of the input matrices A,B,C is given
    if init_given==0
       if ninit==1   % first init by DTLD, if possible
          [A,B,C]=cp3_init(X,R,'dtld');
       else
          [A,B,C]=cp3_init(X,R,'random');  % force random init 
       end
    % Case 2: input matrices A and/or B and/or C were given
    % The ALS loop updates A, then B then C so indeed only B and C have to be
    % initialized but if e.g. A was provided as input and known to be a good 
    % initialization, it has to be used.
    elseif init_given==1 
        if sum(init_vec)==1     % only one matrix was provided
           [A,B,C] = cp_init(X,R,'random');  
            if init_vec(1)==1   % only A was provided
               A = A_init;
               B = (inv((A'*A).*(C'*C))*(kat_rao(A,C)'*X1)).';  
               C = (inv((B'*B).*(A'*A))*(kat_rao(B,A)'*X2)).';
            elseif init_vec(2)==1  % only B was provided
                B = B_init;
                A = (inv((C'*C).*(B'*B))*(kat_rao(C,B)'*X3)).';
                C = (inv((B'*B).*(A'*A))*(kat_rao(B,A)'*X2)).';
            elseif init_vec(3)==1   % only C was provided
                C = C_init;
                A = (inv((C'*C).*(B'*B))*(kat_rao(C,B)'*X3)).';
                B = (inv((A'*A).*(C'*C))*(kat_rao(A,C)'*X1)).';
            end
        elseif sum(init_vec)==2  % two matrices were provided
            if init_vec(1)==0    % B and C are given and A=[]
               B = B_init;C = C_init;
               A = (inv((C'*C).*(B'*B))*(kat_rao(C,B)'*X3)).';
            elseif init_vec(2)==0      % A and C are given in input and B=[]
               A = A_init; C = C_init;
               B = (inv((A'*A).*(C'*C))*(kat_rao(A,C)'*X1)).';   
            elseif init_vec(3)==0  % A and B were given 
               A = A_init; B = B_init; 
               C = (inv((B'*B).*(A'*A))*(kat_rao(B,A)'*X2)).';  
            end
        elseif sum(init_vec)==3  % three matrices were provided
            A = A_init; B = B_init; C = C_init;
        end
    end

    A1=A;A2=A;B1=B;B2=B;C1=C;C2=C;  % useful for Line Search

    % LOOP for alternating updates
    phi=norm(X3-kr(C,B)*A.','fro')^2;
    stop=0;
    phi_vec=[];
    it1=0;

    while stop==0
       it1=it1+1;
       phi_old=phi;  
       % Line Search
       dA=A1-A2;  % search direction for A
       dB=B1-B2;  % search direction for B
       dC=C1-C2;  % search direction for C
       [A,B,C] = cp3_lsearch(A2,B2,C2,dA,dB,dC,X3,lsearch,it1,Pa,Ja,Pb,Jb); 
       
       % Perform ALS updates
         A = (inv((C'*C).*(B'*B))*(kr(C,B)'*X3)).';
         B = (inv((A'*A).*(C'*C))*(kr(A,C)'*X1)).';
         C = (inv((B'*B).*(A'*A))*(kr(B,A)'*X2)).';

        % Normalization to avoid overflow
        normA=sqrt(sum(A.*conj(A)));  % Frobenius norm of each column of A
        normB=sqrt(sum(B.*conj(B)));
        normC=sqrt(sum(C.*conj(C)));
        prod_norm=normA.*normB.*normC;
        Scale_mat=diag(prod_norm.^(1/3)); 
        % equal repartition of power of each rank-1 tensor over the 3 vectors:
        A=A*diag(1./normA)*Scale_mat;
        B=B*diag(1./normB)*Scale_mat;
        C=C*diag(1./normC)*Scale_mat;
        % Calculate the new fit to the model and decide to stop or not
        phi=norm(X3-kr(C,B)*A.','fro')^2;
        if nargout==7
           phi_vec=[phi_vec, phi]; 
        end
        % Define a stop criterion
        if (abs((phi-phi_old)/phi_old)<Tol1) || (it1==MaxIt1) || (phi<Tol1)
           stop=1;
        end
        % Store matrices to prepare next Line Search step
         A2=A1;B2=B1;C2=C1;
         A1=A;B1=B;C1=C;
    end    % Algorithm has converged for this initialization

    % Select this initialization if it is better than previous one
     if (phi < phi_best)
        phi_best=phi;
        A_best=A;
        B_best=B;
        C_best=C;
        niter_best=it1;
        phi_vec_best=phi_vec;
        phi_vec=[];
     end

end   % END of ALL init
     
% Output arguments corresponding to the best initialization
A=A_best;
B=B_best;
C=C_best;
it1=niter_best;
phi=phi_best;
phi_vec=phi_vec_best;
            
%  STEP3: COME BACK to ORIGINAL SPACE and PERFORM A FEW MORE ITERATIONS IN ORIGINAL SPACE
it2=0;
if strcmp(comp,'on')==1
    A=U1*A;
    B=U2*B;
    C=U3*C;
    [A,B,C,phi,it2]=cp3_alsls(X_input,R,lsearch,'off',Tol2,MaxIt2,[],[],1,A,B,C);
end            
end

%*******************************************************************************
function [U1,U2,U3,S,S1,S2,S3] = mlsvd3(X,size_core)
%MLSVD3 Multilinear singular value decomposition of a third-order tensor.
[I1,I2,I3]=size(X);
[U1,S1,temp]=svd(reshape(X,I1,I3*I2),'econ'); S1=diag(S1);
[U2,S2,temp]=svd(reshape(permute(X,[2 3 1]),I2,I1*I3),'econ'); S2=diag(S2);
[U3,S3,temp]=svd(reshape(permute(X,[3 1 2]),I3,I2*I1),'econ'); S3=diag(S3);
if nargin==2
    U1=U1(:,1:min(size_core(1),I2*I3));
    U2=U2(:,1:min(size_core(2),I1*I3));
    U3=U3(:,1:min(size_core(3),I1*I2));
end
S=tmprod(tmprod(tmprod(X,U1',1),U2',2),U3',3);
end

%*******************************************************************************
function X_out = tmprod(X,U,mode)
%TMPROD mode-n tensor-matrix product.
[I,J,K]=size(X);
[M,N]=size(U);
if (mode~=1) && (mode~=2) && (mode~=3)
    error('The input variable mode should be 1, 2 or 3')
end
if N~=size(X,mode) 
    error(['The number of columns of the input matrix should be equal to dimension ',int2str(mode),' of the input tensor'])
end    
if mode==1
    X_out = reshape(U*reshape(X,I,J*K) ,M,J,K);      
elseif mode==2
    X_out = permute(reshape (U*reshape(permute(X,[2 1 3]),J,I*K), M,I,K),[2 1 3]);        
elseif mode==3
    X_out = permute(reshape (U*reshape(permute(X,[3 1 2]),K,I*J), M,I,J),[2 3 1]);
end
end

%*******************************************************************************
function [X_mat]=tens2mat(X,mode)
%TENS2MAT Matrix unfoldings of a 3rd order tensor X along a given mode
% INPUTS: - X : tensor of size (IxJxK)
%         - mode = 1 or 2 or 3
% OUTPUTS: X_mat is the matrix unfolding representation of X
% if mode==1:  X_mat is IKxJ  (i=1,...,I, is the slowly varying index)
% if mode==2:  X_mat is JIxK  (j=1,...,J, is the slowly varying index)
% if mode==3   X_mat is KJxI  (k=1,...,K  is the slowly varying index)
[I,J,K]=size(X);
if mode==1
    X_mat=reshape(permute(X,[3 1 2]),I*K,J);
elseif mode==2
    X_mat=reshape(X,J*I,K);
elseif mode==3
    X_mat=reshape(permute(X,[2 3 1]),K*J,I);
else
    error('Input argument mode must be 1, 2 or 3');
end
end

%******************************************************************************
function [A,B,C]=cp3_init(X,R,init_type)
%CP3_INIT Initialization of the loading matrices for the CP3 decomposition of X
% INPUTS: 
% - X: 3rd order tensor of size (IxJxK)
% - R: Number of rank-1 components in the PARAFAC model
% - init_type='dtld' to initialize by Direct Trilinear Decomposition if possible,
%   otherwise randomly,
% - init_type='random' to enforce random initialization.
% OUTPUTS:
% A(IxR) B(JxR) and C(KxR): Loading Matrices

[I,J,K]=size(X);
% Check if 2 dimensions are greater than the rank to see if DTLD can be used
[size_sort,perm_vec]=sort([I J K],'descend');
if ((size_sort(2))<R) 
    init_type='random';    % DTLD can not be used so use random init
end
        
if strcmp(init_type,'dtld')==1
    [A,B,C] = cp3_dtld(X,R);    
elseif strcmp(init_type,'random')==1
    if isreal(X); A=rand(I,R);B=rand(J,R);C=rand(K,R);
    else;A=rand(I,R)+j*rand(I,R);B=rand(J,R)+j*rand(J,R);C=rand(K,R)+j*rand(K,R);
    end
    % Orthogonal matrices if possible else well conditioned matrices
    if I>=R;A=orth(A);else;A=orth(A')';end
    if J>=R;B=orth(B);else;B=orth(B')';end
    if K>=R;C=orth(C);else;C=orth(C')';end
end
end

%*******************************************************************************
function [A,B,C] = cp3_dtld(X,R)
%CP3_DTLD CANDECOMP/PARAFAC by Direct Trilinear Decomposition (DTLD).
%   [A,B,C]=cp3_dtld(X,R) computes the DTLD of a third-order tensor X.
%   DTLD is a method to decompose X as a sum of R rank-1 tensors, i.e., it
%   (non-optimally) fits a rank-R candecomp/parafac (cp) model to X.
%   A, B, C are estimates of the mode-1, mode-2 and mode-3 loading matrices
%   of the cp model, respectively.
%
%   DTLD can only be used if at least two dimensions of the input tensor X 
%   are greater than or equal to the rank R of the decomposition.
%
%   DTLD already gives the exact solution in the case of an exact cp3
%   model, given that the loading matrices in the two long modes are full
%   column rank.
%
%   DTLD computes a truncated mlsvd of X with a core tensor of dimensions 
%   R by R by 2, after which A, B and C are obtained from the generalized 
%   eigenvalue decomposition (GEVD) of the matrix pencil formed by the 
%   two slices of this tensor.
%
%   DTLD can be used as a good starting point for optimization algorithms
%   designed to (optimally) fit the cp model.

%   Copyright 2010
%   Version: 12-07-10
%   Authors: Dimitri Nion (dimitri.nion@gmail.com)
%
%   References:
%   [1] E. Leurgans, R.T. Ross, and R.B. Abel, "A decomposition for 
%       three-way arrays", SIAM J. Matrix Anal. Appl., 14 (1993), 
%       pp. 1064-1083.

[I,J,K]=size(X);
% Check if DTLD is applicable
[size_sort,perm_vec]=sort([I J K],'descend');
if ((size_sort(2))<R) 
    error('cp3_init_dtld:InvalidDimensions',['The input tensor ' ...
          'must have at least two dimensions greater than or equal to R']);
end
% Now Rotate X such that the shortest dimension is in the 3 rd mode
X_perm=permute(X,perm_vec);
[I1,J1,K1]=size(X_perm);
% Truncated MLSVD
[U1,U2,U3,X_core] = mlsvd3(X_perm,[R,R,2]);
% Generalized EVD of the matrix pencil formed by these two slices
[Bti,S]=eig(X_core(:,:,1),X_core(:,:,2));
B=pinv(Bti).';
s=diag(S);
[a ind]=sort(abs(s),'descend');
s=s(ind);
B = B(:,ind);
% If X is real but B is complex (pairs of conjugate eigenvalues)
% then transform B with a similarity matrix to get a real-valued matrix
if isreal(s)~=1 && isreal(X)==1  
     r=1;
     T=[0.5 -1i*0.5;0.5 1i*0.5];  % 2 by 2 block: similarity matrix
     while r<=R
          if isreal(s(r))~=1
              B(:,r:r+1)= B(:,r:r+1)*...
                          diag(exp(-1i*sum(angle(B(1,r:r+1)))/2))*T;   
              r=r+2;
          else
              r=r+1;
          end
     end       
     B=real(B);
end
% find A
A=X_core(:,:,1)/B.';
% expand to the original space
A=U1*A;
B=U2*B;
% find C in original space
C = (((B'*B).*(A'*A))\(kr(B,A)'*reshape(X_perm,I1*J1,K1))).';
% Now place the matrices A, B and C in the right mode from perm_vec
L=cell(1,3);L{1}=A;L{2}=B;L{3}=C;
L(:,perm_vec)=L;
A=L{1};B=L{2};C=L{3};
end

%*******************************************************************************       
function C = kr(A,B)
%KR Khatri-Rao product.
[I R1]=size(A); J=size(B,1); 
C=zeros(I*J,R1);
for j=1:R1
    C(:,j)=reshape(B(:,j)*A(:,j).',I*J,1);
end
end

%*******************************************************************************      
function [A,B,C] = cp3_lsearch(A,B,C,dA,dB,dC,X,lsearch,it,Pa,Ja,Pb,Jb) 
%CP3_LSEARCH Line search for CANDECOMP/PARAFAC order 3.
%   New loading matrices A (IxR), B (JxR), C(KxR) of cp3 are computed
%   from their previous values and from the search directions dA, dB, dC,
%   as follows:
%  
%      A <- A + mu * dA
%      B <- B + mu * dB
%      C <- C + mu * dC
%
%   Line Search for cp3 can for instance be used in gradient-based
%   algorithms or with Alternating Least Squares (ALS), in order to speed 
%   up convergence.
%
%   For instance, if Line Search is used with ALS, the search directions
%   may be defined outside this function by dA=A1-A, dB=B1-B, dC=C1-C and
%   the input matrices by A=A1, B=B1 and C=C1 where (A1,B1,C1) denote
%   conditional least squares estimates at the current ALS iteration it,
%   whereas (A,B,C) denote old estimates, such that
% 
%      A <- A + mu*(A1-A)
%      B <- B + mu*(B1-B)
%      C <- C + mu*(C1-C)
%
%   This means that, in the context of ALS, setting mu=1 annihilates the
%   line search, since A, B and C will be set to their current values A1,
%   B1 and C1.
%
%   For instance, if Line Search is used with gradient-based algorithms,
%   the search directions dA, dB, dC, could be (- gradient) of the cost 
%   function wrt A, B and C, respectively, in a gradient descent 
%   algorithm, or other more sophisticated search directions like in a
%   conjugate gradient algorithm.
%
%   Several choices of mu are possible:
%
%      - if lsearch = 'none' then no interpolation is done (in the sense of 
%        lsearch for als), i.e. we enforce mu=1.
%      - if lsearch = 'lsh' (Line Search proposed by Harshman), then mu is 
%        fixed to 1.25 and if the interpolated matrices do not decrease 
%        the cost function, then mu=1 is enforced.
%      - if lsearch = 'lsb' (Line Search proposed by Bro), then mu is fixed
%        to it^(1/3) and if the interpolated matrices do not decrease 
%        the cost function, then mu=1 is enforced.
%      - if lsearch = 'elsr' (Exact Line Search with real step), then we 
%        seek for the optimal real-valued step mu that minimizes the cost 
%        function. This amounts to minimizing a polynomial of degree 6 in
%        mu. It can also be used for complex-valued data.
%      - if lsearch = 'elsc' (Exact Line Search with complex step), then we
%        seek for the optimal complex-valued step mu that minimizes the
%        cost function. If the data are real-valued, the step has to be
%        real-valued, and in this case elsc is replaced by elsr.
%
%   INPUTS: 
%   
%      - A,B,C: estimates of the loading matrices A,B,C at it-1
%      - dA,dB,dC: search directions
%      - X: KJxI matrix unfolding of the observed tensor X 
%      - lsearch: = 'none', 'lsh','lsb', 'elsr' or 'elsc'
%      - it is the iteration step number       
%      - Pa,Ja,Pb,Jb: matrices used only by elsc in order to alternate
%        between updates of real(mu) and imag(mu) in an efficient way.
%       These matrices are generated outside this function as follows:
%           Pa2=fliplr(pascal(4));Pa=zeros(4);
%    		for n=1:4;Pa(n,n:end)=diag(Pa2,n-1).';end     
%    		Ja=repmat([1i^3,1i^2,1i,1],4,1);
%    		Pb2=fliplr(pascal(4));Pb=zeros(4);
%           for n=1:4;Pb(n,n:end)=diag(Pb2,n-1).';end 
%           Jb=toeplitz([1,0,0,0],[1,1i,1i^2,1i^3]);
%
%   OUTPUTS: 
%
%      - Updated loading matrices A,B,C

%   Copyright 2010
%   Version: 09-07-10
%   Authors: Dimitri Nion (dimitri.nion@gmail.com),
%
%   References:
%   [1] R.A. Harshman, "Foundations of the PARAFAC procedure: Model and 
%       Conditions for an explanatory Multi-mode Factor Analysis", UCLA 
%       Working Papers in Phonetics, vol.16, pp 1-84, 1970.
%   [2] R. Bro, "Multi-Way Analysis in the Food Industry: Models, 
%       Algorithms, and Applications", PhD. dissertation, University of 
%       Amsterdam, 1998.
%   [3] M. Rajih, P. Comon and R.A. Harshman, "Enhanced Line Search: A 
%       Novel Method to Accelerate  PARAFAC", SIAM Journal Matrix Anal. and
%       Appl. (SIMAX), Volume 30 , Issue 3, pp. 1128-1147, Sept 2008.
%   [4] D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%       Complex-Valued Tensor Decompositions. Application in DS-CDMA", 
%       Signal Processing, vol. 88, issue 3, pp. 749-755, March 2008.

it_start=3; % Line Search will start after it_start iterations (at least 2)
if isreal(X) && strcmp(lsearch,'elsc')==1; lsearch='elsr';end

%   The cp3 least-squares cost function is  
%   phi = norm (X - kr(C,B)*A.' , 'fro')^2, 
%   where kr is the Khatri-Rao product and X is the mode-1 KJ by I matrix 
%   unfolding of the IxJxK tensor that is decomposed. Note that, in the 
%   row indexing of X, the index k=1,..,K has to vary more slowly than 
%   j=1,..,K, since the cost function involves kr(C,B) and not kr(B,C).
%
%   The crucial point is to choose a good step size mu such that
%      phinew = norm (X - kr((C+mu*dC),(B+mu*dC))*(A+mu*dA).' , 'fro')^2,   
%   is smaller than phi.

% lsearch='none', i.e., standard ALS            
if strcmp(lsearch,'none')==1
     mu=1; 
     A=A+mu*dA;
     B=B+mu*dB;
     C=C+mu*dC;

% lsearch='lsh'
elseif strcmp(lsearch,'lsh')==1
    if it<it_start
        mu=1;
        A=A+mu*dA;
        B=B+mu*dB;
        C=C+mu*dC;   
    else
        % Compute phi with mu=1
        A=A+dA;
        B=B+dB;
        C=C+dC;
        phi=norm(X-kr(C,B)*A.','fro');
        % Compute phi with mu=1.25
        mu=1.25;
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
        phim=norm(X-kr(Cm,Bm)*Am.','fro');
        % Accept or Reject Interpolation                   
        if phim < phi     % accept
            A=Am;
            B=Bm;
            C=Cm;              
        end
    end   
% lsearch='lsb'         
elseif strcmp(lsearch,'lsb')==1
    if it<it_start
        mu=1;
        A=A+mu*dA;
        B=B+mu*dB;
        C=C+mu*dC;  
    else
        % Compute phi with mu=1
        A=A+dA;
        B=B+dB;
        C=C+dC;
        phi=norm(X-kr(C,B)*A.','fro');
        % Compute phi with mu=it^(1/3)
        mu=it^(1/3);
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
        phim=norm(X-kr(Cm,Bm)*Am.','fro');
        % Accept or Reject Interpolation                   
        if phim < phi     % accept
            A=Am;
            B=Bm;
            C=Cm;              
        end
    end
% lsearch='elsr', compute optimal real-valued mu
elseif strcmp(lsearch,'elsr')==1
    if it<it_start
        mu=1;  
    else
        KdCdB=kr(dC,dB);
        KdCB=kr(dC,B);
        KCdB=kr(C,dB);
        KCB=kr(C,B);
        Mat3=KdCdB*dA.';
        Mat2=KdCdB*A.' + (KdCB+KCdB)*dA.';
        Mat1=KCB*dA.' + (KdCB+KCdB)*A.';
        Mat0=KCB*A.'-X;
        M=[Mat3(:) Mat2(:) Mat1(:) Mat0(:)]; 
        H_mat=real(M'*M);  
        % Now we define the coefficients of the 6th order polynomial
        d6=H_mat(1,1);
        d5=2*H_mat(1,2);
        d4=2*H_mat(1,3)+H_mat(2,2);
        d3=2*(H_mat(1,4)+H_mat(2,3));
        d2=2*H_mat(2,4)+H_mat(3,3);
        d1=2*H_mat(3,4);
        d0=H_mat(4,4);
        pol=[d6 d5 d4 d3 d2 d1 d0];
        pol_der=[6*d6 5*d5 4*d4 3*d3 2*d2 d1];
        sqrts=roots(pol_der);
        sqrts=sqrts(imag(sqrts)==0); % real roots
        sqrts=[sqrts;1];   
        % Choice of optimal mu
        extremum=polyval(pol,sqrts);
        mu=sqrts(find(extremum==min(extremum),1));
        mu=mu(1);
    end
    A=A+mu*dA;
    B=B+mu*dB;
    C=C+mu*dC;   
% lsearch='elsc'
% Alternate between updates of real and imaginary part of mu
elseif strcmp(lsearch,'elsc')==1
    if it<it_start
        mu=1; 
    else
        Tolelsc=1e-4;
        Niterelsc=50;
        KdCdB=kr(dC,dB);
        KdCB=kr(dC,B);
        KCdB=kr(C,dB);
        KCB=kr(C,B);
        Mat3=KdCdB*dA.';
        Mat2=KdCdB*A.' + (KdCB+KCdB)*dA.';
        Mat1=KCB*dA.' + (KdCB+KCdB)*A.';
        Mat0=KCB*A.'-X;
        M=[Mat3(:) Mat2(:) Mat1(:) Mat0(:)]; 
        H_mat=M'*M;  
        % Initialization
        mu=1;
        a=1;  % real part of step
        u=[mu^3 mu^2 mu 1].';
        phi_new=abs(u'*H_mat*u);
        phi0=phi_new;   % initial value of the cost function (with mu=1)
        phi_diff=phi_new;
        it_in=0;
        % Alternate between updates of real part a and imag part b of mu 
        while (phi_diff > Tolelsc) && (it_in < Niterelsc)
            it_in=it_in+1;
            phi_old=phi_new;
            % Update imag part b      
            Da=Pa.*Ja.*toeplitz([1,0,0,0],[1 a a^2 a^3]);
            Ha_mat=real(Da'*H_mat*Da);
            % Now we define the coefficients of the 6th order
            % polynomial in b
            d6=Ha_mat(1,1);
            d5=Ha_mat(1,2)+Ha_mat(2,1);
            d4=Ha_mat(1,3)+Ha_mat(2,2)+Ha_mat(3,1);
            d3=Ha_mat(1,4)+Ha_mat(2,3)+Ha_mat(4,1)+Ha_mat(3,2);
            d2=Ha_mat(2,4)+Ha_mat(3,3)+Ha_mat(4,2);
            d1=Ha_mat(3,4)+Ha_mat(4,3);
            d0=Ha_mat(4,4);
            pol=[d6 d5 d4 d3 d2 d1 d0];
            pol_der=[6*d6 5*d5 4*d4 3*d3 2*d2 d1];
            sqrts=roots(pol_der);
            sqrts=sqrts(imag(sqrts)==0);
            sqrts=[sqrts;1];      
            % Choice of optimal b
            extremum=polyval(pol,sqrts);
            b=sqrts(find(extremum==min(extremum),1));
            b=b(1);

            % Update real part a                   
            Db=Pb.*Jb.*toeplitz([1,0,0,0],[1 b b^2 b^3]);
            Hb_mat=real(Db'*H_mat*Db);
            % Now we define the coefficients of the 6th order
            % polynomial in a
            d6=Hb_mat(1,1);
            d5=Hb_mat(1,2)+Hb_mat(2,1);
            d4=Hb_mat(1,3)+Hb_mat(2,2)+Hb_mat(3,1);
            d3=Hb_mat(1,4)+Hb_mat(2,3)+Hb_mat(4,1)+Hb_mat(3,2);
            d2=Hb_mat(2,4)+Hb_mat(3,3)+Hb_mat(4,2);
            d1=Hb_mat(3,4)+Hb_mat(4,3);
            d0=Hb_mat(4,4);
            pol=[d6 d5 d4 d3 d2 d1 d0];           
            pol_der=[6*d6 5*d5 4*d4 3*d3 2*d2 d1];
            sqrts=roots(pol_der);
            sqrts=sqrts(imag(sqrts)==0);
            sqrts=[sqrts;1];   % we add 1 in the set of possible values    
            % Choice of optimal a
            extremum=polyval(pol,sqrts);
            a=sqrts(find(extremum==min(extremum),1));
            a=a(1);

            % update mu
            mu=a+1i*b;
            u=[mu^3 mu^2 mu 1].';
            phi_new=abs(u'*H_mat*u);
            phi_diff=abs(phi_new-phi_old);         
        end
        if phi_new>phi0; mu=1;end
    end
    % Interpolate
    A=A+mu*dA;
    B=B+mu*dB;
    C=C+mu*dC;
end
end