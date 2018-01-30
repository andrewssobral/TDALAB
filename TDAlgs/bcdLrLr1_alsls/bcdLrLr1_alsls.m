function [A,B,C,phi,it1,it2,phi_vec]=bcdLrLr1_alsls(X,L_vec,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A,B,C)
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

%---------------------------------------------------------------------
% Check inputs and define default parameters
%---------------------------------------------------------------------
if length(size(X))~=3;error('The input argument X has to be a third-order tensor');end
if nargin<2;error('At least 2 input arguments are required');end
if exist('lsearch')~=1 || isempty(lsearch);    lsearch='elsr';    end
if exist('comp')~=1 || isempty(comp);        comp='on';          end
if exist('Tol1')~=1 || isempty(Tol1);        Tol1=1e-6;       end
if exist('MaxIt1')~=1 || isempty(MaxIt1);    MaxIt1=5e3;      end
if exist('Tol2')~=1 || isempty(Tol2);        Tol2=1e-4;       end
if exist('MaxIt2')~=1 || isempty(MaxIt2);    MaxIt2=5e1;      end
if exist('Ninit')~=1 || isempty(Ninit);      Ninit=3;         end

N=sum(L_vec);
R=length(L_vec);

%-----------------------------------------------------------------------
%  Deal with initialization
%-----------------------------------------------------------------------
init_vec=[0 0 0];  % [0 0 0] if no input matrices is given to initialize
init_given=0;
if exist('A')==1 && isempty(A)~=1
    if (size(A,1)==size(X,1)) && (size(A,2)==N)
    init_vec(1)=1; 
    A_init=A;
    else
    error('The input matrix A does not have the right size')
    end
end
if exist('B')==1 && isempty(B)~=1
    if (size(B,1)==size(X,2)) && (size(B,2)==N)
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
    error('The input matrix B does not have the right size')
    end
end
if sum(init_vec)~=0  % at least one matrix is given
    init_given=1;
    Ninit=1;   
end

% Dimension of the problem + initialization of variables
[I,J,K]=size(X);     % Size of the problem
phi_best=inf;              % Useful to select the best initialization
IndL=[0, cumsum(L_vec)];   % Pattern of the partition
Ps=zeros(N,R);            % Create pattern matrix
for r=1:R
    Ps(IndL(r)+1:IndL(r+1),r)=ones(L_vec(r),1);
end

% STEP1: DIMENSION REDUCTION by truncating the MLSVD
if strcmp(comp,'on')==1
   X_input=X;                 % Store initial input tensor (we will need it in the final refinement stage)
   size_core=[I,J,K];
   if (K*J>N+1) && (I>N); size_core(1)=N+1; end   % first dimension can be reduced
   if (K*I>N+1) && (J>N); size_core(2)=N+1; end   % second dimension can be reduced
   if (J*I>R+1)  && (K>R); size_core(3)=R+1; end      % third dimension can be reduced
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
    
% STEP 2: Fit the BCD-LL1 model with several initializations
for ninit=1:Ninit   % Try max_init different initializations
    % INITIALIZATION
    % Case 1: none of the input matrices A,B,C are given so call generate them
    if init_given==0
        if ninit==1
        [A,B,C]= bcdLrLr1_init(X,L_vec,'gevd');  
        else
        [A,B,C]= bcdLrLr1_init(X,L_vec,'random');  % force random init
        end
    % Case 2: input matrices A and/or B and/or C were given    
    % The ALS loop updates A, then B then C so indeed only B and C have to be 
    % initialized but if e.g. A was provided as input and known to be a good 
    % initialization, it has to be used.
    elseif init_given==1 
        if sum(init_vec)==1     % only one matrix was provided
            [A,B,C]= bcdLrLr1_init(X,L_vec,'ranom');
            if init_vec(1)==1   % only A was provided
               A = A_init;
               B = X1.'/kr_part(A,C,L_vec,ones(1,R)).';  
               C = X2.'/(kr(B,A)*Ps).';
            elseif init_vec(2)==1  % only B was provided
                B = B_init;
                A = X3.'/kr_part(C,B,ones(1,R),L_vec).';
                C = X2.'/(kr(B,A)*Ps).';
            elseif init_vec(3)==1   % only C was provided
                C = C_init;
                A = X3.'/kr_part(C,B,ones(1,R),L_vec).';
                B = X1.'/kr_part(A,C,L_vec,ones(1,R)).';
            end
        elseif sum(init_vec)==2  % two matrices were provided
            if init_vec(1)==0      
               B = B_init; C = C_init;
               A = X3.'/kr_part(C,B,ones(1,R),L_vec).';
               B = X1.'/kr_part(A,C,L_vec,ones(1,R)).';   
               C = X2.'/(kr(B,A)*Ps).';
            elseif init_vec(2)==0      % A and C are given in input and B=[]
               A = A_init; C = C_init;
               B = X1.'/kr_part(A,C,L_vec,ones(1,R)).';   
               A = X3.'/kr_part(C,B,ones(1,R),L_vec).';
               C = X2.'/(kr(B,A)*Ps).';
            elseif init_vec(3)==0  % A and B were given 
               A = A_init; B = B_init;
               C = X2.'/(kr(B,A)*Ps).';  
               A = X3.'/kr_part(C,B,ones(1,R),L_vec).';
               B = X1.'/kr_part(A,C,L_vec,ones(1,R)).';
            end
        elseif sum(init_vec)==3  % three matrices were provided
            A = A_init; B = B_init; C = C_init;
        end
    end  
    A1=A;A2=A;B1=B;B2=B;C1=C;C2=C;  % useful for Line Search

    % LOOP for alternating updates
    phi=norm(X3-kr_part(C,B,ones(1,R),L_vec)*A.','fro');
    stop=0;
    phi_vec=[];
    it1=0;

    while stop==0
        it1=it1+1;
        phi_old=phi;
        % Perform Line Search
        dA=A1-A2;  % search direction for A
        dB=B1-B2;  % search direction for B
        dC=C1-C2;  % search direction for C
        [A,B,C] = bcdLrLr1_lsearch(A2,B2,C2,dA,dB,dC,X1,lsearch,it1,L_vec); 

        % Perform ALS updates
        % Update A
        A=X3.'/kr_part(C,B,ones(1,R),L_vec).';
        % QR factorization for better stability of the algo
          for r=1:R
             [QA RA]=qr(A(:,IndL(r)+1:IndL(r+1)),0);
              A(:,IndL(r)+1:IndL(r+1))=QA;
          end
          % Update B
         B=X1.'/kr_part(A,C,L_vec,ones(1,R)).';
         % Update C
         C=X2.'/(kr(B,A)*Ps).';
         for r=1:R
            C(:,r)=C(:,r)/norm(C(:,r));
         end
         % Calculate the new fit to the model and decide to stop or not
         phi=norm(X3-kr_part(C,B,ones(1,R),L_vec)*A.','fro');
         phi_vec=[phi_vec, phi]; 
         % Define a stop criterion
         if (abs((phi-phi_old)/phi_old)<Tol1) || (it1==MaxIt1) || (phi<Tol1)
            stop=1;
          end
         % Store matrices to prepare next Line Search step
         A2=A1;B2=B1;C2=C1;
         A1=A;B1=B;C1=C;
    end    % Algorithm has converged for this initialization

% Select this init if it is better than previous one
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
            
%    STEP3: COME BACK to ORIGINAL SPACE and PERFORM A FEW MORE ITERATIONS IN ORIGINAL SPACE
it2=0;
if strcmp(comp,'on')==1
    A=U1*A;
    B=U2*B;
    C=U3*C;
    [A,B,C,phi,it2]=bcdLrLr1_alsls(X_input,L_vec,lsearch,'off',Tol2,MaxIt2,[],[],1,A,B,C);
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
function Mat = kr_part(B,C,partB,partC)
%KR_PART Partition-Wise Kronecker product      
[J M]=size(B);
[K N]=size(C);
if (sum(partB)~=M) 
    error(['Error: a matrix with ',int2str(M),' columns can not be partitioned in such a way'])
end
if (sum(partC)~=N) 
    error(['Error: a matrix with ',int2str(N),' columns can not be partitioned in such a way'])
end
if length(partB)~=length(partC)
     error('Error: the 2 input matrices do not have the same number of blocks')
end

indB=[0 cumsum(partB)];
indC=[0 cumsum(partC)];
indMat=[0 cumsum(partB.*partC)];

 Mat=zeros(J*K,sum(partB.*partC));
 for i=1:length(partC)  
     Mat(:,indMat(i)+1:indMat(i+1))=fast_kron( B(:,indB(i)+1:indB(i+1)) , C(:,indC(i)+1:indC(i+1)));
 end
end
         
%*******************************************************************************         
function C = fast_kron (A,B)
%FAST_KRON fast kronecker product
[I,L1]=size(A);
[J,L2]=size(B);

if (L1==1) && (L2==1)
     C=reshape(B*A.',I*J,1);
elseif (L1==1) && (L2>1)
    Bt=B.'; 
    C=reshape(Bt(:)*A.',L2,I*J).';
elseif (L2==1) && (L1>1)
     C=reshape(B*A(:).',I*J,L1);
else
     C=reshape(permute(reshape(B(:)*A(:).',[J,L2,I,L1]),[1 3 2 4]),[I*J,L1*L2]);
end
end

%*******************************************************************************        
function [A,B,C]= bcdLrLr1_init(X,L_vec,init_type)
%BCDLrLr1_INIT Initialization of the loading matrices for the BCD-(Lr,Lr,1)
% INPUTS: 
% - X: 3rd order tensor of size (IxJxK)
% - R: Number of rank-(Lr,Lr,1) components
% - init_type='gevd' to initialize by generalized EVD if possible,
%   otherwise randomly,
% - init_type='random' to enforce random initialization.
% OUTPUTS:
% A(IxN) B(JxN) and C(KxR): Loading Matrices
[I,J,K]=size(X);
N=sum(L_vec);
R=length(L_vec);
IndL=[0, cumsum(L_vec)];  
% Check if I and J are both greater than RL
if min(I,J)<N
    init_type='random';    %  use random init
end

if strcmp(init_type,'gevd')==1
    [A,B,C] = bcdLrLr1_gevd(X,L_vec);
elseif strcmp(init_type,'random')==1
   if isreal(X);
    A=rand(I,N);B=rand(J,N);C=rand(K,R);
   else
    A=rand(I,N)+j*rand(I,N);B=rand(J,N)+j*rand(J,N);C=rand(K,R)+j*rand(K,R);
   end
   % Orthogonalize matrices or submatrices
   if (I>=N);A=orth(A);
   else ; for r=1:R; A(:,IndL(r)+1:IndL(r+1))=orth(A(:,IndL(r)+1:IndL(r+1))); end; end
   if (J>=N);B=orth(B);
   else ; for r=1:R; B(:,IndL(r)+1:IndL(r+1))=orth(B(:,IndL(r)+1:IndL(r+1))); end; end
   if K>=R; C=orth(C); else; C=orth(C');end
end
end

%*******************************************************************************
function [A,B,C] = bcdLrLr1_gevd(X,L_vec)
% BCD_LrLr1_GEVD BCD-(Lr,Lr,1) by generalized EVD of two slices.
[I,J,K]=size(X);
N=sum(L_vec);
R=length(L_vec);
IndL=[0, cumsum(L_vec)];  
% Check if I and J are both greater than N
if min(I,J)<N
       error('bcdLrLr1_gevd:InvalidDimensions',['The input tensor ' ...
          'must have its first two dimensions greater than or equal to N']);
end
% Truncated MLSVD
[U1,U2,U3,X_core] = mlsvd3(X,[N,N,2]);
% Exploit the two slices of X_core to find A (NxN) and B(NxN)
% This is a generalized eigenvalue problem
X1=X_core(:,:,1);
X2=X_core(:,:,2); 
X12=X1/X2;
% Eigenvalue Decomposition
[A S]=eig(X12); 
s=diag(S);
[a b]=sort(abs(s),'descend');
s=s(b(1:N));
S = diag(s);
% In absence of noise, the eigenvalues consist
% of R subsets, each subset consisting of Lr equal eigenvalues, 
% to be associated with one of the R components

% Find A
A = A(:,b(1:N));  % b is the permutation vector that associates the eigenvectors L by L

% If X is real but A is complex (because X12 has one or several pairs of conjugate eigenvalues)
% then transform the columns of A with a similarity matrix to get a real-valued one
if isreal(s)~=1 && isreal(X)==1  
     r=1;
     T=[0.5 -i*0.5;0.5 i*0.5];  % 2 by 2 block: similarity matrix
     while r<=N
          if isreal(s(r))~=1
              A(:,r:r+1)= A(:,r:r+1)*diag(exp(-i*sum(angle(A(1,r:r+1)))/2))*T;   % Produces 2 real columns
              r=r+2;
          else
              r=r+1;
          end
     end
end

% find B  
B = X1.'/A.';
% expand to original space
A=U1*A;
B=U2*B;
% find C in original space
Ps=zeros(N,R);            % Create pattern matrix
for r=1:R
    Ps(IndL(r)+1:IndL(r+1),r)=ones(L_vec(r),1);
end
X2=tens2mat(X,2);
C=X2.'/(kr(B,A)*Ps).';
end    

%*******************************************************************************
function [A,B,C] = bcdLrLr1_lsearch(A,B,C,dA,dB,dC,X,lsearch,it,L_vec)
%BCDLL1_lsearch Line Search for BCD-(Lr,Lr,1)
%   New loading matrices A (IxN), B (JxN), C(KxR) are computed
%   from their previous values and from the search directions dA, dB, dC,
%   as follows:
%  
%      A <- A + mu * dA
%      B <- B + mu * dB
%      C <- C + mu * dC
%
%   Line Search for bcdLrLr1 can for instance be used in gradient-based
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
%      - X: IKxJ matrix unfolding of the observed tensor X 
%      - lsearch: = 'none', 'lsh','lsb', 'elsr' or 'elsc'
%      - it is the iteration step number   
%      - L_vec: the way A and B are partitioned
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


it_start=3;   % Line Search will start after niter_start iterations (at least 2)
if isreal(X) && strcmp(lsearch,'elsc')==1; lsearch='elsr';end
R=length(L_vec);
N=sum(L_vec);

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
        phi=norm(X-kr_part(A,C,L_vec,ones(1,R))*B.','fro');
        % Compute phi with mu=1.25
        mu=1.25;
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
        phim=norm(X-kr_part(Am,Cm,L_vec,ones(1,R))*Bm.','fro');
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
        phi=norm(X-kr_part(A,C,L_vec,ones(1,R))*B.','fro');
       % Compute phi with mu=it^(1/3)
        mu=it^(1/3);
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
        phim=norm(X-kr_part(Am,Cm,L_vec,ones(1,R))*Bm.','fro');
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
        KdAdC=kr_part(dA,dC,L_vec,ones(1,R));
        KAdC=kr_part(A,dC,L_vec,ones(1,R));
        KdAC=kr_part(dA,C,L_vec,ones(1,R));
        KAC=kr_part(A,C,L_vec,ones(1,R));
        Mat3=KdAdC*dB.';
        Mat2=KdAdC*B.' + (KAdC+KdAC)*dB.';
        Mat1=KAC*dB.' + (KAdC+KdAC)*B.';
        Mat0=KAC*B.'-X;
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
    
    
elseif strcmp(lsearch,'elsc')==1
    % Alternate between updates of modulus and argument of mu
     if it<it_start
        mu=1; 
     else
        Tolelsc=1e-4;
        Niterelsc=50;
        KdAdC=kr_part(dA,dC,L_vec,ones(1,R));
        KAdC=kr_part(A,dC,L_vec,ones(1,R));
        KdAC=kr_part(dA,C,L_vec,ones(1,R));
        KAC=kr_part(A,C,L_vec,ones(1,R));
        Mat3=KdAdC*dB.';
        Mat2=KdAdC*B.' + (KAdC+KdAC)*dB.';
        Mat1=KAC*dB.' + (KAdC+KdAC)*B.';
        Mat0=KAC*B.'-X;
        M=[Mat3(:) Mat2(:) Mat1(:) Mat0(:)]; 
        H_mat=M'*M;
        M1=real(H_mat);
        M2=imag(H_mat);
        % Initialization
        mu=1;
        r=abs(mu);
        b=angle(mu);
        t=tan(b/2);
        u=[mu^3 mu^2 mu 1].';
        phi_new=abs(u'*H_mat*u);
        phi0=phi_new;   % initial value of the cost function (with mu=1)
        phi_diff=phi_new;
        it_in=0;
        % Alternate between updates of modulus and angle
        while (phi_diff > Tolelsc) && (it_in < Niterelsc)
            it_in=it_in+1;
            phi_old=phi_new;
        
            % Polynomial expression as a function of r
            h6=M1(1,1);
            h5=2*M1(1,2)*cos(b)+2*M2(1,2)*sin(b);
            h4=M1(2,2)+2*M1(1,3)*cos(2*b)+2*M2(1,3)*sin(2*b);
            h3=2*M1(1,4)*cos(3*b)+2*M1(2,3)*cos(b)+2*M2(1,4)*sin(3*b)+2*M2(2,3)*sin(b);
            h2=M1(3,3)+2*M1(2,4)*cos(2*b)+2*M2(2,4)*sin(2*b);
            h1=2*M1(3,4)*cos(b)+2*M2(3,4)*sin(b);
            h0=M1(4,4);
            pol=[h6 h5 h4 h3 h2 h1 h0];
            pol_der=[6*h6 5*h5 4*h4 3*h3 2*h2 h1];
            sqrts=roots(pol_der);
            sqrts=sqrts(imag(sqrts)==0);
            sqrts=[sqrts;1];      
            extremum=polyval(pol,sqrts);
            r=sqrts(find(extremum==min(extremum),1));
            r=r(1);
            
            % Polynomial expression as a function of t=tan(b/2)
            a1=2*r^3*M1(1,4);
            a2=2*r^4*M1(1,3)+2*r^2*M1(2,4);
            a3=2*r^5*M1(1,2)+2*r^3*M1(2,3)+2*r*M1(3,4);
            a4=r^6*M1(1,1)+r^4*M1(2,2)+r^2*M1(3,3)+ M1(4,4);
            b1=2*r^3*M2(1,4);
            b2=2*r^4*M2(1,3)+2*r^2*M2(2,4);
            b3=2*r^5*M2(1,2)+2*r^3*M2(2,3)+2*r*M2(3,4);
            % Numerator coefficients
            d6=-a1+a2-a3+a4;
            d5=6*b1-4*b2+2*b3;
            d4=15*a1-5*a2-a3+3*a4;
            d3=-20*b1+4*b3;
            d2=-15*a1-5*a2+a3+3*a4;
            d1=6*b1+4*b2+2*b3;
            d0=a1+a2+a3+a4;
             % Denominator
            e6=-3*b1+2*b2-b3;
            e5=-18*a1+8*a2-2*a3;
            e4=45*b1-10*b2-b3;
            e3=60*a1-4*a3;
            e2=-45*b1-10*b2+b3;
            e1=-18*a1-8*a2-2*a3;
            e0=3*b1+2*b2+b3;
            Q=[e6 e5 e4 e3 e2 e1 e0];
            sqrts=roots(Q);
            sqrts=sqrts(imag(sqrts)==0);
            extremum=(1./(1+sqrts.^2).^3).*(d6*sqrts.^6+d5*sqrts.^5+d4*sqrts.^4+...
                  d3*sqrts.^3+d2*sqrts.^2+d1*sqrts+d0*ones(length(sqrts),1));
            b=sqrts(find(extremum==min(extremum)));
            b=b(1);
            b=2*atan(b);

            % updated mu
             mu=r*exp(1i*b);
             t=tan(b/2);
             u=[mu^3 mu^2 mu 1].';
             phi_new=abs(u'*H_mat*u);
             phi_diff=abs(phi_new-phi_old);
        end 
        if phi_new>phi0; mu=1;end       
     end
       A=A+mu*dA;
       B=B+mu*dB;
       C=C+mu*dC;  
end
end
