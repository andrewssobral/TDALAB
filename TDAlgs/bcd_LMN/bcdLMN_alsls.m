function [A,B,C,D,phi,it1,it2,phi_vec]=bcdLMN_alsls(X,R,L,M,N,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A,B,C,D)
%BCDLMN_ALSLS Block Component Decomposition (BCD) in rank-(L,M,N) terms.
%   [A,B,C,D]=bcdLMN_alsls(X,R,L,M,N) computes the BCD-(L,M,N) of a third-order tensor.
%   The alternating least squares (ALS) algorithm is used, possibly coupled
%   with different choices of line search to speed up convergence.
%
% DESCRIPTION:
% This function computes the Block Component Decomposition (BCD) of the tensor X
% of size (I,J,K) in a sum of R terms of rank-(L,M,N): 
% X = (D1 x1 A1 x2 B1 x3 C1)  + (D2 x1 A2 x2 B2 x3 C2) + ... + (DR x1 AR x2 BR x3 CR) where xn is the mode-n product
% A = {A1, A2, ..., AR} is a cell of R matrices Ar, each of size (IxL),
% B = {B1, B2, ..., BR} is a cell of R matrices Br, each of size (JxM),
% C = {C1, C2, ..., CR} is a cell of R matrices Cr, each of size (KxN),
% D = {D1, D2, ..., DR} is a cell of R core tensors Dr, each of size (LxMxN),
%
% INPUTS: 
% X: tensor of size (IxJxK)  
% R: number of rank-(L,M,.) components
% L: number of columns of the matrices Ar, r=1...R
% M: number of columns of the matrices Br, r=1...R
% N: number of columns of the matrices Cr, r=1...R
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
% - A,B,C,D   : Initial loadings in cell format.
%               If one or more given, loadings used to initialize the algorithm; 
%               in this case, only this initialization is tried (Ninit=1)
%
% OUTPUTS:
% - Estimates: 
%   A : cell of R blocks of size (IxL), 
%   B : cell of R blocks of size (JxM)
%   C : cell of R blocks of size (KxN)
%   D : cell of R tensors of size (LxMxN)
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
if nargin<5;error('At least 5 input arguments are required');end
if exist('lsearch')~=1 || isempty(lsearch);    lsearch='elsr';    end
if exist('comp')~=1 || isempty(comp);        comp='on';          end
if exist('Tol1')~=1 || isempty(Tol1);        Tol1=1e-6;       end
if exist('MaxIt1')~=1 || isempty(MaxIt1);    MaxIt1=5e3;      end
if exist('Tol2')~=1 || isempty(Tol2);        Tol2=1e-4;       end
if exist('MaxIt2')~=1 || isempty(MaxIt2);    MaxIt2=5e1;      end
if exist('Ninit')~=1 || isempty(Ninit);      Ninit=3;         end

%-----------------------------------------------------------------------
%  Deal with initialization
%-----------------------------------------------------------------------
init_vec=[0 0 0 0];  % [0 0 0 0] if no loading is given to initialize
init_given=0;
if exist('A')==1 && isempty(A)~=1
    if (size(A{1},1)==size(X,1)) && (size(A,2)==R)
    init_vec(1)=1; 
    A_init=A;
    else
    error('The input cell A does not have the right size')
    end
end
if exist('B')==1 && isempty(B)~=1
    if (size(B{1},1)==size(X,2)) && (size(B,2)==R)
    init_vec(2)=1;
    B_init=B;
    else
    error('The input cell B does not have the right size')
    end
end
if exist('C')==1 && isempty(C)~=1
    if (size(C{1},1)==size(X,3)) && (size(C,2)==R)
    init_vec(3)=1;
    C_init=C;
    else
    error('The input cell C does not have the right size')
    end
end
if exist('D')==1 && isempty(D)~=1
    if (( size(D,2)==R && size(D,1)==1) || (size(D,1)==R && size(D,2)==1)) 
    init_vec(4)=1;
    D_init=D;
    else
    error('The input cell D does not have the right size')
    end
end
if sum(init_vec)~=0  % at least one matrix is given
    init_given=1;
    Ninit=1;   
end

% Dimension of the problem + initialization of variables
[I,J,K]=size(X);     % Size of the problem
phi_best=inf;              % Useful to select the best initialization
RL=R*L;                    % Number of columns of A 
RM=R*M;                    % Number of columns of B
RN=R*N;
RLMN=R*L*M*N;
L_vec=repmat(L,1,R);       % Pattern of the blocks in A: L columns for each sub-matrix
M_vec=repmat(M,1,R);       % Pattern of the blocks in B: M columns for each sub-matrix
N_vec=repmat(N,1,R);       % Pattern of the blocks in C: N columns for each sub-matrix

% STEP1: DIMENSION REDUCTION by truncating the MLSVD
if strcmp(comp,'on')==1
   X_input=X;                 % Store initial input tensor (we will need it in the final refinement stage)
   size_core=[I,J,K];
    if min(K*J,I)>=min(R*L,R*N*M); size_core(1)=min(R*L,R*N*M); end        
    if min(I*K,J)>=min(R*M,R*L*N); size_core(2)=min(R*M,R*L*N); end    
    if min(J*I,K)>=min(R*N,R*M*L); size_core(3)=min(R*N,R*M*L); end      
     % Call HOSVD, also named MLSVD (Multilinear Singular Value Decomposition)
   [U1, U2, U3, X] = mlsvd3(X,size_core);  
   [I,J,K]=size(X);   % we now work with the new size of X 
   if init_given==1 % then compress the matrices given to initialize accordingly
        if init_vec(1)==1
          A_init=mat2cell(U1'*cell2mat(A_init),I,L_vec);
        end
        if init_vec(2)==1
          B_init=mat2cell(U2'*cell2mat(B_init),J,M_vec);        
        end
        if init_vec(3)==1
          C_init=mat2cell(U3'*cell2mat(C_init),K,N_vec);        
        end     
   end
end

% Matrix Unfoldings of X
X1=tens2mat(X,1);  % X1 is IKxJ (or smaller dimensions if compression was done)
X2=tens2mat(X,2);  % X2 is JIxK (idem)
X3=tens2mat(X,3);  % X3 is KJxI (idem)
X_vec=X3(:);       % X_vec is IKJx1
    
% STEP 2: Fit the decomposition with several initializations
for ninit=1:Ninit   % Try max_init different initializations
    % INITIALIZATION
    % Case 1: none of the input cells A,B,C,D is given so call generate them
    if init_given==0
       [A_cell,B_cell,C_cell,D_cell]= bcdLMN_init(X,L_vec,M_vec,N_vec);   % random initialization in cell format
           
    % Case 2: input A and/or B and/or C and/or D were given    
    elseif init_given==1 
        if (sum(init_vec)==1) || (sum(init_vec)==2)     % only one or two matrices were provided
            [A_cell,B_cell,C_cell,D_cell]= bcdLMN_init(X,L_vec,M_vec,N_vec);  % generate random initializations and then overwrite by the one given in input
            if init_vec(1)==1   % A was provided
               A_cell = A_init; 
               [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
               [C_cell,C]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec);               
               [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);
            elseif init_vec(2)==1  % B was provided
                B_cell = B_init; 
                [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec);
                [C_cell,C]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec); 
                [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);
            elseif init_vec(3)==1   % C was provided
                C_cell = C_init; 
                [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec);
                [B_cell,B]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
                [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);
            elseif init_vec(4)==1   % D was provided
                D_cell = D_init; 
                [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec);
                [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
                [C_cell,C,D_cell]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec); 
            end
        elseif sum(init_vec)==3  % three matrices were provided
            if init_vec(1)==0    % B,C and D are given and A=[]
               B_cell=B_init;C_cell=C_init;D_cell=D_init;
               [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec); 
               [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
               [C_cell,C]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec); 
               [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);

            elseif init_vec(2)==0      % A,C and D are given in input and B=[]
               A_cell=A_init;C_cell=C_init;D_cell=D_init; 
               [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);   % so we update B to exploit A and C
               [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec); 
               [C_cell,C]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec);  
               [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);

            elseif init_vec(3)==0  % A, B and D were given 
               A_cell=A_init;B_cell=B_init;D_cell=D_init;
               [C_cell,C,D_cell]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec);   % so we update C first to exploit A and B
               [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec); 
               [B_cell,B]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
               [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);

            elseif init_vec(4)==0  % A, B and C were given
               A_cell=A_init;B_cell=B_init;C_cell=C_init;
               [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);
               [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec); 
               [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
               [C_cell,C,D_cell]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec);   

            end
        elseif sum(init_vec)==4   % all matrices were given
            A_cell=A_init;B_cell=B_init;C_cell=C_init;D_cell=D_init;
        end
    end  
            
    A=cell2mat(A_cell);
    B=cell2mat(B_cell);
    C=cell2mat(C_cell);
     
    % Initialization of matrices useful for Line Search
    A1_cell=A_cell;A2_cell=A_cell;
    B1_cell=B_cell;B2_cell=B_cell;
    C1_cell=C_cell;C2_cell=C_cell;  

    % LOOP for alternating updates
    phi=norm(X3-kr_part(C,B,N_vec,M_vec)*blockdiag(tens2mat_block(D_cell,3),R)*A.','fro');
    stop=0;
    phi_vec=[];
    it1=0;
   
    while stop==0
        it1=it1+1;
        phi_old=phi;
        % Perform Line Search
        dA_cell=mat2cell(cell2mat(A1_cell)-cell2mat(A2_cell),I,L_vec);  % search direction for A
        dB_cell=mat2cell(cell2mat(B1_cell)-cell2mat(B2_cell),J,M_vec);  % search direction for B
        dC_cell=mat2cell(cell2mat(C1_cell)-cell2mat(C2_cell),K,N_vec);  % search direction for C
        [A_cell,B_cell,C_cell] = bcdLMN_lsearch(A2_cell,B2_cell,C2_cell,dA_cell,dB_cell,dC_cell,D_cell,X2,lsearch,it1,L_vec,M_vec,N_vec); 
        % Perform ALS updates
        % Update A 
        [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec);
        % Update  B 
        [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec);
        % Update  C
        [C_cell,C]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec);      
        % Update  D
        [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec);      
        
        % Calculate the new fit to the model and decide to stop or not
         phi=norm(X3-kr_part(C,B,N_vec,M_vec)*blockdiag(tens2mat_block(D_cell,3),R)*A.','fro');
         phi_vec=[phi_vec, phi]; 
         % Define a stop criterion
         if (abs((phi-phi_old)/phi_old)<Tol1) || (it1==MaxIt1) || (phi<Tol1)
            stop=1;
          end
        % Store matrices to prepare next Line Search step
        A2_cell=A1_cell;B2_cell=B1_cell;C2_cell=C1_cell;
        A1_cell=A_cell;B1_cell=B_cell;C1_cell=C_cell;
       
    end    % Algorithm has converged for this initialization

% Select this init if it is better than previous one
     if (phi < phi_best)
        phi_best=phi;
        A_best=A;
        B_best=B;
        C_best=C;
        D_best=D_cell;
        niter_best=it1;
        phi_vec_best=phi_vec;
        phi_vec=[];
     end
end   % END of ALL init
     
% Output arguments corresponding to the best initialization
A=A_best;
B=B_best;
C=C_best;
D_cell=D_best;
it1=niter_best;
phi=phi_best;
phi_vec=phi_vec_best;
            
%    STEP3: COME BACK to ORIGINAL SPACE and PERFORM A FEW MORE ITERATIONS IN ORIGINAL SPACE
it2=0;
if strcmp(comp,'on')==1
    A=U1*A;
    B=U2*B;
    C=U3*C;
end
A_cell=mat2cell(A,size(A,1),L_vec);
B_cell=mat2cell(B,size(B,1),M_vec);
C_cell=mat2cell(C,size(C,1),N_vec);

if strcmp(comp,'on')==1
    [A_cell,B_cell,C_cell,D_cell,phi,it2]=bcdLMN_alsls(X_input,R,L,M,N,lsearch,'off',Tol2,MaxIt2,[],[],1,A_cell,B_cell,C_cell,D_cell);
end
A=A_cell;
B=B_cell;
C=C_cell;
D=D_cell;
end

%*******************************************************************************    
function [A_cell,A,D_cell]=update_A(B_cell,C_cell,D_cell,X3,L_vec,M_vec,N_vec)
 % Least Squares Update of A, given B, C and D
R=length(L_vec);
D_mat = tens2mat_block(D_cell,3);     % size  MLxRN
B = cell2mat(B_cell);
C = cell2mat(C_cell);
A = X3.'/(kr_part(C,B,N_vec,M_vec)*blockdiag(D_mat,R)).';
A_cell=mat2cell(A,size(A,1),L_vec);  
for r=1:R;[A_cell{r} RA_cell{r}]=qr(A_cell{r},0);end
A=cell2mat(A_cell);
if nargout==3  % Incorporate the RA_cell matrices in D_cell
 for r=1:R
     [Lr Mr,Nr]=size(D_cell{r});
     D_cell{r}=mat2tens(tens2mat(D_cell{r},3)*RA_cell{r}.',3,[Lr,Mr,Nr]);
 end
end
end

%*******************************************************************************
function [B_cell,B,D_cell]=update_B(A_cell,C_cell,D_cell,X1,L_vec,M_vec,N_vec)   
R=length(L_vec);
D_mat = tens2mat_block(D_cell,1);     % size  LNxRM
A = cell2mat(A_cell);
C = cell2mat(C_cell);
B = X1.'/(kr_part(A,C,L_vec,N_vec)*blockdiag(D_mat,R)).';
B_cell=mat2cell(B,size(B,1),M_vec);  
for r=1:R;[B_cell{r} RB_cell{r}]=qr(B_cell{r},0);end
B=cell2mat(B_cell);
if nargout==3 % Incorporate the RB_cell matrices in D_cell
 for r=1:R
     [Lr Mr,Nr]=size(D_cell{r});
     D_cell{r}=mat2tens(tens2mat(D_cell{r},1)*RB_cell{r}.',1,[Lr,Mr,Nr]);
 end
end
end
     
%*******************************************************************************    
function [C_cell,C,D_cell]=update_C(A_cell,B_cell,D_cell,X2,L_vec,M_vec,N_vec)
R=length(L_vec);
D_mat = tens2mat_block(D_cell,2);     % size  MLxRN
A = cell2mat(A_cell);
B = cell2mat(B_cell);
C = X2.'/(kr_part(B,A,M_vec,L_vec)*blockdiag(D_mat,R)).';
C_cell=mat2cell(C,size(C,1),N_vec);  
for r=1:R;[C_cell{r} RC_cell{r}]=qr(C_cell{r},0);end
C=cell2mat(C_cell);
if nargout==3  % Incorporate the RC_cell matrices in D_cell
 for r=1:R
     [Lr Mr,Nr]=size(D_cell{r});
     D_cell{r}=mat2tens(tens2mat(D_cell{r},2)*RC_cell{r}.',2,[Lr,Mr,Nr]);
 end     
end
end  

%*******************************************************************************    
function [D_cell]=update_D(A_cell,B_cell,C_cell,X_vec,L_vec,M_vec,N_vec)
    
% X_vec is of size IKJx1
R = length(L_vec);
A=cell2mat(A_cell);
B=cell2mat(B_cell);
C=cell2mat(C_cell);
D = kr_part(kr_part(A,C,L_vec,N_vec),B,L_vec.*N_vec,M_vec)\X_vec; % size RLNMx1
LMN_vec=[0 L_vec].*[0 M_vec].*[0 N_vec];
D_cell=cell(1,R);
for r=1:R
  Dr_vec=D(sum(LMN_vec(1:r))+1:sum(LMN_vec(1:r+1)));  % size LNMx1
  D_cell{r}=permute(reshape(Dr_vec,M_vec(r),N_vec(r),L_vec(r)),[3 1 2]);
end
end

%*******************************************************************************
function C_mat = tens2mat_block(C,mode)
% This function rewrites a cell C of R tensors C{r}=Cr (LxMxN) into a matrix C_mat
% of size: 
% - L*N rows by R*M columns if mode==1
% - M*L rows by R*N columns if mode==2
% - N*M rows by R*L columns if mode==3
if min(size(C))~=1
   error('The cell C must have one dimension equal to 1')
end
R=max(size(C));    % The R tensors are stacked in the cell C (either in column of row format)

if mode==1
   for r=1:R
     D{1,r}=tens2mat(C{r},1);   % This transforms LxMxN tensor into LNxM matrix
   end
   C_mat=cell2mat(D);           % Of size LNxRM

elseif mode==2
   for r=1:R
     D{1,r}=tens2mat(C{r},2);   % This transforms LxMxN tensor into MLxN matrix
   end
   C_mat=cell2mat(D);           % Of size MLxRN

elseif mode==3
   for r=1:R
     D{1,r}=tens2mat(C{r},3);   % This transforms LxMxN tensor into NMxL matrix
   end
   C_mat=cell2mat(D);           % Of size NMxRL
end
end

%*******************************************************************************
function D_block=blockdiag(D_mat,R)
% This function creates a block-diagonal matrix from R sumatrices 
% The input matrix D_mat is a stack of these R submatrices on a row
%Each submatrix has size H*L
H=size(D_mat,1);
L=size(D_mat,2)/R;

for r=1:R
    D_block((r-1)*H+1:r*H , (r-1)*L+1:r*L) = D_mat(1:H , (r-1)*L+1:r*L);
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
function [X]=mat2tens(X_mat,mode,size_vec)
%MAT2TENS Tensorization of a matrix (reciprocal of tens2mat)     
if mode<1 || mode >3
    error ('Input argument mode must be a 1 , 2 or 3')
end
I=size_vec(1);
J=size_vec(2);
K=size_vec(3);
if mode==1
    X=permute(reshape(X_mat,K,I,J),[2 3 1]);
elseif mode==2
    X=reshape(X_mat,I,J,K);
elseif mode==3
    X=permute(reshape(X_mat,J,K,I),[3 1 2]);
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
function [A,B,C,D]= bcdLMN_init(X,L_vec,M_vec,N_vec)
%BCDLMN_INIT Random initialization for the BCD-(L,M,N)
[I,J,K]=size(X); 
R=length(L_vec);
A=cell(1,R);
B=cell(1,R);
C=cell(1,R);
D=cell(1,R);
%  generate random matrices
if isreal(X)
    for r=1:R
    A{r}=randn(I,L_vec(r));
    B{r}=randn(J,M_vec(r));
    C{r}=randn(K,N_vec(r));
    end
else
    for r=1:R
     A{r}=randn(I,L_vec(r))+j*randn(I,L_vec(r));
     B{r}=randn(J,M_vec(r))+j*randn(J,M_vec(r));
     C{r}=randn(K,N_vec(r))+j*randn(K,N_vec(r));
    end
end

% Orthogonalize matrices or submatrices
 if (I>=sum(L_vec))   % A is tall or square
    A=mat2cell(orth(cell2mat(A)),I,L_vec);
 else 
     for r=1:R;
        A{r}=orth(A{r});
     end
 end
 if (J>=sum(M_vec))   % B is tall or square
     B=mat2cell(orth(cell2mat(B)),J,M_vec);
 else
      for r=1:R;
          B{r}=orth(B{r});
      end
 end
 if (K>=sum(N_vec))   % C is tall or square
     C=mat2cell(orth(cell2mat(C)),K,N_vec);
 else
      for r=1:R;
          C{r}=orth(C{r});
      end
 end 
 X3=tens2mat(X,3);  % X3 is KJxI (idem)
 X_vec=X3(:);       % X_vec is IKJx1
 [D]=update_D(A,B,C,X_vec,L_vec,M_vec,N_vec);
end

%*******************************************************************************
function [A,B,C] = bcdLMN_lsearch(A,B,C,dA,dB,dC,D,X,lsearch,it,L_vec,M_vec,N_vec)     
%BCDLMN_lsearch Line Search for BCD-(L,M,N)
%   New loading A, B, C are computed
%   from their previous values and from the search directions dA, dB, dC,
%   as follows:
%  
%      A <- A + mu * dA
%      B <- B + mu * dB
%      C <- C + mu * dC
%
%   Line Search for bcdLMN can for instance be used in gradient-based
%   algorithms or with Alternating Least Squares (ALS), in order to speed 
%   up convergence.
%   The core tensors Dr are kept constant, although they could in principle
%   also be interpolated.
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
%      - A,B,C: estimates in cell format of the loading matrices A,B,C at it-1
%      - D is the cell that holds the set of R current core tensors Dr,     
%      - dA,dB,dC: search directions in cell format
%      - X: JIxK matrix unfolding of the observed tensor X 
%      - lsearch: = 'none', 'lsh','lsb', 'elsr' or 'elsc'
%      - it is the iteration step number   
%      - L_vec: the way A is partitioned
%      - M_vec: the way B is partitioned
%      - N_vec: the way C is partitioned
%
%   OUTPUTS: 
%      - Updated loadings A,B,C in cell format

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

 % Convert inputs from cell format to matrix format
A=cell2mat(A);dA=cell2mat(dA);
B=cell2mat(B);dB=cell2mat(dB);
C=cell2mat(C);dC=cell2mat(dC);
R=length(L_vec);
Dm=blockdiag(tens2mat_block(D,2),R);  % size RMLxRN (block diagonal matrix)
        

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
        phi=norm(X-kr_part(B,A,M_vec,L_vec)*Dm*C.','fro');
        % Compute phi with mu=1.25
        mu=1.25;
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
        phim=norm(X-kr_part(Bm,Am,M_vec,L_vec)*Dm*Cm.','fro');
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
        phi=norm(X-kr_part(B,A,M_vec,L_vec)*Dm*C.','fro');
       % Compute phi with mu=it^(1/3)
        mu=it^(1/3);
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
        phim=norm(X-kr_part(Bm,Am,M_vec,L_vec)*Dm*Cm.','fro');
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
        KdBdA=kr_part(dB,dA,M_vec,L_vec);
        KBdA=kr_part(B,dA,M_vec,L_vec);
        KdBA=kr_part(dB,A,M_vec,L_vec);
        KBA=kr_part(B,A,M_vec,L_vec);
        DmdC=Dm*dC.';
        DmC=Dm*C.';
        Mat3=KdBdA*DmdC;
        Mat2=KdBdA*DmC + (KBdA+KdBA)*DmdC;
        Mat1=KBA*DmdC + (KBdA+KdBA)*DmC;
        Mat0=KBA*DmC-X;
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
        KdBdA=kr_part(dB,dA,M_vec,L_vec);
        KBdA=kr_part(B,dA,M_vec,L_vec);
        KdBA=kr_part(dB,A,M_vec,L_vec);
        KBA=kr_part(B,A,M_vec,L_vec);
        DmdC=Dm*dC.';
        DmC=Dm*C.';
        Mat3=KdBdA*DmdC;
        Mat2=KdBdA*DmC + (KBdA+KdBA)*DmdC;
        Mat1=KBA*DmdC + (KBdA+KdBA)*DmC;
        Mat0=KBA*DmC-X;
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

% Convert to cell format
A=mat2cell(A,size(A,1),L_vec);
B=mat2cell(B,size(B,1),M_vec);
C=mat2cell(C,size(C,1),N_vec);
end
