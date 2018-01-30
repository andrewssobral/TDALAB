function [err,Ae2,perm_block,perm_col,P,D]=solve_blockperm(Ae,A,L_vec,method)
%SOLVE_BLOCKPERM find block permutation and scaling for blocks of Se to match blocks of S
%   [err,Ae2,perm_block,perm_col,P,D]=solve_blockperm(Ae,A,L_vec,method)
% Let A=[A1,A2,...,AR] be a partitioned matrix, where the R blocks Ar (IxLr), 
% r=1,...,R, may have a different number of columns, the partition being given 
% by L_vec=[L1,L2,...,LR].
% A may be fat or tall, but each block Ar is tall and full column rank. 
%
% Let Ae be an estimate of A which is linked to A as follows,
%               Ae = A * D * P  +  Noise
%               A = Ae * inv(P) * inv(D) + Noise
% where 
% - P is a block-wise permutation matrix that permutes the R blocks of Ae
% - D a block-diagonal matrix with R blocks, each of size LrxLr, on the diagonal
%
% This function does the following: 
% - Estimation of P 
% - Estimation of D 
% - Builds  Ae2 = Ae * inv(P) * inv(D), such that Ae2=A in case of perfect 
%   estimation (Noise=0).
% - Computation of the relative error err = norm(Ae2-A,'fro')/norm(A,'fro')
%-------------------------------------------------------------------------------
% Inputs: 
% - Ae (I by N)   is an estimate of A
% - A  (I by N)   is the true matrix
% - L_vec, with N=sum(L_vec) holds the way A is partitioned
%
% Optional input
% - method     : ='subspace' (default)  to find the permutation by measuring 
%                   angles between subspaces
%               ='LS'   to use a least squares approach, by computing pinv(Ae)*A
%-------------------------------------------------------------------------------
% Outputs:
% - err     norm of the residual
% - Ae2     = Ae * inv(P) * inv(D)
% - perm_block    : the way the R blocks of A  are permuted  
%                   (vector holding a permutation of 1:R)
% - perm_col      : the way the N columns of A are permuted 
% - P             : N by N permutation matrix, built according to the 
%                   permutation pattern in perm_col and perm_block
% - D             : block-diagonal matrix
%-------------------------------------------------------------------------------
% Example:
%     clear all
%     % select problem
%     SNR=10;
%     L_vec=[3 3 3];
%     I=10;
%     method='LS';
%     R=length(L_vec);
%     N=sum(L_vec);
%     CumL=[0,cumsum(L_vec)];
%     % Generate A
%     A=randn(I,N)+j*randn(I,N);
%     % Generate block diagonal matrix D
%     D=zeros(N,N);
%     for r=1:R;
%         D(CumL(r)+1:CumL(r+1),CumL(r)+1:CumL(r+1))=randn(L_vec(r),L_vec(r))+j*randn(L_vec(r),L_vec(r));
%     end
%     % generate permutation
%     perm_block=randperm(R);    % Permutation matrix L_vec(P) is the way Ae is partitioned
%     perm_col=mat2cell(1:N,1,L_vec);
%     perm_col=cell2mat(perm_col(:,perm_block));
%     IN=eye(N);
%     P=IN(:,perm_col);
%     % build Ae = A * D * P
%     Ae=A * D * P;        
%     % Add noise
%     Noise=randn(size(Ae))+j*randn(size(Ae));
%     sigma=(10^(-SNR/20))*(norm(Ae(:),'fro')/norm(Noise(:),'fro'));
%     Ae=Ae+sigma*Noise;    
%     % Estimate permutation and compute residual Mean Square Error
%     [err,Ae2, perm_block_est, perm_col_est, P_est,D_est]=solve_blockperm(Ae,A,L_vec,method);
%     err
%     norm(D-D_est)
%     nnz(perm_col-perm_col_est)
%     nnz(perm_block-perm_block_est)
%     nnz(P-P_est)
%-----------------------------------------------------------------------------------------------------    
% Store dimensions
R=length(L_vec);
I=size(Ae,1);
N=size(Ae,2);
if N ~= sum(L_vec)
    error('The number of columns of Q and Qe has to be equal to sum(L_vec)')
end
if I~=size(A,1) || N~=size(A,2)
    error('A and Ae must have the same dimensions')
end
if exist('method')~=1 || isempty(method);
    method='subspace';
end
CumL=[0,cumsum(L_vec)];


%-----------------------------------------------------------------------------
% STEP 1: find the permutation matrix and build Ae2=Ae*inv(P)=A*D
%-----------------------------------------------------------------------------
[Ae2, perm_block, perm_col, P]=perm_bdiag(Ae,A,L_vec,method);


%-----------------------------------------------------------------------------
% STEP 2: find the block diagonal matrix D and build Ae2=Ae2*inv(D)
%-----------------------------------------------------------------------------
D=zeros(N,N);
for r=1:R
    invDr=pinv(Ae2(:,CumL(r)+1:CumL(r+1)))*A(:,CumL(r)+1:CumL(r+1));
    Ae2(:,CumL(r)+1:CumL(r+1))=Ae2(:,CumL(r)+1:CumL(r+1))*invDr;
    D(CumL(r)+1:CumL(r+1),CumL(r)+1:CumL(r+1))=pinv(invDr);
end

%-----------------------------------------------------------------------------
% STEP 3: Compute relative Mean Square Error
%-----------------------------------------------------------------------------
err=norm(A-Ae2,'fro')/norm(A,'fro');

% version with normalization of each column of A and Ae2
%err=norm(A*diag(1./sqrt(diag(A'*A)))-Ae2*diag(1./sqrt(diag(Ae2'*Ae2))),'fro')...
%/norm(A*diag(1./sqrt(diag(A'*A))),'fro');


%*******************************************************************************
function [Ae2, perm_block, perm_col, P]=perm_bdiag(Ae,A,L_vec,method)
%PERM_BDIAG Find block-wise permutation between two matrices.
% Blocks may be of different size and matrices may be fat or tall.
%
% Let A=[A1,A2,...,AR] be a partitioned matrix, where the R blocks Ar (IxLr), 
% r=1,...,R, may have a different number of columns, the partition being given 
% by L_vec=[L1,L2,...,LR].
% A may be fat or tall, but each block Ar is tall and full column rank. 
%
% Let Ae be an estimate of A which is linked to A as follows,
%               Ae = A * D * P  +  Noise
%               A = Ae * inv(P) * inv(D) + Noise
% where 
% - P is a block-wise permutation matrix that permutes the R blocks of A
% - D a block-diagonal matrix with R blocks, each of size LrxLr, on the diagonal
%
% This function does the following: 
% - Estimation of P
% - Builds  Ae2 = Ae * inv(P), such that Ae2 matches A, only up to a block 
% diagonal matrix (inv(D)), in case where Noise=0
%-------------------------------------------------------------------------------
% Inputs: 
% - Ae (I by N)   is an estimate of A
% - A  (I by N)   is the true matrix
% - L_vec, with N=sum(L_vec) holds the way A is partitioned
%
% Optional input
% - method     : ='subspace' (default)  to find the permutation by measuring 
%                 angles between subspaces
%               ='LS'   to use a least squares approach, by computing pinv(Ae)*A
%-------------------------------------------------------------------------------
% Outputs:
% - Ae2           : Ae2 = Ae * inv(P), is the block-wise permuted version of Ae
% - perm_block    : the way the R blocks of A are permuted
% - perm_col      : the way the RL columns of A are permuted,
% - P             : permutation matrix, built according to perm_col and perm_block
%-------------------------------------------------------------------------------
% Method: 
%  - 'LS' : only used if all values of L_vec are identical, i.e., all subblocks 
%           Ar have the same size, e.g., L_vec=[2 2 2 2]
%           * CASE 1.1: A is tall or square (and assumed full column rank)
%     ---> then find the permutation by examination of pinv(Ae)*A, which is a 
%           block-wise permuted block diagonal matrix
%           * CASE 1.2: A is fat (and assumed full row rank)
%     ---> then find the permutation by examination of pinv(KAe)*KA, where 
%          KA=A(kr)A and KAe=Ae(kr)Ae, where (kr) denotes the Khatri Rao product
%         (block-wise Kronecker product). It is assumed that KA and KAe are full
%         column rank, such that one can proceed as in case 1.1. 
%         So we must have I^2 > sum(L_vec.^2).
% - 'subspace' : works in all cases. A may be fat or tall, with blocks of the 
%   same size or not. We use angle between subspaces. For instance, if 
%   L_vec=[2 1 2 1 3], we first select the two columns of A with Lr=1 and for 
%   each one, we measure the angle between the selected column and all columns 
%   of Ae and select the columns of Ae with the minimum angle. Then we do the 
%   same for the two blocks Lr=2 and so on
%-------------------------------------------------------------------------------
% Example: uncomment the code below to test perm_bdiag
%  clear all
%     % select problem
%     pb='22';        % '11': A is tall with blocks of same size
%                     % '12': A is fat with blocks of same size
%                     % '21': A is tall with blocks of different size
%                     % '22': A is fat with blocks of different size
% 
%     if strcmp(pb,'11')==1
%         R=6;
%         L=4;
%         L_vec=L*ones(1,R);
%         N=sum(L_vec);
%         I=N+1;   % to make A tall
%     elseif strcmp(pb,'12')==1
%         R=4;
%         L=3;
%         L_vec=L*ones(1,R);
%         N=sum(L_vec);
%         I=N-L;    % to make A fat    
%     elseif strcmp(pb,'21')==1
%         L_vec=[4 2 3 1];
%         R=length(L_vec);
%         N=sum(L_vec);
%         I=N+1;   % to make A tall
%     elseif strcmp(pb,'22')==1
%         L_vec=[2 1 3 4 1 3];
%         R=length(L_vec);
%         N=sum(L_vec);
%         I=N-R;   % to make A fat
%     end
%     CumL=[0,cumsum(L_vec)];
%     A=randn(I,N)+j*randn(I,N);
%     % Generate block diagonal matrix D
%     D=zeros(N,N);
%     for r=1:R;
%         D(CumL(r)+1:CumL(r+1),CumL(r)+1:CumL(r+1))=randn(L_vec(r),L_vec(r))+j*randn(L_vec(r),L_vec(r));
%     end
%     % generate permutation
%     perm_block=randperm(R);    % Permutation matrix L_vec(P) is the way Ae is partitioned
%     perm_col=mat2cell(1:N,1,L_vec);
%     perm_col=cell2mat(perm_col(:,perm_block));
%     IN=eye(N);
%     P=IN(:,perm_col);
%     % build Ae = A * D * P
%     Ae=A * D * P;        
% 
%     % Call permutation finder
%     [Ae2, perm_block_est, perm_col_est, P_est]=perm_bdiag(Ae,A,L_vec);
%     norm(A*D-Ae2)
%     nnz(perm_col-perm_col_est)
%     nnz(perm_block-perm_block_est)
%     nnz(P-P_est)   
%-------------------------------------------------------------------------------
% Store dimensions
R=length(L_vec);
I=size(Ae,1);
N=size(Ae,2);
if N ~= sum(L_vec)
    error('The number of columns of A and Ae has to be equal to sum(L_vec)')
end
if I~=size(A,1) || N~=size(A,2)
    error('A and Ae must have the same dimensions')
end
if exist('method')~=1 || isempty(method);
    method='subspace';
end

if (mean(L_vec)~=max(L_vec))
    method='subspace';
end
    
%*******************************************************************************
% CASE 1 : method='LS', only use if all values of L_vec are identical
%*******************************************************************************
if  strcmp(method,'LS')==1
    % CASE 1.1:  A is tall full column rank
    if  (I>=N)  
        [Ae2, perm_block, perm_col, P]=permSb(Ae,A,R);
    
    % CASE 2.2:  A is fat full row rank
    elseif (I<N) 
        Aein=Ae;
        Ae = kr_part(Ae,Ae,L_vec,L_vec);
        A  = kr_part(A,A,L_vec,L_vec);
        [Ae2, perm_block]=permSb(Ae,A,R);
        perm_col=mat2cell(1:N,1,L_vec);
        perm_col=cell2mat(perm_col(:,perm_block));                                 
        IN=eye(N);
        P=IN(:,perm_col);
        Ae2=Aein*inv(P);
    end
    
    
%******************************************************************************
% CASE 3: not all values of L are identical and A  can be fat or tall
%******************************************************************************  
elseif  strcmp(method,'subspace')==1
    
        CumL=[0,cumsum(L_vec)];
        Acell=cell(1,R);
        for r=1:R
            Acell{r}=A(:,CumL(r)+1:CumL(r+1));
        end
        [Lv,pL]=sort(L_vec,'ascend');       % we will associate subspaces from the smallest to the highest dimension
        
        Acell=Acell(:,pL);
        Aein=Ae;
        Lremain=1:N;
        Lvecp=[];
        % Now associate the matrices in Acell with one matrix in Ae, by measuring angle between subspaces
        for r=1:R
            Ar=Acell{r};  % size I by Lv(r)
            % measure the angle between Ar and all possible matrices of Ae of size I by Lv(r)
            ang=zeros(1,size(Ae,2)-Lv(r)+1);
            for m=1:size(Ae,2)-Lv(r)+1
                ang(m)=subspace(Ar/norm(Ar),Ae(:,m:m+Lv(r)-1)/norm(Ae(:,m:m+Lv(r)-1)));
            end
            pm=find(ang==min(ang));
            Ae(:,pm:pm+Lv(r)-1)=[];
            Lvecp=[Lvecp Lremain(pm:pm+Lv(r)-1)];
            Lremain(pm:pm+Lv(r)-1)=[];
        end
        
        uu=mat2cell(1:N,1,L_vec);
        uu=cell2mat(uu(:,pL));
        
        
        CumLv=[0,cumsum(Lv)];
        Lperm=cell(1,R);
        for r=1:R
            Lperm{r}=Lvecp(CumLv(r)+1:CumLv(r+1));    
        end
        Firstc=zeros(1,R);
        for r=1:R
            Firstc(r)=Lperm{r}(1);    % gives the first column of each block of L_vec(P)
        end
        [ff,Ls]=sort(Firstc,'ascend'); 
        LP=Lv(:,Ls);    % is equal to  L_vec(perm_block) gives the way Ae is partitioned
        perm_block=pL(:,Ls);     % vector holding the way the R blocks of A have been permuted

        perm_col=mat2cell(1:N,1,L_vec);
        perm_col=cell2mat(perm_col(:,perm_block)); 
        
        % Compute Ae * inv(P)
        Ae2=mat2cell(Aein,I,LP);
        Ae2(:,perm_block)=Ae2;       
        Ae2=cell2mat(Ae2);
    
        % Build matrix P
        IN=eye(N);
        P=IN(:,perm_col);    
    
end 

%*************************************************************************************************************************
function [Ae2, perm_block, perm_col, P]=permSb(Ae,A,R)
% Find block-wise permutation between two matrices where all blocks have the same size and matrices are full column rank
%------------------------------------------------------------------------------------------------------------------------
% Consider a partitioned matrix A=[A1,A2,...,AR], where all R blocks have the same size IxL,
% with I>=(RL), i.e., A is tall and an estimate Ae of A which is linked to A as follows,
%               Ae = A * D * P  +  Noise
%               A = Ae * inv(P) * inv(D) + Noise
% where 
% - P is a block-wise permutation matrix that permutes the R blocks 
% - D a block-diagonal matrix with R blocks, each of size LxL, on the diagonal
%
% This function estimates P 
% and builds the new estimate Ae2 = Ae * inv(P), such that Ae2 matches A only up to a block diagonal matrix (in case where Noise=0)
%-------------------------------------------------------------------------------------------------------------------------------
% Inputs: 
% - Ae (I by RL)   is an estimate of A, Ae=[Ae1, Ae2, ..., AeR]
% - A  (I by RL)   is the true matrix, A=[A1, A2, ...., AR]
% - R              is the number of IxL blocks in A and Ae
%------------------------------------------------------------------------------------------------------
% Outputs:
% - Ae2           : Ae2 = Ae * inv(P), is the block-wise permuted version of Ae
% - perm_block    : the way the R blocks of Ae are permuted
% - perm_col      : the way the RL columns of Ae are permuted,
% - P             : permutation matrix, built according to perm_col and perm_block
%------------------------------------------------------------------------------------------------------
% Example
%     clear all
%     L=2;
%     R=4;
%     I=24;
%     L_vec=repmat(L,1,R);
%     CumL=[0,cumsum(L_vec)];
%     N=sum(L_vec);
%     A=randn(I,N)+j*randn(I,N);
%     D=zeros(N,N);
%     for r=1:R;
%         D(CumL(r)+1:CumL(r+1),CumL(r)+1:CumL(r+1))=randn(L_vec(r),L_vec(r))+j*randn(L_vec(r),L_vec(r));
%     end
%     
%     perm_block=randperm(R);                     % vector holding the way the blocks of A*D will be permuted
%     Ind=reshape(1:R*L,L,R);
%     perm_col=reshape(Ind(:,perm_block),1,R*L);  % vector holding the way the columns of of A*D will be permuted
%     IN=eye(N);
%     P=IN(:,perm_col);
%     % build Ae = A * D * P  or  A = Ae*inv(P)*inv(D)
%     Ae=A * D * P;            
%     % estimate P
%     [Ae2, perm_block_est,perm_col_est,P_est]=permSb(Ae,A,R);   % Ae2=Ae*P
%     norm(D-pinv(A)*Ae2)
%     nnz(P-P_est)
%     nnz(perm_block-perm_block_est)
%     nnz(perm_col_est-perm_col)
%---------------------------------------------------------------------------------------------------------------------------

[I,N]=size(A);
L=N/R;
if mod(N,R)~=0
    error('error: wrong input argument R ')
end

V=pinv(Ae)*A;    % =P*D
V=diag(1./sqrt(diag(V*V')))*V;    % normalize each row
prod_block=abs(V);                                % the (LxL) blocks that have the highest norm

prod_init=zeros(R,R);
for r=1:R
    for k=1:R
        prod_init(r,k)=norm(prod_block((r-1)*L+1:r*L,(k-1)*L+1:k*L),'fro');   % substitute each LxL block by its frobenius norm 
    end
end
      
% prod_init is a permuted diagonal matrix so now the task is to find the position of the elements of this matrix with
% highest magnitude, which corresponds to the position of the L by L blocks.
prod=prod_init;
ordre1=1:R;
ordre2=1:R;
perm_opt1=ordre1;
perm_opt2=ordre2;
for p=1:R
    [row,col]=find(prod==max(max(prod)));
    row=row(1);col=col(1);
    perm_opt1(p)=ordre1(col);
    perm_opt2(p)=ordre2(row);
    ordre1(col)=[];
    ordre2(row)=[];
    prod(row,:)=[];  
    prod(:,col)=[];
end
[asce perm2]=sort(perm_opt2,'ascend');
perm_block=perm_opt1(perm2);
Ind=reshape(1:R*L,L,R);
perm_col=reshape(Ind(:,perm_block),1,R*L);

% Build matrix P
IN=eye(N);
P=IN(:,perm_col);

% Build Ae2=Ae*inv(P);
Ae2=Ae*inv(P);

%*************************************************************************************************************************
function Mat = kr_part(B,C,partB,partC)
        %---------------------------------------------------------------------------------------------
        % This function performs block wise kronecker product between the partitioned matrices B and C
        % Example: if B=[B1 B2 B3] and C=[C1 C2 C3], then M=[kron(B1,C1) kron(B2,C2) kron(B3,C3)]
        %---------------------------------------------------------------------------------------------
        % INPUTS:  
        % - B: of size JxM 
        % - C: of size KxN
        % - partB : vector that holds the number of columns of each block B1,B2,...
        % - partC : vector that holds the number of columns of each block C1,C2,...
        %   Example: partB=[L1 L2] means that M=L1+L2
        %-----------------------------------------------------------------------------------------------
        % OUTPUT
        % Mat of size JKxsum(partB.*partC) that contains the kronecker product of the blocks
        %------------------------------------------------------------------------------------------------
        
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
        
%*********************************************************************************         
function C = fast_kron (A,B)

        %--------------------------------------------------------------------------
        % Compute the Kronecker Product of A and B, equivalent to kron.m function
        %--------------------------------------------------------------------------
        [I,L1]=size(A);
        [J,L2]=size(B);

        if (L1==1) & (L2==1)
             C=reshape(B*A.',I*J,1);
        elseif (L1==1) & (L2>1)
            Bt=B.'; 
            C=reshape(Bt(:)*A.',L2,I*J).';
        elseif (L2==1) & (L1>1)
             C=reshape(B*A(:).',I*J,L1);
        else
             C=reshape(permute(reshape(B(:)*A(:).',[J,L2,I,L1]),[1 3 2 4]),[I*J,L1*L2]);
        end       