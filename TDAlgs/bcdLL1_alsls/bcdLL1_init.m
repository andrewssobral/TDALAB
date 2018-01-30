function [A,B,C]= bcdLL1_init(X,R,L,init_type)
%BCDLL1_INIT Initialization of the loading matrices for the BCD-(L,L,1)
% INPUTS: 
% - X: 3rd order tensor of size (IxJxK)
% - R: Number of rank-(L,L,1) components
% - init_type='gevd' to initialize by generalized EVD if possible,
%   otherwise randomly,
% - init_type='random' to enforce random initialization.
% OUTPUTS:
% A(IxRL) B(JxRL) and C(KxR): Loading Matrices
[I,J,K]=size(X);
RL=R*L;
% Check if I and J are both greater than RL
if min(I,J)<RL
    init_type='random';    %  use random init
end

if strcmp(init_type,'gevd')==1
    [A,B,C] = bcdLL1_gevd(X,R,L);
elseif strcmp(init_type,'random')==1
   if isreal(X);
    A=rand(I,RL);B=rand(J,RL);C=rand(K,R);
   else
    A=rand(I,RL)+j*rand(I,RL);B=rand(J,RL)+j*rand(J,RL);C=rand(K,R)+j*rand(K,R);
   end
   % Orthogonalize matrices or submatrices
   if (I>=RL);A=orth(A);
   else ; for r=1:R; A(:,(r-1)*L+1:r*L)=orth(A(:,(r-1)*L+1:r*L)); end; end
   if (J>=RL);B=orth(B);
   else ; for r=1:R; B(:,(r-1)*L+1:r*L)=orth(B(:,(r-1)*L+1:r*L)); end; end
   if K>=R; C=orth(C); else; C=orth(C');end
end
end

%*******************************************************************************
function [A,B,C] = bcdLL1_gevd(X,R,L)
% BCD_LL1_GEVD BCD-(L,L,1) by generalized EVD of two slices.
[I,J,K]=size(X);
RL=R*L;
% Check if I and J are both greater than RL
if min(I,J)<RL
       error('bcdLL1_gevd:InvalidDimensions',['The input tensor ' ...
          'must have its first two dimensions greater than or equal to RL']);
end
% Truncated MLSVD
[U1,U2,U3,X_core] = mlsvd3(X,[RL,RL,2]);
% Exploit the two slices of X_core to find A (LRxLR) and B(LRxLR)
% This is a generalized eigenvalue problem
X1=X_core(:,:,1);
X2=X_core(:,:,2); 
X12=X1/X2;
% Eigenvalue Decomposition
[A S]=eig(X12); 
s=diag(S);
[a b]=sort(abs(s),'descend');
s=s(b(1:RL));
S = diag(s);
% In absence of noise, the eigenvalues consist
% of R subsets, each subset consisting of L equal eigenvalues, 
% to be associated with one of the R components

% Find A
A = A(:,b(1:RL));  % b is the permutation vector that associates the eigenvectors L by L

% If X is real but A is complex (because X12 has one or several pairs of conjugate eigenvalues)
% then transform the columns of A with a similarity matrix to get a real-valued one
if isreal(s)~=1 && isreal(X)==1  
     r=1;
     T=[0.5 -i*0.5;0.5 i*0.5];  % 2 by 2 block: similarity matrix
     while r<=RL
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
Ps=kron(eye(R),ones(L,1));
X2=tens2mat(X,2);
C=X2.'/(kr(B,A)*Ps).';
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
