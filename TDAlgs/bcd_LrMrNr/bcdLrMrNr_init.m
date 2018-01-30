function [A,B,C,D]= bcdLrMrNr_init(X,L_vec,M_vec,N_vec)
%BCDLMN_INIT Random initialization for the BCD-(Lr,Mr,Nr)
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

