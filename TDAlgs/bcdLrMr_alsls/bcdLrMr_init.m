function [A,B,C]= bcdLM_init(X,L_vec,M_vec)
%BCDLM_INIT Random initialization for the BCD-(L,M,.)
[I,J,K]=size(X); 
R=length(L_vec);
A=cell(1,R);
B=cell(1,R);
C=cell(1,R);
%  generate random matrices
if isreal(X)
    for r=1:R
    A{r}=randn(I,L_vec(r));
    B{r}=randn(J,M_vec(r));
    C{r}=randn(L_vec(r),M_vec(r),K);
    end
else
    for r=1:R
     A{r}=randn(I,L_vec(r))+j*randn(I,L_vec(r));
     B{r}=randn(J,M_vec(r))+j*randn(J,M_vec(r));
     C{r}=randn(L_vec(r),M_vec(r),K)+j*randn(L_vec(r),M_vec(r),K);
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
end
