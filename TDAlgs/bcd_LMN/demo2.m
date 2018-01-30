% demo script to test bcdLMN_alsls algorithm.
% The purpose of this demo is to show how to exploit several initializations.
% The function bcdLMN_alsls has 2 kinds of optional input arguments to deal with initialization:
% - Ninit : the algorithm will generate Ninit different starting points and for
%   each starting point, it will stop when the stop criterion is satisfied 
%  (input parameters Tol1 and MaxIt1 to select tolerance and max number of iterations).
%  Output parameters are the ones associated to the best starting point.
% - A,B,C: if one or more matrices are provided, the latter are used as a starting point
%   and this is the only starting point that will be tried, no matter the value of Ninit
%  (it is for instance useful when one already knows estimates close to the solution).

% When one uses e.g. Ninit=10, it can take time before the algorithm stops because
% each starting point will be tried until the convergence criterion is satisfied.
% Thus, it can be useful in practice to try Ninit=10 starting points for a limited
% number of iterations, say MaxItInit (less than MaxIt1) and then go on until convergence
% with only the best initialization obtained among the Ninit tried.
% This demo show how to to this.

clear all
close all
clc

%**********************************************
%--- Choose PARAMETERS of the DEMO
%**********************************************
    %--- Data parameters ----
    data_type='complex';       % choose 'real' or 'complex' to select the kind of data to generate
    I=12;                      % Dimensions  of the tensor
    J=11;
    K=20;
    R=3;                       % number of components
    L=3;                       % rank of the R component matrices Ar
    M=2;                       % rank of the R component matrices Br
    N=2;                       % rank of the R component matrices Cr
    L_vec=L*ones(1,R);
    M_vec=M*ones(1,R);
    N_vec=N*ones(1,R);
    SNR=inf;                

    
%***************************************************
%---- Build Loading matrices and observed tensor-----
%**************************************************** 
% Generate loading matrices
A=cell(1,R);
B=cell(1,R);
C=cell(1,R);
D=cell(1,R);
if strcmp(data_type,'real')==1
    for r=1:R
    A{r}=randn(I,L_vec(r));
    B{r}=randn(J,M_vec(r));
    C{r}=randn(K,N_vec(r));
    D{r}=randn(L_vec(r),M_vec(r),N_vec(r));
    end
elseif strcmp(data_type,'complex')==1
    for r=1:R
     A{r}=randn(I,L_vec(r))+j*randn(I,L_vec(r));
     B{r}=randn(J,M_vec(r))+j*randn(J,M_vec(r));
     C{r}=randn(K,N_vec(r))+j*randn(K,N_vec(r));
     D{r}=randn(L_vec(r),M_vec(r),N_vec(r))+j*randn(L_vec(r),M_vec(r),N_vec(r));
    end
end
% Create observed tensor that follows the BCD-(L,M,N) model
X=zeros(I,J,K);
for r=1:R
    X=X+tmprod(tmprod(tmprod(D{r},A{r},1),B{r},2),C{r},3);
end

% Add noise
if strcmp(data_type,'real')==1
   Noise_tens=randn(I,J,K);
else
   Noise_tens=randn(I,J,K)+j*randn(I,J,K);
end
sigma=(10^(-SNR/20))*(norm(reshape(X,J*I,K),'fro')/norm(reshape(Noise_tens,J*I,K),'fro'));
X=X+sigma*Noise_tens;


%-------------------------------------------------------------------------------
% Compute the BCD-(L,M,N) with a limited number of iterations for each starting point
% and then use the best starting point to go on until convergence
%-------------------------------------------------------------------------------
Ninit=50;
MaxItInit=10;
% Perform at most MaxItInit iterations for each starting point
[A_init,B_init,C_init,D_init,phi,itInit]=bcdLMN_alsls(X,R,L,M,N,[],[],[],MaxItInit,[],[],Ninit);
% Then go on with the best starting point
[A_est,B_est,C_est,D_init,phi,it1]=bcdLMN_alsls(X,R,L,M,N,[],[],[],[],[],[],[],A_init,B_init,C_init);
% compute relative errors
[err_A]=solve_blockperm(cell2mat(A_est),cell2mat(A),L_vec);
[err_B]=solve_blockperm(cell2mat(B_est),cell2mat(B),M_vec);
[err_C]=solve_blockperm(cell2mat(C_est),cell2mat(C),M_vec);
disp(['err on A:         ',num2str(err_A)])
disp(['err on B:         ',num2str(err_B)])
disp(['err on C:         ',num2str(err_B)])
disp(['phi:              ',num2str(phi)])
disp(['itInit:           ',num2str(itInit)])
disp(['it1:              ',num2str(it1)])
