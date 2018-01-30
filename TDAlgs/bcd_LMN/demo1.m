% demo script to test bcdLMN_alsls algorithm.
% The purpose of this demo is to observe the different convergence curves for 
% different line search method and different random initializations.
clear all
close all
clc

%**********************************************
%--- Choose PARAMETERS of the DEMO
%**********************************************
    %--- Data parameters ----
    data_type='complex';       % choose 'real' or 'complex' to select the kind of data to generate
    I=16;                      % Dimensions  of the tensor
    J=17;
    K=18;
    R=2;                       % number of components
    L=3;                       % rank of the R component matrices Ar
    M=2;                       % rank of the R component matrices Br
    N=2;                       % rank of the R component matrices Cr
    L_vec=L*ones(1,R);
    M_vec=M*ones(1,R);
    N_vec=N*ones(1,R);
    SNR=inf;                    % SNR [dB], choose SNR=inf for a noise-free model
       
    %--- Algorithm parameters
    comp='on';          % ='on' or ='off' to perform or not dimensionality reduction 
    Tol1=1e-6;          % Tolerance 
    MaxIt1=300;        % Max number of iterations
    Tol2=1e-4;          % tolerance in refinement stage (after decompression)
    MaxIt2=50;         % Max number of iterations in refinement stage
    Ninit=5;            % Number of initializations

    
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
% Compute the BCD-(L,M,N) for several initializations
%-------------------------------------------------------------------------------
for ninit=1:Ninit
    disp('-----------------------------------------------------')
    disp(['Initialization ',num2str(ninit)])
    disp('-----------------------------------------------------')
    % generate starting point
    [A_init,B_init,C_init,D_init]= bcdLMN_init(X,L_vec,M_vec,N_vec);
    %-----------------------------------------------------
    %  COMPUTE THE DECOMPOSITION for this initialization
    %-----------------------------------------------------
    disp('Algorithm 1: als without line search')
    [A_est,B_est,C_est,D_est,phi,it1,it2,phi_als]=bcdLMN_alsls(X,R,L,M,N,'none',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B),M_vec);
    [errC]=solve_blockperm(cell2mat(C_est),cell2mat(C),N_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp(['Error on C:  ',num2str(errC)])
    disp('  ')
    
    disp('Algorithm 2: als with line search proposed by Harshman')
    [A_est,B_est,C_est,D_est,phi,it1,it2,phi_lsh]=bcdLMN_alsls(X,R,L,M,N,'lsh',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B),M_vec);
    [errC]=solve_blockperm(cell2mat(C_est),cell2mat(C),N_vec);   
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp(['Error on C:  ',num2str(errC)])    
    disp('  ')

    disp('Algorithm 3: als with line search proposed by Bro')
    [A_est,B_est,C_est,D_est,phi,it1,it2,phi_lsb]=bcdLMN_alsls(X,R,L,M,N,'lsb',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B),M_vec);
    [errC]=solve_blockperm(cell2mat(C_est),cell2mat(C),N_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp(['Error on C:  ',num2str(errC)])
    disp('  ')

    disp('Algorithm 4: als with exact line search and real step')
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B),M_vec);
    [errC]=solve_blockperm(cell2mat(C_est),cell2mat(C),N_vec);
    [A_est,B_est,C_est,D_est,phi,it1,it2,phi_elsr]=bcdLMN_alsls(X,R,L,M,N,'elsr',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp(['Error on C:  ',num2str(errC)])
    disp('  ')

    disp('Algorithm 5: als with exact line search and complex step')
    [A_est,B_est,C_est,D_est,phi,it1,it2,phi_elsc]=bcdLMN_alsls(X,R,L,M,N,'elsc',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B),M_vec);
    [errC]=solve_blockperm(cell2mat(C_est),cell2mat(C),N_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp(['Error on C:  ',num2str(errC)])
    disp('  ')

    
    figure
    loglog(phi_als,'b','LineWidth',2);hold on;
    loglog(phi_lsh,'y','LineWidth',2);
    loglog(phi_lsb,'r','LineWidth',2);
    loglog(phi_elsr,'m','LineWidth',2);
    loglog(phi_elsc,'c','LineWidth',2);
    xlabel('Iterations')
    ylabel('phi')
    title(['Evolution of phi for initialization ',num2str(ninit)])
    legend('als','als+lsh','als+lsb','als+elsr','als+elsc')
end