% demo script to test bcdLM_alsls algorithm.
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
    I=14;                      % Dimensions  of the tensor
    J=9;
    K=30;
    R=5;                       % number of components
    L=3;                       % rank of the R component matrices Ar
    M=2;                       % rank of the R component matrices Br
    L_vec=L*ones(1,R);
    M_vec=M*ones(1,R);
    SNR=40;                    % SNR [dB], choose SNR=inf for a noise-free model
    power_max=3;               % Each of the R terms is normalized and then multiplied by a value between 1 and power_max. 
                               % The highest this value, the more difficult the problem is (swamps are more likely)
    power_vec=linspace(1,power_max,R);   % Holds the set of R amplitude values
       
    %--- Algorithm parameters
    comp='on';          % ='on' or ='off' to perform or not dimensionality reduction 
    Tol1=1e-6;          % Tolerance 
    MaxIt1=1000;        % Max number of iterations
    Tol2=1e-4;          % tolerance in refinement stage (after decompression)
    MaxIt2=50;         % Max number of iterations in refinement stage
    Ninit=5;            % Number of initializations

    
%***************************************************
%---- Build Loading matrices and observed tensor-----
%**************************************************** 
A_cell=cell(1,R);
B_cell=cell(1,R);
C_cell=cell(1,R);
if strcmp(data_type,'real')==1
    for r=1:R
    A_cell{r}=randn(I,L);B_cell{r}=randn(J,M);C_cell{r}=randn(L,M,K);
    end
elseif strcmp(data_type,'complex')==1
    for r=1:R
    A_cell{r}=randn(I,L)+j*randn(I,L);B_cell{r}=randn(J,M)+j*randn(J,M);C_cell{r}=randn(L,M,K)+j*randn(L,M,K);
    end
end
% Create observed tensor that follows the BCD-(L,M,.) model
Xr=cell(1,R);
X=zeros(I,J,K);
for r=1:R
    XX=zeros(I,J,K);
    for k=1:K
        XX(:,:,k)=A_cell{r}*C_cell{r}(:,:,k)*B_cell{r}.';
    end
    Xr{r}=power_vec(r)*XX/norm(reshape(XX,I*J,K),'fro');   %--> normalize and weight the rth tensor
    X=X+Xr{r};   % --> Build the observed tensor
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
% Compute the BCD-(L,M,.) for several initializations
%-------------------------------------------------------------------------------
for ninit=1:Ninit
    disp('-----------------------------------------------------')
    disp(['Initialization ',num2str(ninit)])
    disp('-----------------------------------------------------')
    % generate starting point
      [A_init,B_init,C_init]=bcdLM_init(X,L_vec,M_vec);     
         
    %-----------------------------------------------------
    %  COMPUTE THE DECOMPOSITION for this initialization
    %-----------------------------------------------------
    disp('Algorithm 1: als without line search')
    [A_est,B_est,C_est,phi,it1,it2,phi_als]=bcdLM_alsls(X,R,L,M,'none',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A_cell),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B_cell),M_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp('  ')
    
    disp('Algorithm 2: als with line search proposed by Harshman')
    [A_est,B_est,C_est,phi,it1,it2,phi_lsh]=bcdLM_alsls(X,R,L,M,'lsh',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A_cell),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B_cell),M_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp('  ')

    disp('Algorithm 3: als with line search proposed by Bro')
    [A_est,B_est,C_est,phi,it1,it2,phi_lsb]=bcdLM_alsls(X,R,L,M,'lsb',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A_cell),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B_cell),M_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp('  ')

    disp('Algorithm 4: als with exact line search and real step')
    [A_est,B_est,C_est,phi,it1,it2,phi_elsr]=bcdLM_alsls(X,R,L,M,'elsr',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A_cell),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B_cell),M_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
    disp('  ')

    disp('Algorithm 5: als with exact line search and complex step')
    [A_est,B_est,C_est,phi,it1,it2,phi_elsc]=bcdLM_alsls(X,R,L,M,'elsc',comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
    [errA]=solve_blockperm(cell2mat(A_est),cell2mat(A_cell),L_vec);
    [errB]=solve_blockperm(cell2mat(B_est),cell2mat(B_cell),M_vec);
    disp(['Final phi :  ',num2str(phi)])
    disp(['It1       :  ',num2str(it1)])
    disp(['It2       :  ',num2str(it2)])
    disp(['Error on A:  ',num2str(errA)])
    disp(['Error on B:  ',num2str(errB)])
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