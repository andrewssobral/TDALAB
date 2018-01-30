% demo to test  bcdLrLr1_alsls algorithm.

clear all
close all
clc

%**********************************************
%--- Choose PARAMETERS of the DEMO
%**********************************************
    %--- Data parameters ----
    data_type='complex';       % choose 'real' or 'complex' to select the kind of data to generate
    I=16;                      % Dimensions  of the tensor
    J=21;
    K=10;
    L_vec=[2 1 3];
    R=length(L_vec);
    N=sum(L_vec);
    SNR=70;                    % SNR [dB], choose SNR=inf for a noise-free model
    power_max=3;               % Each of the R terms is normalized and then multiplied by a value between 1 and power_max. 
                               % The highest this value, the more difficult the problem is (swamps are more likely)
    power_vec=linspace(1,power_max,R);   % Holds the set of R amplitude values
       
    %--- Algorithm parameters
    lsearch='elsr';
    comp='on';          % ='on' or ='off' to perform or not dimensionality reduction 
    Tol1=1e-8;          % Tolerance 
    MaxIt1=5000;        % Max number of iterations
    Tol2=1e-6;          % tolerance in refinement stage (after decompression)
    MaxIt2=500;         % Max number of iterations in refinement stage
    Ninit=5;            % Number of initializations

    
%***************************************************
%---- Build Loading matrices and observed tensor-----
%**************************************************** 
IndL=[0, cumsum(L_vec)];   % Pattern of the partition
if strcmp(data_type,'real')==1
    A=randn(I,N);B=randn(J,N);C=randn(K,R);
elseif strcmp(data_type,'complex')==1   
    A=randn(I,N)+j*randn(I,N);B=randn(J,N)+j*randn(J,N);C=randn(K,R)+j*randn(K,R);
end

% Create observed tensor that follows the BCD-(L,L,1) model
X_mat=zeros(I*J,K);
for r=1:R
    Xr=reshape(A(:,IndL(r)+1:IndL(r+1))*B(:,IndL(r)+1:IndL(r+1)).',I*J,1)*C(:,r).';
    X_mat=X_mat + power_vec(r)*Xr/norm(Xr,'fro');
end
X=reshape(X_mat,I,J,K);

%---- Possibly Add noise---------
if strcmp(data_type,'real')==1
   Noise_tens=randn(I,J,K);
elseif strcmp(data_type,'complex')==1
   Noise_tens=randn(I,J,K)+j*randn(I,J,K);
end
sigma=(10^(-SNR/20))*(norm(reshape(X,J*I,K),'fro')/norm(reshape(Noise_tens,J*I,K),'fro'));
X=X+sigma*Noise_tens;


%-------------------------------------------------------------------------------
% Compute the BCD-(Lr,Lr,1)
%-------------------------------------------------------------------------------
tic
[A_est,B_est,C_est,phi,it1,it2,phi_vec]=bcdLrLr1_alsls(X,L_vec,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit);
% To use all default parameters:
% [A_est,B_est,C_est,phi,it1,it2,phi_vec]=bcdLrLr1_alsls(X,L_vec);
time=toc;
% compute relative errors
[err_A]=solve_blockperm(A_est,A,L_vec);
[err_B]=solve_blockperm(B_est,B,L_vec);
[err_C]=solve_blockperm(C_est,C,ones(1,R));
disp(['err on A:         ',num2str(err_A)])
disp(['err on B:         ',num2str(err_B)])
disp(['err on C:         ',num2str(err_C)])
disp(['phi:              ',num2str(phi)])
disp(['it1:              ',num2str(it1)])
disp(['it2:              ',num2str(it2)])


loglog(phi_vec)
title('Evolution of phi for the best initialization')
xlabel('Iterations')
ylabel('phi')