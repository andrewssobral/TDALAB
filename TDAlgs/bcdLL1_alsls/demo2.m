% demo 2 to test  bcdLL1_alsls algorithm.
% The purpose of this demo is to show how to use all default parameters

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
    R=4;                       % number of components
    L=4;                       % rank of the R component matrices A and B
    L_vec=L*ones(1,R);
    SNR=70;                    % SNR [dB], choose SNR=inf for a noise-free model
    power_max=3;               % Each of the R terms is normalized and then multiplied by a value between 1 and power_max. 
                               % The highest this value, the more difficult the problem is (swamps are more likely)
    power_vec=linspace(1,power_max,R);   % Holds the set of R amplitude values
     
%***************************************************
%---- Build Loading matrices and observed tensor-----
%**************************************************** 
if strcmp(data_type,'real')==1
    A=randn(I,R*L);B=randn(J,R*L);C=randn(K,R);
elseif strcmp(data_type,'complex')==1   
    A=randn(I,R*L)+j*randn(I,R*L);B=randn(J,R*L)+j*randn(J,R*L);C=randn(K,R)+j*randn(K,R);
end

% Create observed tensor that follows the BCD-(L,L,1) model
X_mat=zeros(I*J,K);
for r=1:R
    Xr=reshape(A(:,(r-1)*L+1:r*L)*B(:,(r-1)*L+1:r*L).',I*J,1)*C(:,r).';
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
% Compute the BCD-(L,L,1) with default parameters
%-------------------------------------------------------------------------------
tic
[A_est,B_est,C_est]=bcdLL1_alsls(X,R,L);
time=toc;
% compute relative errors
[err_A]=solve_blockperm(A_est,A,L_vec);
[err_B]=solve_blockperm(B_est,B,L_vec);
[err_C]=solve_blockperm(C_est,C,ones(1,R));
disp(['err on A:         ',num2str(err_A)])
disp(['err on B:         ',num2str(err_B)])
disp(['err on C:         ',num2str(err_C)])
    