% demo2 to test  bcdLrLr1_alsls algorithm.
% The purpose of this demo is to show how to exploit several initializations.
% The function bcdLrLr1_alsls has 2 kinds of optional input arguments to deal with initialization:
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
    K=30;
    L_vec=[2 1 3 2];
    R=length(L_vec);
    N=sum(L_vec);
    SNR=inf;                    % SNR [dB], choose SNR=inf for a noise-free model
    power_max=3;               % Each of the R terms is normalized and then multiplied by a value between 1 and power_max. 
                               % The highest this value, the more difficult the problem is (swamps are more likely)
    power_vec=linspace(1,power_max,R);   % Holds the set of R amplitude values
       
    
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
% Compute the BCD-(Lr,Lr,1) with a limited number of iterations for each starting point
% and then use the best starting point to go on until convergence
%-------------------------------------------------------------------------------
Ninit=50;
MaxItInit=10;
% Perform at most MaxItInit iterations for each starting point
[A_init,B_init,C_init,phi,itInit]=bcdLrLr1_alsls(X,L_vec,[],[],[],MaxItInit,[],[],Ninit);
% Then go on with the best starting point
[A_est,B_est,C_est,phi,it1]=bcdLrLr1_alsls(X,L_vec,[],[],[],[],[],[],[],A_init,B_init,C_init);

% compute relative errors
[err_A]=solve_blockperm(A_est,A,L_vec);
[err_B]=solve_blockperm(B_est,B,L_vec);
[err_C]=solve_blockperm(C_est,C,ones(1,R));
disp(['err on A:         ',num2str(err_A)])
disp(['err on B:         ',num2str(err_B)])
disp(['err on C:         ',num2str(err_C)])
disp(['phi:              ',num2str(phi)])
disp(['itInit:           ',num2str(itInit)])
disp(['it1:              ',num2str(it1)])
