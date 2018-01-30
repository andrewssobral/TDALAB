function [S,w]=PMFwasobi(x,opts)
%% call wasobi
% 
defopts=struct('NumOfComp',[],'AR_order',1,'rmax',1);
if ~exist('opts','var')
    opts=struct;
end
[r AR_order rmax]=scanparam(defopts,opts);

M=size(x,2);
if isempty(r)
    r=M;
end
x0=x;


[z temp1 temp2]=svds(x,r,'L');
z=z';

S=wasobi(z,AR_order,rmax);
S=S';
w=S\x;

function [signals]= wasobi(x,AR_order,rmax);
%
% implements algorithm WASOBI for blind source separation of
% AR sources.
% INPUT:    x .... input data matrix d x N
%           d .... dimension of the data
%           N .... length of the data
%           ARmax .. maximum AR order of the separated sources
%           rmax ... a constant that may help to stabilize the algorithm.
%                    it has the meanong of maximum magnitude of poles of the
%                    AR sources. The choice rmax=1 means that no stabilization 
%                    is applied.
% OUTPUT: Wwasobi ...... estimated de-mixing matrix
%         Wsobi ........ initial estimate of the matrix obtained by SOBI
%         ISR .......... estimated ISR matrix which represents approximate accuracy 
%                        of the separation provided that there is no additive 
%                        noise in the model.
%
num_of_iterations = 3;
[d N]=size(x);

Xmean=mean(x,2);
x=x-Xmean*ones(1,N);  %%%%%%%%%  removing the sample mean
[AOL_1 Rx_est] = sobi(x,AR_order);
AOL_init = AOL_1;
AOL = AOL_1;
Rx_est = [];
for in = 1:num_of_iterations
    [AOL_1,Rx_est,AR]=newasobi(x,AR_order+1,AOL_1,Rx_est,rmax);
    AOL = AOL*AOL_1;
end
ISR=CRLB4(AR)/N;
Wwasobi=inv(AOL);
Wsobi=inv(AOL_init);
poradi=cluster(ISR);       %%% performs "spectral clustering"
Wwasobi=Wwasobi(poradi,:);
Wsobi=Wsobi(poradi,:);
ISR=ISR(poradi,poradi);
signals=Wwasobi*x+(Wwasobi*Xmean)*ones(1,N);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

function [A_est_SOBI Rx_est]=sobi(x,AR_order);
format long 
[d,N]=size(x);
T=length(x(1,:))-AR_order;
Rx_est=corr_est(x,T,0);

% SOBI algorithm
[H1,D]=eig(Rx_est);
W_est=[];
for k=1:d
    W_est=[W_est;D(k,k)^(-.5)*H1(:,k)'/(H1(:,k)'*H1(:,k))];
end
Rz_est=corr_est(W_est*x,T,AR_order);
U_est=eye(d);
md=d*(AR_order+1);
for in=1:2*d
    for k=1:d
        for j=[1:k-1 k+1:d]
         %%   G=squeeze([Rz_est(k,k,s)-Rz_est(j,j,s); Rz_est(k,j,s)+Rz_est(j,k,s)])';
            G=[Rz_est(k,k+d:d:md)-Rz_est(j,j+d:d:md); Rz_est(k,j+d:d:md)+Rz_est(j,k+d:d:md)]';
            [V D]=eig(G'*G);
            [index1,index2]=max(diag(D));
            v=V(:,index2)/(V(:,index2)'*V(:,index2));
            teta=-.5*asin(v(2));
            if (abs(cos(2*teta)-v(1))>10^-8)
                teta=pi/2-teta;
            end
            U_est0=[cos(teta) sin(teta); -sin(teta) cos(teta)];
            U_est(:,[k j])=U_est(:,[k j])*U_est0;
            aux=Rz_est(:,k+d:d:md)*cos(teta)-Rz_est(:,j+d:d:md)*sin(teta);
            Rz_est(:,j+d:d:md)=Rz_est(:,k+d:d:md)*sin(teta)+Rz_est(:,j+d:d:md)*cos(teta);
            Rz_est(:,k+d:d:md)=aux;
            aux=Rz_est(k,d+1:md)*cos(teta)-Rz_est(j,d+1:md)*sin(teta);
            Rz_est(j,d+1:md)=Rz_est(k,d+1:md)*sin(teta)+Rz_est(j,d+1:md)*cos(teta);
            Rz_est(k,d+1:md)=aux;
        end
    end
end
A_est_SOBI=inv(W_est)*U_est;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

function R_est=corr_est(x,T,q,A_pre)
%
[K,KK]=size(x);
if ~exist('A_pre')
    A_pre_inv = eye(K);
else
    A_pre_inv = inv(A_pre);
end
R_est=[];
for index=1:q+1
%   R_est(:,:,index)=1/T*x(:,1:T)*x(:,index:T+index-1)';
%   R_est(:,:,index) = A_pre_inv*R_est(:,:,index)*A_pre_inv';
%   r_est=[r_est;reshape(R_est(:,:,index),1,K^2)];
   R_aux=1/T*A_pre_inv*(x(:,1:T)*x(:,index:T+index-1)')*A_pre_inv';
   R_est=[R_est R_aux];
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

function [A_est_WASOBI,Rx_est,AR]=newasobi(x,M,A_pre,Rx_est1,rmax);
%
%   Computes one iteration needed for wasobi
%
%   d ... number of independent components
%   M ... number of lags = AR_order+1
%   Vnorm ... estimated demixing matrix obtained in preprocessing
%   Wnorm ...  estimated demixing matrix obtained by the advanced tuning
%

status=1;
[d,N]=size(x);
d2=d*(d+1)/2;
T=length(x(1,:))-M+1;
if isempty(Rx_est1)
    Rx_est=corr_est(x,T,M-1,A_pre);
else
   A_pre_inv = inv(A_pre);
   for index=1:M
       id=(index-1)*d;
       Rx_est(:,id+1:id+d)=A_pre_inv*Rx_est1(:,id+1:id+d)*A_pre_inv';
   end
end
R=zeros(M,d);
for index=1:M
    id=(index-1)*d;
    R(index,:)=diag(Rx_est(:,id+1:id+d)).'; %%% columns of R will contain 
                           %%% covariance function of the separated components
end
%
[AR,sigmy]=armodel(R,rmax);      %%% compute AR models of estimated components
%
Y=zeros(d2,M);      %%%% will be transformed data - lower triangle of lagged covariance
%CVinv=zeros(M,M,d2); %%%% will be covariance matrices of rows of Y
CVinv=[];                %%%%  Now, it will have dimension zeros(M,M*d2).
im=1; Md=M*d; condnu=[];
for i=1:d
  for k=1:i
%     Y(im,:)=squeeze((Rx_est(i,k,:)+Rx_est(k,i,:))/2)';  %%% computing Y and CV
     Y(im,:)=(Rx_est(i,k:d:Md)+Rx_est(k,i:d:Md))/2;  
     CV=chyba2(AR(:,i),AR(:,k),sigmy(1,i),sigmy(1,k),d,M);
     [V,D] = eig(CV);
     dd=diag(D);
     if max(dd)>min(dd)*1e10
        dd=dd+max(dd)*1e-10;
     end   
     CV=V*diag(1./dd)*V';
     CVinv=[CVinv CV*(2-(i==k))];
     im=im+1;
  end
end  
% condnu
%%A_est_WASOBI=serial5(Y, CVinv, d, M); %% TO MINIMIZE THE WEIGHTED
%%                                                         CRITERION
% function A0=serial5(Y,CVinv,d,M); 
A0=eye(d);
H=zeros(M*d,M*d);
g=zeros(M*d,1);
gx=0;
m=0;
for id=1:d        %%%%%% sums d2 contributions to H, g, gx
    for id2=1:id
        m=m+1; im=(m-1)*M;
%       Wm=squeeze(CVinv(:,:,m)); % optimum weighting matrix = inverse of covariance m.
        Wm=CVinv(:,im+1:im+M);
        aux=kron(eye(M),(A0(id,:).*A0(id2,:))');
        H=H+aux*Wm*aux';                         %% size Md x Md
        g=g+aux*Wm*Y(m,:)';
        gx=gx+Y(m,:)*Wm*Y(m,:)';
    end
end  
Lambda=H\g;        % Optimum parameter vector lambda in (3.7) for given A
value=gx-g'*Lambda;
clear H g
%
para=[A0(:); Lambda(d+1:end,1)];  %%% initial parameter for iteration
md=(M-1)*d;
%
ff=[];
num_of_iteration = 3;
while (length(ff)<num_of_iteration)
    ff=[ff value];
    B11=zeros(d^2,d^2);
    B12=zeros(d^2,md);
    B22=zeros(md,md);
    C1=zeros(d^2,1);
    C2=zeros(md,1);
    LaMat=reshape(Lambda,d,M);
    m=0;
    for id=1:d        
      for id2=1:id
          m=m+1; im=(m-1)*M;
          %%% Wm=squeeze(CVinv(:,:,m)); 
          Wm=CVinv(:,im+1:im+M);
          aux=kron(eye(M),(A0(id,:).*A0(id2,:))');
          Hkl=zeros(M,d^2); 
          Hkl(:,id:d:end)=LaMat'.*(ones(M,1)*A0(id2,:));
          Hkl(:,id2:d:end)=Hkl(:,id2:d:end)+LaMat'.*(ones(M,1)*A0(id,:));
          Gkl=aux(d+1:end,:)';
          rkl=Y(m,:)-Lambda'*aux;
          HW=Hkl'*Wm;
          B11=B11+HW*Hkl;
          B12=B12+HW*Gkl;
          B22=B22+Gkl'*Wm*Gkl;
          C1=C1+HW*rkl';
          C2=C2+Gkl'*Wm*rkl';
      end
    end
    para=para+[B11 B12; B12' B22]\[C1; C2];
    A0=reshape(para(1:d^2,1),d,d);
    Lambda(d+1:end,1)=para(d^2+1:end,1);
 %   value=criter2(A0, Y, CVinv, d, M);
end
A_est_WASOBI=A0;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

function [AR,sigmy]=armodel(R,rmax)
%
% to compute AR coefficients of the sources given covarience functions 
% but if the zeros have magnitude > rmax, th ezaros are pushed back.
%
[M,d]=size(R); v1=zeros(1,d); v2=v1; Rs=R;
for id=1:d
    AR(:,id)=[1; -toeplitz(R(1:M-1,id),R(1:M-1,id)')\R(2:M,id)];
   % AR(:,id)=(polystab(AR(:,id)'))';  %%%% STABILIZE THE RESULT
    v=roots(AR(:,id)); %%% mimicks the matlab function "polystab"
%    v1(1,id)=max(abs(v));
    vs=0.5*(sign(abs(v)-1)+1);
    v=(1-vs).*v+vs./conj(v);
    vmax=max(abs(v));
%    v2(1,id)=max(abs(v));
    if vmax>rmax
       v=v*rmax/vmax;
    end   
    AR(:,id)=poly(v)'; %%% reconstructs back the covariance function
    A_1 = toeplitz(AR(:,id),[1 zeros(1,M-1)]);
    A_2 = zeros(M);
    A_2(1:M-1,2:M) = hankel(AR(2:M,id));
    Rs(:,id) = inv(A_1+A_2)*[1;zeros(M-1,1)]; %% computes covariance functions
end 
% sigmy=sum(AR.*R)
sigmy=R(1,:)./Rs(1,:);
% [v1; v2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi=chyba2(AR1,AR2,sigma1,sigma2,d,M);
%
% to compute approximate asymptotic covariance matrix of lagged products of
% two autoregressive processes by residue calculus via Eran Doron. 
%
p_den=2*M-2;
sigma=sigma1*sigma2;
den=conv(AR1',AR2');
A_1 = toeplitz(den.',[1 zeros(1,p_den)]);
A_2 = zeros(p_den+1);
A_2(1:p_den,2:p_den+1) = hankel(den(2:p_den+1));
phi = (A_1+A_2)\[sigma;zeros(p_den,1)];
Phi=zeros(M,M);
for k=1:M
    ind1=abs((0:M-1)-k+1)+1;
    ind2=(0:M-1)+k;
    Phi(k,:)=phi(ind1')'+phi(ind2')'; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ISR = CRLB4(AR);
%
% CRLB4(AR) generates the CRLB for gain matrix elements (in term 
% of ISR) for blind separation of K Gaussian autoregressive sources 
% whose AR coefficients (of the length M, where M-1 is the AR order)
% are stored as columns in matrix AR.

[M K]=size(AR);

for k=1:K
    A_1 = toeplitz(AR(:,k),[1 zeros(1,M-1)]);
    A_2 = zeros(M);
    A_2(1:M-1,2:M) = hankel(AR(2:M,k));
    Rs(:,k) = inv(A_1+A_2)*[1;zeros(M-1,1)]; %% computes covariance functions
end 

A_inv=eye(K);

sum_Rs_s=zeros(K,K);

for k=1:K
    for s=0:M-1
        for t=0:M-1
            sum_Rs_s(k,:)=sum_Rs_s(k,:)+AR(s+1,k)*AR(t+1,k)*Rs(abs(s-t)+1,:);
        end
    end
end

denom=sum_Rs_s'.*sum_Rs_s+eye(K)-1;
ISR=sum_Rs_s'./denom.*(ones(K,1)*Rs(1,:))./(Rs(1,:)'*ones(1,K));
ISR(eye(K)==1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function poradi=cluster(ISR)
%
% performs spectral clustering
%
ISRmax=max(sum(ISR,2));
M=ISR/(2*ISRmax);
M=M+diag(-sum(M,2)+1);
[V D]=eig(M);
[y is]=sort(-abs(diag(D)));
evec=V(:,is(2));
[y poradi]=sort(real(evec)+imag(evec));

