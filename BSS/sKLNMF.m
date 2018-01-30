function [ B H t] = sKLNMF( X,r,alpha)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
TolX=1e-6;
MaxIt=2e+3;
default_alpha=0.1;

ts=clock;
[M L]=size(X);
if nargin==1
    r=M;
    alpha=default_alpha;
elseif nargin==2
    alpha=default_alpha;
end

%% initialization
% [B H]=nndsvd(X,r,0);
B=rand(M,r);
H=rand(r,L);

for i=1:MaxIt
    B0=B;
    H=H.*(B'*(X./(B*H+eps)));
    H=H/(1+alpha);
    
    C=B.*((X./max(B*H,eps))*H');
    B=C./repmat(sum(H,2)',M,1);
    sumB=max(sum(B),eps);
    B=B./sumB(ones(M,1),:);
    sumB=sumB';
    H=H.*sumB(:,ones(1,L));
    if norm(B0-B,1)<TolX
        break;
    end
end
t=etime(clock,ts);
