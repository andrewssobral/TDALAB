function [P  res ] = call_APTLD( Y,opts )
% Usage: [P res]=aptld(Y,J,opts);
%     P: ktensor
%     opts. [Tol] [MaxIter] [q] [p] [r]
%          : p,q,r: 10^20 in default.
%
% Ref:
% Alternating penalty trilinear decomposition algorithm
% for second-order calibration with application to
% interference-free analysis of excitation–emission
% matrix fluorescence data

%% !!! This huge penality parameter is proposed and set by the author
penality=10^20;


defopts = struct('NumOfComp',0,'Tol',1e-10,'MaxIter',1000,'p',penality,'q',penality,'r',penality);
if ~exist('opts','var')
    opts = struct;
end
[J,tol,maxIter,p,q,r]=scanparam(defopts,opts);

if strcmp(class(Y),'ktensor')||strcmp(class(Y),'ttensor')||strcmp(class(Y),'tensor')
    Y=double(Y);
end
if strcmp(class(Y),'cell')
    if length(Y)~=3
        error('only works for PARAFAC/CP-3 model.');
    end
    
    II=size(Y{1});
    JJ=size(Y{2},1);
    KK=size(Y{3},1);
    X=[];
    for k=1:KK
        X=[X Y{1}*diag(Y{3}(k,:))*Y{2}'];
    end
elseif strcmp(class(Y),'double')&&ndims(Y)==3
    [II JJ KK]=size(Y);
    X=reshape(Y,II,JJ*KK);
else
    error('Unsuported input data type.');
end
NN=J;
[A{1}, A{2}, A{3}, res]=aptld(X,KK,NN,tol,maxIter,p,q,r);
P=ktensor(A(:));
%% Source code of aptld
%XPI=[X..1, X..2,..., X..k], size: (I * J )* K
%epsilon is the tolerance
%I is the row number
%J is the column number
%K is the channel number
%N is the estimated component number
%LFT is the loss function
%M is the iterative number
function [A, B, C, LFT, M]=aptld (XPI,K,N,epsilon,maxIter,p,q,r)
if nargin<4
    epsilon=10*eps*norm(XPI,1)*max(size(XPI));
end
[I, JK]=size(XPI); J=JK/K;
XIJK=reshape (XPI, I, J, K); %cut X along K direction
XJKI=shiftdim(XIJK,1); %cut X along I
XKIJ=shiftdim(XIJK,2); %cut X along J
% initialization of A & B
A=rand(I,N);B=rand(J,N);C=zeros(K,N);
TOL=10; M=0; LFT=[ ]; LF=0.01; 

%%  Definition of p q r
% p=10^20;q=10^20;r=10^20;
%   end of definition


% initialization of C
CD1=0;CD2=0;aa=0;bb=0;CD3=0;CD4=0;PB=pinv(B',epsilon);
PA=pinv(A',epsilon);Da=diag(ones(N,1)./diag(A'*A));
Db=diag(ones(N,1)./diag(B'*B));
for i=1:I
    CD1=CD1+XJKI(:,:,i)'* PB*Da*diag(A(i,:));
    aa=aa+diag(A(i,:))*Da*diag(A(i,:));
end
for j=1:J
    CD2=CD2+XKIJ(:,:,j)*PA*Db*diag(B(j,:));
    bb=bb+diag(B(j,:))*Db*diag(B(j,:));
    CD3=CD3+XKIJ(:,:,j)*A*diag(B(j,:));
    CD4=CD4+diag(B(j,:))*(A'*A)*diag(B(j,:));
end
C=(p*(CD1+CD2)+CD3)*pinv(p*aa+p*bb+CD4, epsilon);
%start to caculate LFT and do iteration
while TOL>epsilon && M<maxIter
    %estimation of B
    BD1=0;BD2=0; BD3=0;BD4=0;cc=0;aa=0;PC=pinv(C',epsilon);
    Dc=diag(ones(N,1)./diag(C'*C));
    for k=1:K
        BD1=BD1+XIJK(:,:,k)'*PA*Dc*diag(C(k,:));
        cc=cc+diag(C(k,:))*Dc*diag(C(k,:));
    end
    for i=1:I
        BD2=BD2+XJKI(:,:,i)*PC*Da*diag(A(i,:));
        aa=aa+diag(A(i,:))*Da*diag(A(i,:));
        BD3=BD3+XJKI(:,:,i)*C*diag(A(i,:));
        BD4=BD4+diag(A(i,:))*(C'*C)*diag(A(i,:));
    end
    B=(q*(BD1+BD2)+BD3)*pinv(q*cc+q*aa+BD4, epsilon);
    %estimation of A
    AD1=0;AD2=0;AD3=0;AD4=0;bb=0;cc=0;PB=pinv(B',epsilon);
    Db=diag(ones(N,1)./diag(B'*B));
    for j=1:J
        AD1=AD1+XKIJ(:,:,j)'*PC*Db*diag(B(j,:));
        bb=bb+diag(B(j,:))*Db*diag(B(j,:));
    end
    for k=1:K
        AD2=AD2+XIJK(:,:,k)*PB*Dc*diag(C(k,:));
        cc=cc+diag(C(k,:))*Dc*diag(C(k,:));
        AD3=AD3+XIJK(:,:,k)*B*diag(C(k,:));
        AD4=AD4+diag(C(k,:))*(B'*B)*diag(C(k,:));
    end
    A=(r*(AD1+AD2)+AD3)* pinv (r*bb+r*cc+AD4, epsilon);
    % normalization of A & B
    for n=1:N
        A(:,n)=A(:,n)./sqrt(sum(A(:,n).*A(:,n))); B(:,n)=B(:,n)./...
        sqrt(sum(B(:,n).*B(:,n)));
    end
    %estimation of C
    CD1=0;CD2=0;aa=0;bb=0;CD3=0;CD4=0;
    PA=pinv(A',epsilon);
    Da=diag(ones(N,1)./diag(A'*A));
    for i=1:I
        CD1=CD1+XJKI(:,:,i)'*PB*Da*diag(A(i,:));
        aa=aa+diag(A(i,:))*Da*diag(A(i,:));
    end
    for j=1:J
        CD2=CD2+XKIJ(:,:,j)*PA*Db*diag(B(j,:));
        bb=bb+diag(B(j,:))*Db*diag(B(j,:));
        CD3=CD3+XKIJ(:,:,j)*A*diag(B(j,:));
        CD4=CD4+diag(B(j,:))*(A'*A)*diag(B(j,:));
    end
    C=(p*(CD1+CD2)+CD3)*pinv(p*aa+p*bb+CD4, epsilon);
    %stopping criterion
    %calculate loss function
    LFTT=0;
    for k=1:K
        XXX(:,:,k)=A*diag(C(k,:))*B';
    end
    for k=1:K
        LFTT=LFTT+trace((XIJK(:,:,k)-XXX(:,:,k))'*(XIJK(:,:,k)- ...
        XXX(:,:,k)));
    end
    TOL=abs((LFTT-LF)/LF);
    LFT=[LFT,LFTT];
    LF=LFTT;
    M=M+1;
end
%post-processing to keep sign convention
[maxa,inda]=max(abs(A)); [maxb,indb]=max(abs(B));
asign=ones(N,1);bsign=ones(N,1);
for n=1:N
    asign(n)=sign(A(inda(n),n));
    bsign(n)=sign(B(indb(n),n));
end
A=A*diag(asign);B=B*diag(bsign);C=C*diag(asign)*diag(bsign);