function [P res ] = swatld( Y,opts )
% Usage: [P res ] = swatld( Y,opts )
%  Opts: [NumOfComp],[Tol],[MaxIter];
%
% Ref: 
% A novel trilinear decomposition algorithm for second-order linear calibration
%   Zeng-Ping Chen, Hai-Long Wu, Jian-Hui Jiang, et al.
%   Chemometrics and Intelligent Laboratory Systems Vol.52(1)2000. 75-86
defopts = struct('NumOfComp',0,'Tol',1e-10,'MaxIter',1000);
if ~exist('opts','var')
    opts = struct;
end

[J,tol,maxIter]=scanparam(defopts,opts);

if strcmp(class(Y),'ktensor')||strcmp(class(Y),'ttensor')||strcmp(class(Y),'tensor')
    Y=double(Y);
end

if strcmp(class(Y),'cell')
    if length(Y)~=3
        error('only works for PARAFAC/CP-3 model.');
    end
    [II]=size(Y{1},1);JJ=size(Y{2},1);KK=size(Y{3},1);
    for k=1:KK
        Xk{k}=Y{1}*diag(Y{3}(k,:))*Y{2}';
    end
elseif strcmp(class(Y),'double')&&ndims(Y)==3
    [II,JJ,KK]=size(Y);
    for k=1:KK
        Xk{k}=Y(:,:,k);
    end
    Yi=reshape(Y,II,JJ*KK); % only for the convenience of calculating res
else
    error('Unsuported input data type.');
end

NN=J;

Xj=cell(1,JJ);
Xi=cell(1,II);

for j=1:JJ
    P=cellfun(@(x) x(:,j),Xk,'UniformOutput',false);
    Xj{j}=cell2mat(P)';
end
for i=1:II
    Q=cellfun(@(x) x(i,:)',Xk,'UniformOutput',false);
    Xi{i}=cell2mat(Q);
end

A{1}=randn(II,NN);
A{1}=bsxfun(@rdivide,A{1},max(sum(A{1}.*A{1}),eps).^0.5);
A{2}=randn(JJ,NN);
A{2}=bsxfun(@rdivide,A{2},max(sum(A{2}.*A{2}),eps).^0.5);

for it=1:maxIter
    
    % update A{3}='C'
    dA=max(sum(A{1}.*A{1}),eps)';dB=max(sum(A{2}.*A{2}),eps)';
    iA=pinv(A{1});iB=pinv(A{2});
    cC=cellfun(@(x) 0.5.*(diag(iB*x'*A{1})./dA+...
        diag(iA*x*A{2})./dB),Xk,'UniformOutput',false);
    A{3}=cell2mat(cC)';
    
    % update A{2}='B'
    dC=max(sum(A{3}.*A{3}),eps)';
    iC=pinv(A{3});
    cB=cellfun(@(x) 0.5.*(diag(iA*x'*A{3})./dC+...
        diag(iC*x*A{1})./dA),Xj,'UniformOutput',false);
    A{2}=cell2mat(cB)';
    
    
    % update A{1}='A'    
    dB=max(sum(A{2}.*A{2}),eps)';
    iB=pinv(A{2});
    cA=cellfun(@(x) 0.5.*(diag(iC*x'*A{2})./dB+...
        diag(iB*x*A{3})./dC),Xi,'UniformOutput',false);
    A{1}=cell2mat(cA)';
    
      
    dA=max(sum(A{1}.*A{1}),eps)';
    A{1}=bsxfun(@rdivide,A{1},dA');
    A{2}=bsxfun(@rdivide,A{2},dB');
    A{3}=bsxfun(@times,A{3},dA'.*dB');
    
    if strcmp(class(Y),'cell')
        res(it)=abs(norm(ktensor(Y)-ktensor(A)));
    else
        res(it)=norm(Yi-A{1}*khatrirao(A{3},A{2})','fro');
    end
    
    if (it>10)
        if abs(res(it)-res(it-1))<tol
            break;
        end
    end

    
end % it

P=ktensor(A(:));


end

