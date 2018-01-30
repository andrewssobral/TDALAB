function [ Class sample training] = tuckerDA( sample,training,group,options )
%%Tucker Discriminative Analysis
% coded by Guoxu Zhou
% E-Mail: zhouguoxu@gmail.com

defopts=struct('NumOfComp',[],'maxiter',100,'tol',1e-6,'Eigs','ratio','classifier','lda','distance','euclidean','classify_rule','default');
if ~exist('options','var')
    options=struct();
end

[J maxiter tol proj_mode classifier distance classify_rule]=scanparam(defopts,options);
training=tensor(training);
tsize=size(training);
N=numel(tsize);
nTr=numel(group);
one2nTr=1:nTr;


if tsize(end)~=nTr
    error('The length of GROUP must equal the last dimension of TRAINING.');
end

labels=unique(group);
C=numel(labels);

idx=cell(1,N);
U=cell(1,N-1);
for n=1:N-1
    idx{n}=1:tsize(n);
    U{n}=randn(J(n),tsize(n));
end

%% training stage
training_m=mean(double(training),N);
training_b=cell(1,C);
training_w=cell(1,C);
nc=zeros(1,C);

for c=1:C
    nc(c)=sum(group==labels(c));
    idx{N}=one2nTr(group==labels(c));
    training_mc=tensor(mean(double(training(idx{:})),N));
    training_b{c}=training_mc-training_m;
    
    training_mc=reshape(training_mc,[tsize(1:N-1) 1]);
    training_w{c}=training(idx{:})-ttm(training_mc,ones(nc(c),1),N);
end

%% test
for it=1:maxiter
    err=0;
    for n=1:N-1
        u0=U{n};
        B=zeros(tsize(n),tsize(n));
        W=B;
        
        for c=1:C
            temp=double(tenmat(ttm(training_b{c},U,-n),n));
            B=B+nc(c)*(temp*temp');
            
            temp=double(tenmat(ttm(training_w{c},U([1:n-1 n+1:end]),[1:n-1 n+1:N-1]),n));
            W=W+temp*temp';
        end
        
        switch lower(proj_mode)
            case 'ratio'
            %% get U{n}
                [temp d]=eigs(B,W+B,J(n));
                temp=bsxfun(@rdivide,temp,sum(temp.^2).^.5);   
            case 'diff'
                %% get U{n} -- yet another way
                [~,zeta]=eigs(B,W,1,'LM');
                CoreMat=B-zeta.*W;
                [temp D]=eigs(CoreMat,J(n),'LM');
        end
        
        U{n}=temp';   
        U{n}=bsxfun(@times,U{n},sign(diag(U{n}*u0')));
        
        err=err+abs(trace(U{n}*u0')/J(n)-1);
    end % for n
%     fprintf('%d -- err=%f\n',it,err);
    if err<tol
        break;
    end
end

if it>=maxiter
    warning('Not perfectly converged. You may try a larger number of maxiter.');
end

%% classify
%% project
training=ttm(training,U,1:N-1);
training=double(tenmat(training,N));

sample=ttm(tensor(sample),U,1:N-1);
sample=double(tenmat(sample,N));

switch lower(classifier)
    case 'lda'        
        if strcmpi(classify_rule,'default')
            classify_rule='linear';
        end
        Class=classify(sample,training,group,classify_rule);
    case 'knn'
        if strcmpi(classify_rule,'default')
            classify_rule='nearest';
        end
        Class=knnclassify(sample,training,group,1,distance,classify_rule);
    otherwise
        error('Undefined classifier.');
end

% if sum(J(1:N-1))==N-1 %% without lda
%
% else  %% with lda
%     %% update
%     training=ttm(training,U,1:N-1);
%     training=double(tenmat(training,N));
%
%
%     sample=ttm(tensor(sample),U,1:N-1);
%     sample=double(tenmat(sample,N));
%
%     Class=knnclassify(sample,training,group);
% %     Class=classify(sample,training,group);
%
% %       opts.subgroups=4;
% %       Class=cobe_classify(sample,training,group,opts);
% end

end

