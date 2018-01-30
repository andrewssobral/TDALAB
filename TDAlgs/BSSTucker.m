function [ Y nfit] = BSSTucker(X,opts )
%% CP incorporating BSS extraction / Multilinear BSS
%       based on CP_ALS @ tensor toolbox 2.4
% function [P,U,fit] = bsscpals(X,opts)
% X : tensor / ktensor / full tensor /double
% opts.BSSAlgID: BSS algorithms specification for the modes specified by
%                opts.ModesForBSS. length(.BSSAlgID)=length(.ModesForBSS)
%     .SamplingMode = 'times|percentage|FSTD1|FSTD2';
%     .SamplingPara = p: 1<=p, p*J samples are used if SamplingMode='times';
%                     p: 0<p<1, p*J^(N-1) samples are used if
%                     SamplingMode='Percentage';
%                   = p: interger LARGER than J(n) for FSTD.
%     .Tol:  threshold of stop

Xdim=size(X);
NumOfMode=numel(Xdim);

defopts = struct('NumOfComp',0,'MaxIter',1,'DimOrder',1:NumOfMode,'norm',2,...
    'Verbose',false,'PMFalgIDs',ones(1,NumOfMode),'SamplingMode','times','SamplingPara',10,'SamplingThres',0.01,...
    'update',false); 
if ~exist('opts','var')
    opts = struct;
end

[J,MaxIter,dimorder,nor,verbose,PMFalgIDs,SamplingMode,SamplingPara,SamplingThres,update] ...
    =scanparam(defopts,opts);


Xtype=class(X);
if ~strcmpi(Xtype,'double')&&~strcmpi(Xtype,'tensor')
    error('Full tensor input is expected.');
end

A=cell(1,NumOfMode);iA=A;
freeComps=setdiff(1:NumOfMode,dimorder);
for i=freeComps
    A{i}=speye(Xdim(i),Xdim(i));
    iA{i}=A{i};
end

%% extract in the order of dimorder
PMFalgIDs=PMFalgIDs-1;

PMFalgOrder=zeros(1,NumOfMode);
flag=PMFalgIDs>0;
PMFalgOrder(flag>0)=1:sum(flag);

flag=PMFalgIDs(dimorder)>0;
PMFindices=dimorder(flag);
NPMFindices=dimorder(~flag);


global PMFActParas;
[temp PMFalgs]=PMFalgInit();



if verbose
    fprintf('Begin to extract mode matrices using penalized matrix factorization ...\n');
end

ts=tic;
if strcmpi(SamplingMode,'times')
    for n=PMFindices 
        if verbose
            fprintf('PMF: Mode %d will be extracted by using [%s].\n',n,PMFalgs(PMFalgIDs(n)).details);
            fprintf('Settings of the algorithm:\n');
            disp(PMFActParas{PMFalgOrder(n)});
        end  
        z=RSDR; 
        curralg=str2func(PMFalgs(PMFalgIDs(n)).name);
        
        if ~isempty(J)
            PMFActParas{PMFalgOrder(n)}.NumOfComp=J(n);
        end
        
        %% rPCA;
%         [z Ee]=exact_alm_rpca(z,.1,1e-7,1000);
        
        A{n}=curralg(z,PMFActParas{PMFalgOrder(n)});
        A{n}=datanormalize(A{n},nor);
        
        iA{n}=pinv(A{n});
        
        if update
            X=ttm(X,iA{n},n);
            Xdim(n)=J(n);
        end
    end  %
    
    time=toc(ts);    
    if verbose
        fprintf('PMF complete. Total %s seconds are consumped.\n',time);
    end
elseif strcmpi(SamplingMode,'FSTD1')||strcmpi(SamplingMode,'FSTD2')
    FSTDopts.Rank=round(SamplingPara);
    FSTDopts.Verbose=verbose;
    if strcmpi(SamplingMode,'FSTD1')
        try
            Ysap=FSTDT1opt(X,FSTDopts);
        catch ME
            error('FSTD1 does not work. Please use another sampling algorithm.');
        end
    else
        try
            Ysap=FSTDT2rnd(X,FSTDopts);
        catch ME
            error('FSTD2 does not work. Please use another sampling algorithm.');
        end
    end
    for n=PMFindices   
        curralg=str2func(PMFalgs(PMFalgIDs(n)).name);
        
        if ~isempty(J)
            PMFActParas{PMFalgOrder(n)}.NumOfComp=J(n);
        end
        
        
        A{n}=curralg(Ysap.U{n},PMFActParas{PMFalgOrder(n)});        
        A{n}=datanormalize(A{n},nor);
        iA{n}=pinv(A{n});
        
        if update
            X=ttm(X,iA{n},n);
            Xdim(n)=J(n);
        end
        
        time=toc(ts);
        if verbose
            fprintf('PMF: Mode %d has been estimated by using [%s].\n',n,PMFalgs(PMFalgIDs(n)).details);
            disp(PMFActParas{PMFalgOrder(n)});
        end
    end  %
    if verbose
        fprintf('PMF complete. Total %s seconds are consumped.\n',time);
    end
else    
    fprintf(2,'Unknown sampling parameter in MBSS_Tucker.\n');
    return;
end

if ~isempty(NPMFindices)
    if verbose
        fprintf('Begin to estimate the other modes...\n');
    end
    
    %% initialize for other mode matrices
    for n=NPMFindices
        A{n}=orth(randn(Xdim(n),J(n)));
    end
    
    %iA(PMFindices)=cellfun(@pinv,A(PMFindices),'UniformOutput',false);
% 	iA(NPMFindices)=cellfun(@(x) x',A(NPMFindices),'UniformOutput',false);
    iA(freeComps)=A(freeComps);
    for it=1:MaxIter
        if numel(NPMFindices)>1
            if it==1&&~update&&~isempty(PMFindices)
                X=ttm(X,iA(PMFindices),PMFindices);
            end
            
            for n=NPMFindices
                indNP=setdiff(NPMFindices,n);
                Utilde=ttm(X,A(indNP),indNP,'t');
                
                A{n}=nvecs(Utilde,n,J(n));            
                iA{n}=(A{n})';
                
            end % other parts
            
            G = ttm(Utilde, A, n, 't');
        else
            n=NPMFindices;
            if update
                Utilde=X;
            else
                Utilde=ttm(X,iA,-n);
            end
            A{n}=nvecs(Utilde,n,J(n));            
            iA{n}=(A{n})';
            G = ttm(Utilde, A, n, 't');
            break;
        end
    end % end of it
else
    if update
        G=X;
    else
        G=ttm(X,iA);
    end
end

Y=ttensor(G,A);

    function [z]=RSDR
        Xdim=size(X);
        DIM=prod(Xdim([1:n-1 n+1:end]));
        NumOfSampled=round(SamplingPara*J(n));
        
        if DIM<500||NumOfSampled/DIM>0.7||NumOfSampled>5000
            z=double(tenmat(tensor(X),n));
            return;
        end
        NumOfSampled=max(floor(NumOfSampled),J(n));
        
        z=zeros(NumOfSampled,Xdim(n));
        index=1:NumOfMode;index(n)=[];
        var=cell(1,NumOfMode);
        var{n}=1:Xdim(n);
        
        o=randperm(DIM);
        curri=0;
        for t=1:DIM
            [var{index}]=ind2sub(Xdim([1:n-1 n+1:end]),o(t));
            zt=reshape(double(X(var{:})),1,Xdim(n));
            if max(abs(zt))>SamplingThres
                curri=curri+1;
                z(curri,:)=zt;
            end
            if curri>NumOfSampled
                break;
            end
        end
        z=z(1:curri,:)';
    end

end


