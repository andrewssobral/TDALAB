function [rk gs es]=moderankdect(X,n,opts)
%% mode rank detection based on SORTE using randomly sampled fibers
defopts=struct('mode','PCA','p0',500,'SamplingPara',.5,'SamplingThres',1e-3);
if ~exist('opts','var')
    opts=struct;
end
[mode p0 SamplingPara SamplingThres]=scanparam(defopts,opts);

X=double(X);
Xdim=size(X);
NumOfMode=numel(Xdim);
DIM=prod(Xdim([1:n-1 n+1:end]));
        
if min(DIM,Xdim(n))<p0
    z=double(tenmat(tensor(X),n));
    [rk gs es]=sorte(z);
    return;
end

z=double(tenmat(tensor(X),n));
mode=lower(mode);
switch mode
    case 'pca'
        if Xdim(n)<DIM
            z=eig(cov(double(z)'));
        else
            z=eig(cov(double(z)));
        end
        [rk gs es]=sorte(z);
    case 'random'        
        if Xdim(n)<DIM
            ncol=floor(DIM*SamplingPara);
            o=randperm(DIM);
            z=z(:,o(1:ncol));
            if ncol<Xdim(n)
                z=eig(cov(z));
            else
                z=eig(cov(z'));
            end
        else
            ncol=floor(Xdim(n)*SamplingPara);
            o=randperm(Xdim(n));
            z=z(o(1:ncol),:);
            if ncol<DIM
                z=eig(cov(z'));
            else
                z=eig(cov(z));
            end
        end
        [rk gs es]=sorte(z);    
        
    case 'incremental'
        NumOfSampled=p0;
       
        z=zeros(NumOfSampled,Xdim(n));
        index=1:NumOfMode;index(n)=[];
        var=cell(1,NumOfMode);
        var{n}=1:Xdim(n);
        
        o=randperm(DIM);
        curri=0;rk0=0;
        for t=1:DIM
            [var{index}]=ind2sub(Xdim([1:n-1 n+1:end]),o(t));
            zt=reshape(double(X(var{:})),1,Xdim(n));
            if max(abs(zt))>SamplingThres
                curri=curri+1;
                z(curri,:)=zt;
            end
            if curri<p0
                continue;
            end
            if ~rem(curri+1,50)
                [eigvect eigv]=eig(cov(z(1:curri,:)'));
                [eigv]=sort(diag(eigv),'descend');
                [rk gs es]=sorte(eigv);
                
                if (rk<=rk0&&rk*SamplingPara+p0<curri)||(curri>Xdim(n))
                    z=eigv;
                    break;
                else
                    rk0=rk;
                end
            end
        end
        
        
    otherwise
        error('Unsupported mode in component number estimation.');
end
   
end