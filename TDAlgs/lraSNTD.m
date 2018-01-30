function [ Ydec fit] = lraSNTD( Y, opts )
%% Nonnegative Tucker Decomposition Incorporating Low-rank Approximation
%   Usage [Ydec fit]=lraSNTD(Y,opts);
%    opts.
%        .NumOfComp: vector for the number of cols of each factor
%        .maxiter .maxiniter: maximum of iteration number  (internal max iter))
%        .lra : pca|sampling|randpca low-rank approximation algs. 
%        .alg : apg|hals|mult NMF algs
%        .update: true|false independent extraction or
%        extraction-and-update
%        .nncore: true|false nonnegative core tensor or not
%        .dimorder: order of the dim to be extracted
%        .tol
%        .trackit : step for checking the status
%
% Please cite:
%    Guoxu Zhou; Cichocki, A.; Shengli Xie; , "Fast Nonnegative Matrix/Tensor Factorization Based on Low-Rank Approximation," 
%    IEEE Transactions on Signal Processing, vol.60, no.6, pp.2928-2940, June 2012
%    doi: 10.1109/TSP.2012.2190410
%    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166354&isnumber=6198804
%
%
Ydim=size(Y);
N=numel(Ydim);
defopts=struct('NumOfComp',ones(1,N),'maxiter',500,'maxiniter',10,'dimorder',[],'tol',1e-6,'trackit',20,'norm',1,...
    'lra','pca','alg','apg','update',true,'alpha',0,'SamplingPara',1,'nncore',true);
if ~exist('opts','var')
    opts=struct();
end

[R,maxiter,maxiniter,dimorder,tol,trackit,nor,lra,alg,update,alpha,rs,nncore]=scanparam(defopts,opts);

if isempty(dimorder)
    [temp dimorder]=sort(Ydim,'descend');
end
dims=1:N;
[temp i]=setdiff(dimorder,dims(isinf(R)));
dimorder=dimorder(sort(i,'ascend'));
inforder=setdiff(dims,dimorder);



lra_opts=struct('NumOfComp',[],'maxit',maxiter,'maxiniter',maxiniter,'tol',tol,'alpha',alpha,'alg',alg,'lra',lra,'trackit',trackit,'lra_rank',rs);
Y=tensor(Y);
T=Y;
if update
    for m=dimorder
        if ~isinf(R(m))
            Ym=double(tenmat(Y,m));
            lra_opts.NumOfComp=R(m);
%             U{m}=rslraNMF(Ym,lra_opts);
            U{m}=lraNMF(Ym,lra_opts);

            U{m}=datanormalize(U{m},nor);
            Y=ttm(Y,pinv(U{m}),m); 
        else
            R(m)=Ydim(m);
            U{m}=speye(R(m),R(m));
        end
    end
else
    for m=dimorder
        if ~isinf(R(m))
            Ym=double(tenmat(Y,m));
            lra_opts.NumOfComp=R(m);
            U{m}=rslraNMF(Ym,lra_opts);
            U{m}=datanormalize(U{m},nor);
        else
            R(m)=Ydim(m);
            U{m}=speye(R(m),R(m));
        end
    end
end


if ~isempty(inforder)
    for i=inforder
        R(i)=Ydim(i);
        U{i}=speye(R(i));
    end
end
        

if nncore
    nc=ncore(T,U);
else
    iU=cellfun(@pinv,U,'uni',false);
    nc=ttm(T,iU,dimorder);
end
Ydec=ttensor(tensor(nc),U);


%% compute the fitting error
if nargout==2
    normT=norm(T);
    fit=1-norm(tensor(Ydec)-tensor(T))/normT;
end

    function nc=ncore(Y,U)
        nc=rand(R);
        enum=max(double(ttm(Y,U(dimorder),dimorder,'t')),eps);
        UUT=cellfun(@(x) x'*x,U(dimorder),'uni',false);
        for i=1:200
            nc=nc.*(enum./max(double(ttm(tensor(nc),UUT,dimorder)),eps));
        end
    end

end

