function [ P ] = cp_smbss( Y,opts )
%% CPD based on a single mode blind source separation (CP-SMBSS)
% Usage: P = cp_smbss(Y,opts);
% opts:   -- a structure containing the following fields
%        .NumOfComp:  Rank of the output tensor
%        .BSSmode:    the mode on which BSS runs
%        .BSSalgFile: a -mat file specifying the algorithm and its parameters.
%        .lra:        Low-rank approximation mode (PCA|sampling|none)
%        .lra_rank:   rank of low-rank approximation. It should be no less than NumOfComp
%        .krp_init:   random|svd
%        .krp_maxit:  maximum iteration number of Khatri-Rao product projection
%        .constraints: nonnegative|none
%
% Refs:
% G. Zhou and A. Cichocki, Canonical polyadic decomposition based on a single mode blind source separation, 
%     IEEE Signal Processing Letters, vol. 19, no. 8, pp. 523-526, Aug. 2012.
%
% Last update: 2013.04.22
% E-mail: zhouguoxu@brain.riken.jp
% Guoxu Zhou


Ydim=size(Y);
N=numel(Ydim);
defopts=struct('NumOfComp',[],'BSSmode',1,'BSSalgFile','fileName','lra','pca','lra_rank',[],...
    'krp_init','random','krp_maxit',50,'krp_costol',1e-2,'constraints','none','cons_para',0,'verbose',true);
if ~exist('opts','var')
    opts=struct();
end
[J,pmode,BSSalgFile,lra,lra_rank,krp_init,krp_maxit,krp_costol,cons,cons_para,verbose]=scanparam(defopts,opts);
if isempty(J)
    error('opts.NumOfComp should be specified. Now it is empty.');
end

if strcmp(class(Y),'tensor')~=1
    Y=tensor(Y);
end

if length(pmode)>1
    fprintf('[CP_OMP] warning: only one mode will be extracted by using BSS (PMF).\n');
    pmode=pmode(1);
end


%% Extract number of dimensions and norm of X.


%% Set up and error checking on initial guess for U.
% Observe that we don't need to calculate an initial guess for the
% first index in dimorder because that will be solved for in the first
% inner iteration.
P=cell(N,1);

%% RunPMF
PMFalg=load(BSSalgFile);
if numel(PMFalg.alg)>1
    PMFalg.alg=PMFalg.alg{1};
    PMFalg.algname=PMFalg.algname{1};
    PMFalg.algopts=PMFalg.algopts{1};
    fprintf(2,'PMF: %s contains multiple algorithms. Only [%s] will be used.\n',PMFalg.algname);
end
if verbose
    fprintf('PMF: Mode %d will extracted.\n Algorithm: [%s].\n',pmode,PMFalg.algname);
    fprintf('Parameters of the algorithm:\n');
    disp(PMFalg.algopts);
end  

% if isfield(PMFalg.algopts,'NumOfComp')&&~isempty(J)
if ~isempty(J)    
    PMFalg.algopts.NumOfComp=J;
end
Yp=double(tenmat(Y,pmode));


switch lra
    case {'pca','randpca'}
        if isempty(lra_rank) lra_rank=J; end
        Ypl=lowrankapp(Yp,lra_rank,lra);
    case 'sampling'
        if isempty(lra_rank) lra_rank=J; end
        lra_rank=min(lra_rank,size(Yp,2));
        o=randperm(size(Yp,2));
        Ypl=Yp(:,o(1:lra_rank));
    otherwise
        lra='none';
        Ypl=Yp;
end


%% 2D-NMF
if N==2        
    [P{pmode} temp]=PMFalg.alg(Ypl,PMFalg.algopts);
    resdim=setdiff(1:N,pmode);
    if strcmpi(lra,'none')
        P{resdim}=temp';
        P=ktensor(P);
    else
        if strcmpi(cons,'nonnegative')
            str=horzcat('Nonnegative constraint on mode-',...
                num2str(resdim),' is solved by using ''nlssolver''. Otherwise please set ''lra=none''.');
            warning(str);
            P{resdim}=nlssolver(Yp,P{pmode});
%             P{resdim}=zeros(J,size(Yp,2));
%             for col=1:size(Yp,2)
%                 P{resdim}(:,col)=lsqnonneg(P{pmode},Yp(:,col));
%             end
            P{resdim}=P{resdim}';
            P=ktensor(P);
        else
            Yp=(P{pmode}\Yp);
            P{resdim}=Yp';
            P=ktensor(P);
        end
    end
    return;
end

%% N-way CP
P{pmode}=PMFalg.alg(Ypl,PMFalg.algopts);
clear Ypl;
P{pmode}=datanormalize(P{pmode},2);
Yp=(P{pmode}\Yp)';  %-- Yp=khatrirao(A(-pmode)) --%
lambda=sum(Yp.*Yp).^0.5;
Yp=bsxfun(@rdivide,Yp,lambda);

%% reset the dimension
resOrder=[1:pmode-1 pmode+1:N];
resYdims=Ydim(resOrder);

%% get opts for pikrpapprox
krp_opts=struct('init',krp_init,'powit',krp_maxit,'constraints',cons,'cons_para',cons_para,'costol',krp_costol);
[P(resOrder)]=pikrpapprox(Yp,resYdims,krp_opts);

%% Final stage : get pmode and permute
P=ktensor(lambda',P);



end